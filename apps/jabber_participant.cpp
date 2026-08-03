/**
 * @file jabber_participant.cpp
 * @brief preCICE participant for coupling acoustic forcing with flow
 * simulations.
 *
 * @details Computing the acoustically-disturbed freestream at every
 * grid point every timestep can be computationally expensive. With that,
 * it can be load-balanced by splitting your MPI communicator such that some
 * ranks are solely just computing the acoustic field constantly without
 * blocking and sending the data to the ranks that need it for boundary
 * conditions or forcing sponge regions. For projects with preCICE, this
 * app can make load-balancing seamless by having jabber be a one-way coupled
 * participant.
 */

#include <jabber/jabber.hpp>
#include <cxxopts.hpp>
#include <precice/precice.hpp>

#ifdef JABBER_WITH_MPI
#include <mpi.h>
#endif // JABBER_WITH_MPI

#include <iostream>

/// Simple macro for enclosing code section to occur only for rank 0
#ifdef JABBER_WITH_MPI
#define ROOT if (rank == 0)
#else
#define ROOT
#endif

using namespace jabber;
using namespace jabber::app;

int main(int argc, char *argv[])
{
   int rank=0, size=1;
#ifdef JABBER_WITH_MPI
   MPI_Init(&argc, &argv);
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   MPI_Comm_size(MPI_COMM_WORLD, &size);

#endif // JABBER_WITH_MPI

   ROOT
   {
      PrintBanner(std::cout);
      std::cout << "Jabber preCICE Participant" << std::endl
                << LINE << std::endl;
   }
   // Option parser:
   cxxopts::Options options("jabber_participant",
                            "preCICE participant for coupling acoustic "
                            "forcing with flow simulations.");
   options.add_options()
   ("c,config", "Config file.", cxxopts::value<std::string>())
   ("h,help", "Print usage information.");
   cxxopts::ParseResult result = options.parse(argc, argv);

   std::string args_str = result.arguments_string();
   ROOT std::cout << "Command Line Arguments:\n\n" << args_str << std::endl
                  << LINE << std::endl;

   if (result.count("help"))
   {

      ROOT std::cout << options.help() << std::endl;
      return 0;
   }
   if (result.count("config") == 0)
   {
      ROOT std::cerr << "Error: no config file specified." << std::endl;
      return 1;
   }

   // Parse config file
   std::string config_file = result["config"].as<std::string>();
   std::ostream *os = nullptr;
   ROOT os = &std::cout;
   TOMLConfigInput conf(config_file, os);
   ROOT std::cout << LINE << std::endl;

   // Get the preCICE input
   const PreciceParams &precice_conf = *(conf.precice_params);

   // Use the dim based on the base flow freestream velocity
   const int dim = conf.base_flow_params.U.size();

   // Initialize preCICE participant
   precice::Participant participant(precice_conf.participant_name,
                                    precice_conf.config_file, rank, size);

   // Set mesh access region [-inf,inf]
   const std::vector<double> mesh_access_region =
      [&]()
   {
      std::vector<double> acc(2*dim);
      for (int i = 0; i < acc.size(); i++)
      {
         acc[i] = (i % 2 == 0)
                  ?  std::numeric_limits<double>::lowest()
                  :  std::numeric_limits<double>::max();
      }
      return acc;
   }
   ();
   participant.setMeshAccessRegion(precice_conf.fluid_mesh_name,
                                   mesh_access_region);
   participant.initialize();

   // Get mesh information from fluid participant
   const int vertex_size =
      participant.getMeshVertexSize(precice_conf.fluid_mesh_name);
   std::vector<double> coords(dim*vertex_size);
   std::vector<int> vertex_ids(vertex_size);
   participant.getMeshVertexIDsAndCoordinates(precice_conf.fluid_mesh_name,
                                              vertex_ids, coords);

   // Set L_z and dz in the oblique wave configs!!!
   // Get all unique zs
   const double tol = 1e-9;
   std::vector<double> zs(vertex_size);
   for (int i = 0; i < vertex_size; i++)
   {
      if (std::find_if(zs.begin(), zs.end(),
         [&](const double &z)
         {
            return std::abs(z-coords[i+2]) <= tol;
         }) == zs.end())
      {
         zs.push_back(coords[i+2]);
      }
   }
   std::sort(zs.begin(), zs.end());

   // Get the length
   const double l_z = zs.back();
   const double dz = zs[1]-zs[0];

   ROOT
   {
      std::cout << "Domain z-length: " << l_z << std::endl;
      std::cout << "dz: " << dz << std::endl;
      std::cout << "Periodic length: " << l_z + dz << std::endl;
      std::cout << LINE << std::endl;
   }

   // Loop through the sources.
   // If source is a RandomPeriodicOblique, set the
   for (Source::ParamsVariant &source : conf.sources_params)
   {
      Source::Params<Source::Option::PSD> *sp = 
         std::get_if<Source::Params<Source::Option::PSD>>(&source);
      if (sp)
      {
         std::visit(
         [&]<DiscMethod::Option O>(DiscMethod::Params<O> &dp)
         {
            if constexpr (O == DiscMethod::Option::RandomPeriodicOblique)
            {
               dp.z_length = l_z;
               dp.dz = dz;
            }
         }
         , sp->disc_params);
      }
   }

   ROOT
   {
      std::cout << "Updated Sources: " << std::endl;
      conf.PrintSourceParams(*os);
      std::cout << std::endl;
   }


#ifdef JABBER_WITH_MPI
   std::span<const double> rank_coords;
   std::span<const int> rank_vertex_ids;

   GetRankPartition<double>(coords, dim, rank, size, rank_coords);
   coords = std::vector<double>(rank_coords.begin(), rank_coords.end());

   GetRankPartition<int>(vertex_ids, 1, rank, size, rank_vertex_ids);
   vertex_ids = std::vector<int>(rank_vertex_ids.begin(),
                                 rank_vertex_ids.end());
#endif // JABBER_WITH_MPI

   // Assemble AcousticField object
   ROOT std::cout << "Assembling acoustic field data... ";
   AcousticField field = InitializeAcousticField(conf, coords, dim);
   ROOT std::cout << "Done!" << std::endl;

   double time = conf.comp_params.t0;
   double dt;

   // Compute acoustic forcing
   while (participant.isCouplingOngoing())
   {
      dt = participant.getMaxTimeStepSize();

      // Compute acoustic forcing
      field.Compute(time);

      // Send data
      participant.writeData(precice_conf.fluid_mesh_name, "rho",
                            vertex_ids, field.Density());
      for (int d = 0; d < dim; d++)
      {
         participant.writeData(precice_conf.fluid_mesh_name,
                               "rhoV" + std::to_string(d+1),
                               vertex_ids, field.Momentum(d));
      }
      participant.writeData(precice_conf.fluid_mesh_name, "rhoE",
                            vertex_ids, field.Energy());
      participant.advance(dt);
      time += dt;
   }

   participant.finalize();

   return 0;
}
