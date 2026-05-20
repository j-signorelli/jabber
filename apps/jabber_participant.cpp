/**
 * @file jabber_participant.cpp
 * @brief preCICE participant for coupling acoustic forcing with flow
 * simulations.
 *
 * @details
 *
 * @todo EMPHASIZE IMPORTANCE OF THIS -- computing acoustic forcing every
 * timestep by fluid solver is super expensive
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
      ROOT std::cout << "Error: no config file specified." << std::endl;
      return 1;
   }

   // Parse config file
   std::string config_file = result["config"].as<std::string>();
   std::ostream *os = nullptr;
   ROOT os = &std::cout;
   TOMLConfigInput conf(config_file, os);
   ROOT std::cout << LINE << std::endl;

   // Get the preCICE input
   const PreciceParams &precice_conf = *(conf.Precice());

   // Initialize preCICE participant
   precice::Participant participant(precice_conf.participant_name,
                                    precice_conf.config_file, rank, size);
   participant.setMeshAccessRegion(precice_conf.fluid_mesh_name,
                                   precice_conf.mesh_access_region);
   participant.initialize();

   // Get mesh information from fluid participant
   int dim = participant.getMeshDimensions(precice_conf.fluid_mesh_name);
   int vertex_size = participant.getMeshVertexSize(
                        precice_conf.fluid_mesh_name);
   std::vector<double> coords(dim*vertex_size);
   std::vector<int> vertex_ids(vertex_size);
   participant.getMeshVertexIDsAndCoordinates(precice_conf.fluid_mesh_name,
                                              vertex_ids, coords);

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

   double time = conf.Comp().t0;
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
