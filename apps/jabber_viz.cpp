/**
 * @file jabber_viz.cpp
 * @brief Visualizer of acoustic field from config file using MFEM w/ GLVis.
 */

#include <jabber/jabber.hpp>
#include <cxxopts.hpp>
#include <mfem.hpp>

#include <iostream>
#include <algorithm>

using namespace jabber;
using namespace jabber::app;
using namespace mfem;

/**
 * @brief Simple field visualizer, copied directly from 
 * mfem/miniapps/common/fem_extras.cpp/.hpp
 * 
 */
void VisualizeField(socketstream &sock, const char *vishost, int visport,
                    GridFunction &gf, const char *title,
                    int x = 0, int y = 0, int w = 400, int h = 400,
                    const char *keys = NULL, bool vec = false);

int main(int argc, char *argv[])
{
   PrintBanner(std::cout);
   std::cout << "Jabber Acoustic Field Visualizer" << std::endl 
               << LINE << std::endl;

    // Option parser:
   cxxopts::Options options("jabber_viz", 
      "Visualizer of acoustic field from config file using MFEM w/ GLVis.");

   options.add_options()
      ("c,config", "Config file.", cxxopts::value<std::string>())
      ("d,dim", "Grid dimension (1,2,3).", 
                        cxxopts::value<int>()->default_value("2"))
      ("n,num-points", "Number grid points in each dimension.",
                        cxxopts::value<std::size_t>()->default_value("100"))
      ("e,extent", "Grid extent in each direction (such that domain is "
                     "[0,extent]^dim).",
                        cxxopts::value<double>()->default_value("1.0"))
      ("f,fields", "Fields to visualize with GLVis ('rho', 'rhoV', 'rhoE').",
                        cxxopts::value<std::vector<std::string>>()
                        ->default_value("rho,rhoV,rhoE"))
      ("s,dt", "Timestep to use.", 
                        cxxopts::value<double>()->default_value("0.0"))
      ("t,timesteps", "Number of timesteps to run to.",
                        cxxopts::value<int>()->default_value("500"))
      ("h,help", "Print usage information.");

   cxxopts::ParseResult result = options.parse(argc, argv);
   
   std::string args_str = result.arguments_string();
   std::cout << "Command Line Arguments:\n\n" << args_str << std::endl 
               << LINE << std::endl;

   if (result.count("help"))
   {
      std::cout << options.help() << std::endl;
      return 0;
   }
   if (result.count("config") == 0)
   {
      std::cerr << "Error: no config file specified." << std::endl;
      return 1;
   }

   const int dim = result["dim"].as<int>();
   const std::size_t num_pts_d = result["num_points"].as<std::size_t>();
   const std::size_t num_pts_total = std::pow(num_pts_d,dim);
   const double extent = result["extent"].as<double>();
   const std::vector<std::string> fields = result["fields"]
                                             .as<std::vector<std::string>>();
   const double timestep = result["dt"].as<double>();
   const int num_timesteps = result["timesteps"].as<int>();

   // Parse config file
   std::string config_file = result["config"].as<std::string>();
   TOMLConfigInput conf(config_file, &std::cout);
   std::cout << LINE << std::endl;

   // Output grid information to console
   std::cout << "Grid\n";
   std::cout << "\tDimension: ";
   for (int d = 0; d < dim; d++)
   {
      std::cout << num_pts_d << (d+1==dim ? "" : "x");
   }
   std::cout << std::endl;
   std::cout << "\tExtents: ";
   for (int d = 0; d < dim; d++)
   {
      std::cout << "[0," << extent << "]" << ((d + 1 == dim) ? "" : "x");
   }
   std::cout << std::endl;
   std::cout << "\tNumber of points: " << num_pts_total << std::endl;
   std::cout << "\tSpacing: " << (extent/(num_pts_d-1.0)) << std::endl;

   // Create MFEM mesh
   Mesh mesh;
   if (dim == 1)
   {
      mesh = Mesh::MakeCartesian1D(num_pts_d, extent);
   }
   else if (dim == 2)
   {
      mesh = Mesh::MakeCartesian2D(num_pts_d, num_pts_d, 
                                    Element::Type::QUADRILATERAL, true,
                                    extent, extent);
   }
   else
   {
      mesh = Mesh::MakeCartesian3D(num_pts_d, num_pts_d, num_pts_d,
                                    Element::Type::HEXAHEDRON, extent,
                                    extent, extent);
   }

   // Get mesh coordinates ordered as XYZ XYZ
   H1_FECollection fec(1, dim);
   FiniteElementSpace coords_fespace(&mesh, &fec, dim, Ordering::byVDIM);
   GridFunction coords_gf(&coords_fespace);
   mesh.GetNodes(coords_gf);

   // Copy coordinates into std::vector
   std::vector<double> coords(coords_gf.Size());
   for (std::size_t i = 0; i < coords_gf.Size(); i++)
   {
      coords[i] = coords_gf[i];
   }
   coords_gf.Vector::Destroy(); // clear memory

   // Initialize AcousticField
   AcousticField field = InitializeAcousticField(conf, coords, dim);

   // Initialize solution FE spaces + grid functions
   FiniteElementSpace s_fespace(&mesh, &fec, 1);
   FiniteElementSpace v_fespace(&mesh, &fec, 1, Ordering::byNODES);
   GridFunction rho_gf(&s_fespace, field.Density().data()),
               rhoV_gf(&v_fespace, field.Momentum().data()),
               rhoE_gf(&s_fespace, field.Energy().data());

   const int visport = 19916;
   constexpr std::string_view vishost = "localhost";
   std::string keys;
   if (dim == 2)
   {
      keys = "lARjc//";
   }
   else
   {
      keys ="lAc///[[";
   }
   const int Wx = 500, Wy = 500;
   socketstream rho_sock;
   socketstream rhoV_sock;
   socketstream rhoE_sock;


   double t = conf.Comp().t0;
   for (int i = 0; i < num_timesteps; i++)
   {

      field.Compute(t);

      t += timestep;

      // Write to GLVis windows
      int offset = 0;
      for (const std::string &field : fields)
      {
         if (field == "rho")
         {
            VisualizeField(rho_sock, vishost.data(), visport, rho_gf, "Density",
                           offset, 0, Wx, Wy, keys.data());
         }
         else if (field == "rhoV")
         {
            VisualizeField(rhoV_sock, vishost.data(), visport, rhoV_gf, "Momentum",
                           offset, 0, Wx, Wy, keys.data(), true);
         }
         else if (field == "rhoE")
         {
            VisualizeField(rhoE_sock, vishost.data(), visport, rhoE_gf, "Energy",
                           offset, 0, Wx, Wy, keys.data());
         }
         offset += Wx;
      }
   }
   return 0;
}

void VisualizeField(socketstream &sock, const char *vishost, int visport,
                    GridFunction &gf, const char *title,
                    int x, int y, int w, int h, const char * keys, bool vec)
{
   Mesh &mesh = *gf.FESpace()->GetMesh();

   bool newly_opened = false;
   int connection_failed;

   do
   {
      if (!sock.is_open() || !sock)
      {
         sock.open(vishost, visport);
         sock.precision(8);
         newly_opened = true;
      }
      sock << "solution\n";

      mesh.Print(sock);
      gf.Save(sock);

      if (newly_opened)
      {
         sock << "window_title '" << title << "'\n"
              << "window_geometry "
              << x << " " << y << " " << w << " " << h << "\n";
         if ( keys ) { sock << "keys " << keys << "\n"; }
         else { sock << "keys maaAc\n"; }
         if ( vec ) { sock << "vvv"; }
         sock << std::endl;
      }

      connection_failed = !sock && !newly_opened;
   }
   while (connection_failed);
}
