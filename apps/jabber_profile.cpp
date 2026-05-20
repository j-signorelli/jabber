/**
 * @file jabber_profile.cpp
 * @brief Simple profiler tool to obtain execution times of a given config file
 * and grid.
 */

#include <jabber/jabber.hpp>
#include <cxxopts.hpp>

#ifdef JABBER_WITH_MPI
#include <mpi.h>
#endif // JABBER_WITH_MPI

#include <iostream>
#include <cmath>
#include <chrono>
#include <random>

/// Simple macro for enclosing code section to occur only for rank 0
#ifdef JABBER_WITH_MPI
#define ROOT if (rank == 0)
#else
#define ROOT
#endif

using namespace jabber;
using namespace jabber::app;

// Define duration type alias
using dur_t = std::chrono::duration<double, std::micro>;

/// Create a dummy grid.
void CreateGrid(const int dim, const int num_pts_d, const double extent,
                std::vector<double> &coords);

int main(int argc, char *argv[])
{
#ifdef JABBER_WITH_MPI
   MPI_Init(&argc, &argv);

   int rank, size;
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   MPI_Comm_size(MPI_COMM_WORLD, &size);
#endif // JABBER_WITH_MPI

   ROOT
   {
      PrintBanner(std::cout);
      std::cout << "Jabber Profiler Tool" << std::endl
                << LINE << std::endl;
   }
   // Option parser:
   cxxopts::Options options("jabber_profile",
                            "Simple profiler tool to obtain execution times of a given config file \
        and grid.");

   options.add_options()
   ("c,config", "Config file.", cxxopts::value<std::string>())
   ("d,dim", "Grid dimension (1,2,3).",
    cxxopts::value<int>()->default_value("3"))
   ("n,num_points", "Number grid points in each dimension.",
    cxxopts::value<std::size_t>()->default_value("100"))
   ("e,extent", "Grid extent in each direction (such that domain is "
    "[0,extent]^dim).",
    cxxopts::value<double>()->default_value("1.0"))
   ("p,passes", "Number of passes to Compute() to profile, using "
    "randomized times.",
    cxxopts::value<int>()->default_value("10000"))
   ("w,warmup", "Number of warmup passes to Compute(), using randomized "
    "times.", cxxopts::value<int>()->default_value("1000"))
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

   const int dim = result["dim"].as<int>();
   const std::size_t num_pts_d = result["num_points"].as<std::size_t>();
   const std::size_t num_pts_total = std::pow(num_pts_d,dim);
   const double extent = result["extent"].as<double>();
   const int passes = result["passes"].as<int>();
   const int warmup_passes = result["warmup"].as<int>();

   // Parse config file
   std::string config_file = result["config"].as<std::string>();
   std::ostream *os = nullptr;
   ROOT os = &std::cout;
   TOMLConfigInput conf(config_file, os);
   ROOT std::cout << LINE << std::endl;

   // Output grid information to console
   ROOT
   {
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
   }

   // Create a simple [0,1]^dim grid
   std::vector<double> coords(num_pts_total*dim);
   CreateGrid(dim, num_pts_d, extent, coords);

#ifdef JABBER_WITH_MPI
   // Partition the coordinates.
   std::span<const double> rank_coords;
   GetRankPartition<double>(coords, dim, rank, size, rank_coords);
   coords = std::vector<double>(rank_coords.begin(), rank_coords.end());
#endif // JABBER_WITH_MPI

   // Initialize AcousticField
   AcousticField field = InitializeAcousticField(conf, coords, dim);

   // Create an array of randomized times
   std::mt19937 gen(0);
   std::uniform_real_distribution<double> real_dist(0,1);
   std::vector<double> time_rand(passes+warmup_passes);
   for (int i = 0; i < time_rand.size(); i++)
   {
      time_rand[i] = real_dist(gen);
   }

   // When JABBER_WITH_MPI, this holds max compute time for each iteration,
   // across all ranks.
   std::vector<dur_t> compute_times;
   ROOT compute_times.resize(passes+warmup_passes);

#ifdef JABBER_WITH_MPI
   std::vector<dur_t> local_compute_times(passes+warmup_passes);

   // Rank of max compute time for each iteration - only store on root
   std::vector<int> max_compute_ranks;
   ROOT max_compute_ranks.resize(passes+warmup_passes);

   ROOT std::cout << std::endl
                  << "Running with MPI enabled! Outputted times are the max "
                  "across all ranks." << std::endl << std::endl;
   struct TimeRank
   {
      double t;
      int r;
   } time_rank;
#endif // JABBER_WITH_MPI

   // Run profiling loop
   for (int i = 0; i < warmup_passes+passes; i++)
   {
      const std::chrono::time_point<std::chrono::steady_clock> start =
         std::chrono::steady_clock::now();
      field.Compute(time_rand[i]);
      const std::chrono::time_point<std::chrono::steady_clock> end =
         std::chrono::steady_clock::now();

#ifdef JABBER_WITH_MPI
      local_compute_times[i] = end - start;

      time_rank.t = local_compute_times[i].count();
      time_rank.r = rank;

      // Get the max compute time and rank across all ranks
      MPI_Allreduce(MPI_IN_PLACE, &time_rank, 1, MPI_DOUBLE_INT, MPI_MAXLOC,
                    MPI_COMM_WORLD);

      ROOT
      {
         compute_times[i] = dur_t(time_rank.t);
         max_compute_ranks[i] = time_rank.r;
      }
#else
      compute_times[i] = end - start;
#endif // JABBER_WITH_MPI

      ROOT
      {
         if (i < warmup_passes)
         {
            std::cout << "Warmup Pass #" << (i+1) << ": " << compute_times[i]
                      << std::endl;
         }
         else
         {
            std::cout << "Pass #" << (i+1-warmup_passes) << ": "
                      << compute_times[i] << std::endl;;
         }
      }
   }

   // Compute the average compute time after warmups
   ROOT
   {
      dur_t total_dur = std::accumulate(std::next(compute_times.begin(),
                                                  warmup_passes),
                                        compute_times.end(),
                                        dur_t(0.0));
      dur_t ave_dur = total_dur/passes;

      std::cout << LINE << std::endl
                << "Average Compute() Time: " << ave_dur << std::endl;
   }

#ifdef JABBER_WITH_MPI

   // Have each rank compute their local average duration
   dur_t local_total_dur =
      std::accumulate(std::next(local_compute_times.begin(), warmup_passes),
                      local_compute_times.end(),
                      dur_t(0.0));
   dur_t local_ave_dur = local_total_dur/passes;

   // Gather all mean durations + num coordinates to root
   struct RankData
   {
      double local_ave;
      int num_pts;
   };

   RankData rank_data;
   rank_data.local_ave = local_ave_dur.count();
   rank_data.num_pts = coords.size()/dim;

   std::vector<RankData> all_rank_data;
   ROOT all_rank_data.resize(size);
   MPI_Gather(&rank_data, 1, MPI_DOUBLE_INT, all_rank_data.data(), 1,
              MPI_DOUBLE_INT, 0, MPI_COMM_WORLD);

   ROOT
   {
      // Compute histogram of worst-rank counts to determine worst.
      std::vector<std::size_t> rank_counts(size, 0.0);
      for (std::size_t i = 0; i < passes; i++)
      {
         rank_counts[max_compute_ranks[warmup_passes+i]]++;
      }
      const int worst_rank =
      std::distance(rank_counts.begin(),
                    std::max_element(rank_counts.begin(),rank_counts.end()));

      std::cout << "Rank " << worst_rank << " was the slowest rank for " <<
                rank_counts[worst_rank] << " iterations."
                << std::endl << std::endl;

      // Sort indices of rank_ave_dur based on their value
      std::vector<int> r_idxs(size);
      std::iota(r_idxs.begin(), r_idxs.end(), 0);
      std::sort(r_idxs.begin(), r_idxs.end(),
                [&all_rank_data](std::size_t i1, std::size_t i2)
      {
         return all_rank_data[i1].local_ave >
         all_rank_data[i2].local_ave;
      });

      std::cout << "Printing rank mean compute times in order of "
                "slowest to fastest..." << std::endl;

      for (const int &r : r_idxs)
      {
         std::cout << std::format("Rank {:>3}:   {:<10} with {} points",
                                  r, dur_t(all_rank_data[r].local_ave),
                                  all_rank_data[r].num_pts)
                   << std::endl;
      }
   }
#endif // JABBER_WITH_MPI
   return 0;
}

void CreateGrid(const int dim, const int num_pts_d, const double extent,
                std::vector<double> &coords)
{
   const double h = extent/(num_pts_d-1);
   int num_pts_total = 1;
   for (int d = 0; d < dim; d++)
   {
      num_pts_total *= num_pts_d;
   }
   coords.resize(num_pts_total*dim);

   if (dim == 1)
   {
      for (std::size_t i = 0; i < num_pts_d; i++)
      {
         coords[i] = h*i;
      }
   }
   else if (dim == 2)
   {
      for (std::size_t i = 0; i < num_pts_d; i++)
      {
         const double x = h*i;
         for (std::size_t j = 0; j < num_pts_d; j++)
         {
            const std::size_t completed_pts = i*num_pts_d + j;
            coords[completed_pts*dim] = x;
            coords[completed_pts*dim+1] = h*j;
         }
      }
   }
   else // dim == 3
   {
      for (std::size_t i = 0; i < num_pts_d; i++)
      {
         const double x = h*i;
         for (std::size_t j = 0; j < num_pts_d; j++)
         {
            const double y = h*j;
            for (std::size_t k = 0; k < num_pts_d; k++)
            {
               const std::size_t completed_pts = i*num_pts_d*num_pts_d
                                                 + j*num_pts_d + k;
               coords[completed_pts*dim] = x;
               coords[completed_pts*dim+1] = y;
               coords[completed_pts*dim+2] = h*k;
            }
         }
      }
   }
}
