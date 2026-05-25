/**
 * @file jabber_psd.cpp
 * @brief Compute and plot a PSD from a probe of the exact
 * flowfield computed by Jabber using Welch's method.
 * @details Welch's method here has been verified to match SciPy's
 * implementation of Welch's method.
 */

#include <jabber/jabber.hpp>
#include <cxxopts.hpp>
#include <pocketfft_hdronly.h>

#include <iostream>
#include <cmath>
#include <fstream>

#include <fmt/core.h>
#include <fmt/chrono.h>

using namespace jabber;
using namespace jabber::app;

double HammingWindow(std::size_t N, std::size_t n);

int main(int argc, char *argv[])
{
   PrintBanner(std::cout);
   std::cout << "Jabber PSD" << std::endl
               << LINE << std::endl;

   // Option parser:
   cxxopts::Options options("jabber_psd", 
      "Compute a PSD from a probe of the exact flowfield computed "
      "by Jabber using Welch's method.");

   options.add_options()
      ("c,config", "Config file.", cxxopts::value<std::string>())
      ("d,dt", "Timestep to use.", 
         cxxopts::value<double>()->default_value("3.72961861118742e-7"))
      ("t,timesteps", "Number of timesteps to run to.", 
         cxxopts::value<std::size_t>()->default_value("1000000"))
      ("s,nperseg", "Number of points in each segment.",
         cxxopts::value<std::size_t>()->default_value("256"))
      ("o,noverlap", "Number of point overlap in segments. "
                     " Defaults to nperseg/2 (50%)", 
         cxxopts::value<std::size_t>())
      ("w,write-psd-file", "Filename to write PSD data to (if included) as a CSV.",
         cxxopts::value<std::string>())
      ("r,write-press-file", "Filename to write pressure perturbation data to "
                             "(if included).",
         cxxopts::value<std::string>())
      ("p,plot", "Generate a plot of the computed PSD data.")
      ("l,log", "Plot on a log-log scale.",
         cxxopts::value<bool>()->default_value("false"))
      ("i,input-psd", "Input PSD CSV file to plot computed PSD against",
         cxxopts::value<std::string>())
      ("n,nondim", 
         "Nondimensionalize the PSD using the input base flow pressure.",
         cxxopts::value<bool>()->default_value("false"))
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


   const double dt = result["dt"].as<double>();
   const std::size_t nt = result["timesteps"].as<std::size_t>();
   const std::size_t nperseg = result["nperseg"].as<std::size_t>();
   const std::size_t noverlap = (result.count("noverlap") == 0) ? nperseg/2
                                    : result["noverlap"].as<std::size_t>();
   const bool nd = result["nondim"].as<bool>();
   
   if (nperseg > nt)
   {
      throw std::invalid_argument("nperseg must be less than nt!");
      return 1;
   }

   // Parse config file
   std::string config_file = result["config"].as<std::string>();
   TOMLConfigInput conf(config_file, &std::cout);
   std::cout << LINE << std::endl;
   std::vector<double> coords = {0.1, 0.1, 0.1};

   // Initialize AcousticField
   AcousticField field = InitializeAcousticField(conf, coords, 3);

   double time = 0.0;
   const double c_sq = conf.BaseFlow().gamma*
                           conf.BaseFlow().p/conf.BaseFlow().rho;

   // Initialize vector of all pressures.
   std::vector<double> p_prime(nt, 0.0);

   // Run through all times, storing pressures at each timestep
   for (int i = 0; i < nt; i++)
   {
      field.Compute(time);
      time += dt;
      p_prime[i] = c_sq*(field.Density()[0] - conf.BaseFlow().rho);
      if (nd)
      {
         p_prime[i] /= conf.BaseFlow().p;
      }
   }

   const std::size_t shift = nperseg - noverlap;
   const std::size_t num_segs = 1+(nt-nperseg)/shift;

   // Compute W factor for modified periodogram + assemble window fxn
   std::vector<double> window(nperseg);
   double sum_W_sq = 0.0;
   for (std::size_t w = 0; w < nperseg; w++)
   {
      const double w_val = HammingWindow(nperseg, w);
      window[w] = w_val;
      sum_W_sq += w_val*w_val;
   }

   // Compute the DFT for each segment k
   const double fs = 1.0/dt;
   std::vector<double> xw(nperseg);
   std::vector<std::complex<double>> dft(nperseg/2+1);
   std::vector<double> psd(nperseg/2+1, 0.0);                
   for (std::size_t k = 0; k < num_segs; k++)
   {
      // Compute pressure * window
      for (std::size_t w = 0; w < nperseg; w++)
      {
         xw[w] = p_prime[k*shift+w]*window[w];
      }

      // Compute DFT of each segment
      pocketfft::r2c({nperseg}, {sizeof(double)}, 
                     {sizeof(std::complex<double>)},
                     0, pocketfft::FORWARD, xw.data(), 
                     dft.data(), 1.0);

      // Compute modified periodogram of DFT + add averaged contribution 
      // of this segment k to PSD
      for (std::size_t i = 0; i < nperseg/2+1; i++)
      {
         const double mod_p_i = std::norm(dft[i])/(sum_W_sq*fs);
         psd[i] += mod_p_i/num_segs;
      }
   }

   // Apply one-sided correction
   const std::size_t up_to = (nperseg % 2 == 1) ? nperseg/2+1 : nperseg/2;
   for (std::size_t i = 1; i < up_to; i++)
   {
      psd[i] *= 2;
   }

   // Get the frequencies
   std::vector<double> freqs(psd.size());
   for (std::size_t i = 0; i < nperseg/2+1; i++)
   {
      freqs[i] = i*fs/nperseg;
   }
   
   // Write PSD to file:
   if (result.count("write-psd-file") > 0)
   {
      const std::string file_name = result["write-psd-file"].as<std::string>();
      std::ofstream psd_file(file_name);
      for (std::size_t i = 0; i < nperseg/2+1; i++)
      {
         psd_file << fmt::format("{},{}\n", freqs[i], psd[i]);
      }
   }

   // Write pressure perturbations to file
   if (result.count("write-press-file") > 0)
   {
      const std::string file_name = result["write-press-file"]
                                                         .as<std::string>();
      std::ofstream press_file(file_name);
      for (std::size_t i = 0; i < p_prime.size(); i++)
      {
         press_file << fmt::format("{}\n", p_prime[i]);
      }
   }

   // Plot
   if (result.count("plot") > 0)
   {
      std::FILE* gnuplot = popen("gnuplot", "w");

      if (result.count("log") > 0)
      {
         std::fprintf(gnuplot, "set logscale xy\n");

         // Remove the DC component altogether
         freqs.erase(freqs.begin());
         psd.erase(psd.begin());
      }

      
      std:fprintf(gnuplot, "set xlabel 'Frequency'\n");
      std::fprintf(gnuplot, "set ylabel 'PSD'\n");
      std::fprintf(gnuplot, "plot '-' title 'Computed' with points pt 1");
      if (result.count("input-psd") > 0)
      {
         std::fprintf(gnuplot, ", '-' title 'Input PSD' with line");              
      }
      std::fprintf(gnuplot, "\n");
      for (std::size_t i = 0; i < freqs.size(); i++)
      {
         std::fprintf(gnuplot, "%s", fmt::format("{} {}\n", 
                                          freqs[i], psd[i]).c_str());
      }
      std::fprintf(gnuplot, "e\n");

      // Include source PSD if provided
      if (result.count("input-psd") > 0)
      {
         using enum InputXY::Option;
         InputXY::Params<FromCSV> input;
         input.file = result["input-psd"].as<std::string>();
         std::vector<double> in_freqs, in_psd;
         InputXYVisitor psd_in{in_freqs, in_psd};
         psd_in(input);
         for (std::size_t i = 0; i < in_freqs.size(); i++)
         {
            std::fprintf(gnuplot, "%s", fmt::format("{} {}\n", 
                                             in_freqs[i], in_psd[i]).c_str());
         }
         std::fprintf(gnuplot, "e\n");
         
      }

      std::cout << "Enter to close plot...";
      std::fflush(gnuplot);
      std::cin.get();
      pclose(gnuplot);
   }

   return 0;
}



double HammingWindow(std::size_t N, std::size_t n)
{
   constexpr double a0 = 0.54;
   return a0 - (1-a0)*std::cos(2*M_PI*n/(N-1));
}