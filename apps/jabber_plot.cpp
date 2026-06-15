/**
 * @file jabber_plot.cpp
 * @brief Generate a plot of the final wave spectra from a config file.
 */

#include <jabber/jabber.hpp>
#include <cxxopts.hpp>

#include <iostream>
#include <algorithm>
#include <fstream>

using namespace jabber;
using namespace jabber::app;

int main(int argc, char *argv[])
{
   PrintBanner(std::cout);
   std::cout << "Jabber Spectra Plotter" << std::endl
             << LINE << std::endl;

   // Option parser:
   cxxopts::Options options("jabber_plot",
                            "Generate a plot of the final wave spectra from "
                            " a config file.");

   options.add_options()
   ("c,config", "Config file.", cxxopts::value<std::string>())
   ("l,log", "Plot on a log-log scale.",
    cxxopts::value<bool>()->default_value("false"))
   ("h,help", "Print usage information.")
   ("n,nondim",
    "Nondimensionalize the pressure wave amplitudes using the "
    "input base flow pressure.",
    cxxopts::value<bool>()->default_value("false"));

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

   const bool loglog = result["log"].as<bool>();
   const bool nd = result["nondim"].as<bool>();

   // Parse config file
   std::string config_file = result["config"].as<std::string>();
   TOMLConfigInput conf(config_file, &std::cout);
   std::cout << LINE << std::endl;

   std::vector<Wave> waves;
   for (const Source::ParamsVariant &spv : conf.Sources())
   {
      std::visit(SourceVisitor{conf.BaseFlow(), waves}, spv);
   }

   std::vector<double> freqs(waves.size()), amps(waves.size());
   for (std::size_t i = 0; i < waves.size(); i++)
   {
      freqs[i] = waves[i].frequency;
      amps[i] = waves[i].amplitude / (nd ? conf.BaseFlow().p : 1.0);
   }

   std::FILE *gnuplot = popen("gnuplot", "w");

   if (loglog)
   {
      std::fprintf(gnuplot, "set logscale xy\n");
   }

   std::fprintf(gnuplot, "unset key\n");

std:
   fprintf(gnuplot, "set xlabel 'Frequency'\n");
   std::fprintf(gnuplot, "set ylabel 'Wave Amplitude'\n");
   std::fprintf(gnuplot, "plot '-' with points pt 5\n");
   for (std::size_t i = 0; i < freqs.size(); i++)
   {
      std::fprintf(gnuplot, "%s", std::format("{} {}\n",
                                              freqs[i], amps[i]).c_str());
   }
   std::fprintf(gnuplot, "e\n");
   std::cout << "Enter to close plot...";
   std::fflush(gnuplot);
   std::cin.get();
   pclose(gnuplot);

   return 0;
}
