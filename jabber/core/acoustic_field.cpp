#include "acoustic_field.hpp"
#include "kernels.hpp"

#include <math.h>
#include <numeric>
#include <iostream>
#include <format>
#include <string>
#include <ranges>

namespace jabber
{

void WriteWaves(std::span<const Wave> waves, std::ostream &out)
{
   for (std::size_t i = 0; i < waves.size(); i++)
   {
      out << std::format("{},{},{},{}", waves[i].amplitude,
                                          waves[i].frequency,
                                          waves[i].phase, 
                                          waves[i].speed);
      for (std::size_t d = 0; d < waves[i].k_hat.size(); d++)
      {
         out << std::format(",{}", waves[i].k_hat[d]);
         out << (d+1==waves[i].k_hat.size() ? "\n" : "");
      }
   }
}

void ReadWaves(std::istream &in, std::vector<Wave> &waves)
{
   for (std::string line; std::getline(in, line);)
   {
      auto field_view = std::ranges::views::split(line, 
                                                   std::string_view(","));
      auto range_it = field_view.begin();
      auto field_it = *range_it;

      Wave w;
      w.amplitude = std::stod(std::string(field_it.begin(), field_it.end()));
      field_it = *(++range_it);

      w.frequency = std::stod(std::string(field_it.begin(), field_it.end()));
      field_it = *(++range_it);

      w.phase = std::stod(std::string(field_it.begin(), field_it.end()));
      field_it = *(++range_it);

      w.speed = std::string(field_it.begin(), field_it.end())[0];
      field_it = *(++range_it);
      
      for (;range_it != field_view.end(); field_it = *(++range_it))
      {
         const double val = std::stod(std::string(field_it.begin(), 
                                                   field_it.end()));
         w.k_hat.push_back(val);
      }
      waves.emplace_back(w);
   }
}

AcousticField::AcousticField(int dim, std::span<const double> coords,
                  double p_bar, double rho_bar,
                  const std::vector<double> U_bar, double gamma,
                  Kernel kernel)
: kernel_(kernel),
  dim_(dim),
  num_pts_(coords.size()/dim_),
  p_bar_(p_bar),
  rho_bar_(rho_bar),
  U_bar_(U_bar),
  gamma_(gamma),
  c_bar_(std::sqrt(gamma_*p_bar_/rho_bar_)),
  coords_(dim_)
{
   
   // Store the coordinates in an SoA-style
   for (int d = 0; d < Dim(); d++)
   {
      coords_[d].resize(NumPoints());
      for (int i = 0; i < NumPoints(); i++)
      {
         coords_[d][i] = coords[i*Dim() + d];
      }
   }
}

void AcousticField::Finalize()
{
   // Allocate non-time-varying constants
   kernel_args_.rho_coeffs.resize(NumWaves());
   kernel_args_.rhoV_coeffs.resize(Dim()*NumWaves());
   kernel_args_.rhoE_coeffs.resize(NumWaves());
   kernel_args_.wave_omegas.resize(NumWaves());
   kernel_args_.k_dot_x_p_phi.resize(NumWaves()*NumPoints());

   // Note that performance of below was not carefully considered
   for (int w = 0; w < NumWaves(); w++)
   {
      const Wave &wave = Waves()[w];

      kernel_args_.rho_coeffs[w] = wave.amplitude/(c_bar_*c_bar_);
      kernel_args_.rhoE_coeffs[w] = wave.amplitude/(gamma_ - 1.0);
      kernel_args_.wave_omegas[w] = 2*M_PI*wave.frequency;

      // Compute denom = U·k_hat±c and set rhoV_coeffs
      double denom = (wave.speed == 'S' ? -c_bar_ : c_bar_);
      const int speed_encoder = (wave.speed == 'S' ? -1 : 1);
      for (int d = 0; d < Dim(); d++)
      {
         denom += U_bar_[d]*wave.k_hat[d];
         kernel_args_.rhoV_coeffs[d*NumWaves() + w] = 
            speed_encoder*wave.k_hat[d]*wave.amplitude/(rho_bar_*c_bar_);
      }

      // Compute magnitude of wavelength vector k
      const double k = kernel_args_.wave_omegas[w]/denom;

      // Compute + set k·x+φ
      for (std::size_t i = 0; i < NumPoints(); i++)
      {  
         const std::size_t idx = 
         [&]()
         {
            switch (kernel_)
            {
               case Kernel::GridPoint:
                  return w*NumPoints() + i;
                  break;
               case Kernel::Wave:
                  return i*NumWaves() + w;
                  break;
               default:
                  throw std::logic_error("Unimplemented kernel type!");
            }
         }();
         kernel_args_.k_dot_x_p_phi[idx] = wave.phase;
         for (int d = 0; d < Dim(); d++)
         {
            kernel_args_.k_dot_x_p_phi[idx] += 
                                    wave.k_hat[d]*k*coords_[d][i];
         }
      }
   }

   // Allocate flow solution memory
   rho_.resize(NumPoints());
   rhoV_.resize(NumPoints()*Dim());
   rhoE_.resize(NumPoints());
}

void AcousticField::Compute(double t)
{
   // Dispatch to appropriate kernel
   [&]<std::size_t... Dims>(const std::index_sequence<Dims...>&)
   {
      ([&]()
       {
         if (Dim() == Dims)
         {
            if (kernel_ == Kernel::GridPoint)
            {
               ComputeKernel<Dims, true>(NumPoints(), rho_bar_, p_bar_, 
                                    U_bar_.data(), gamma_, NumWaves(), t,
                                    kernel_args_.rho_coeffs.data(),
                                    kernel_args_.rhoV_coeffs.data(),
                                    kernel_args_.rhoE_coeffs.data(), 
                                    kernel_args_.wave_omegas.data(), 
                                    kernel_args_.k_dot_x_p_phi.data(), 
                                    rho_.data(), rhoV_.data(), rhoE_.data());
            }
            else if (kernel_ == Kernel::Wave)
            {
               ComputeKernel<Dims, false>(NumPoints(), rho_bar_, p_bar_, 
                                    U_bar_.data(), gamma_, NumWaves(), t,
                                    kernel_args_.rho_coeffs.data(),
                                    kernel_args_.rhoV_coeffs.data(),
                                    kernel_args_.rhoE_coeffs.data(), 
                                    kernel_args_.wave_omegas.data(), 
                                    kernel_args_.k_dot_x_p_phi.data(), 
                                    rho_.data(), rhoV_.data(), rhoE_.data());
            }
            else
            {
               throw std::logic_error("Unimplemented kernel type!");
            }
         }
       }(), ...);
   }(std::index_sequence<1,2,3>{});
}



} // namespace jabber
