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

void WriteWaves(const std::span<const Wave> &waves, std::ostream &out)
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

      for (; range_it != field_view.end(); field_it = *(++range_it))
      {
         const double val = std::stod(std::string(field_it.begin(),
                                                  field_it.end()));
         w.k_hat.push_back(val);
      }
      waves.emplace_back(w);
   }
}

double ComputeWavenumber(const std::span<const double> &U_infty,
                         const double &c_infty,
                         const double &wave_freq,
                         const std::span<const double> &wave_k_hat,
                         const char &wave_speed)
{
   // Compute denom = U·k_hat±c
   double denom = (wave_speed == 'S' ? -c_infty : c_infty);
   for (int d = 0; d < U_infty.size(); d++)
   {
      denom += U_infty[d]*wave_k_hat[d];
   }

   if (denom <= 0.0)
   {
      throw std::invalid_argument(
         "Invalid wave orientation for given freestream.");
   }

   return 2*M_PI*wave_freq/denom;
}


AcousticField::AcousticField(int dim, const std::span<const double> &coords,
                             const BaseFlow &base_flow, Kernel kernel)
   : dim_(dim),
     base_flow_(base_flow),
     kernel_(kernel),
     num_pts_(coords.size()/dim_),
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
   kernel_args_.k_dot_x_p_psi.resize(NumWaves()*NumPoints());

   // Note that performance of below was not carefully considered
   for (int w = 0; w < NumWaves(); w++)
   {
      const Wave &wave = Waves()[w];

      kernel_args_.rho_coeffs[w] = wave.amplitude/(base_flow_.c*base_flow_.c);
      kernel_args_.rhoE_coeffs[w] = wave.amplitude/(base_flow_.gamma - 1.0);
      kernel_args_.wave_omegas[w] = 2*M_PI*wave.frequency;

      // Compute rhoV_coeffs
      const int speed_encoder = (wave.speed == 'S' ? -1 : 1);
      for (int d = 0; d < Dim(); d++)
      {
         kernel_args_.rhoV_coeffs[d*NumWaves() + w] =
            speed_encoder*wave.k_hat[d]*wave.amplitude
            /(base_flow_.rho*base_flow_.c);
      }

      // Compute magnitude of wavelength vector k
      const double k = ComputeWavenumber(base_flow_.U, base_flow_.c, wave);

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
         }
         ();
         kernel_args_.k_dot_x_p_psi[idx] = wave.phase;
         for (int d = 0; d < Dim(); d++)
         {
            kernel_args_.k_dot_x_p_psi[idx] +=
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
   [&]<std::size_t... Dims>(const std::index_sequence<Dims...> &)
   {
      ([&]()
      {
         if (Dim() == Dims)
         {
            if (kernel_ == Kernel::GridPoint)
            {
               ComputeKernel<Dims, true>(NumPoints(), base_flow_.rho,
                                         base_flow_.p,
                                         base_flow_.U.data(), base_flow_.gamma,
                                         NumWaves(),
                                         t, kernel_args_.rho_coeffs.data(),
                                         kernel_args_.rhoV_coeffs.data(),
                                         kernel_args_.rhoE_coeffs.data(),
                                         kernel_args_.wave_omegas.data(),
                                         kernel_args_.k_dot_x_p_psi.data(),
                                         rho_.data(), rhoV_.data(),
                                         rhoE_.data());
            }
            else if (kernel_ == Kernel::Wave)
            {
               ComputeKernel<Dims, false>(NumPoints(), base_flow_.rho,
                                          base_flow_.p,
                                          base_flow_.U.data(), base_flow_.gamma,
                                          NumWaves(),
                                          t, kernel_args_.rho_coeffs.data(),
                                          kernel_args_.rhoV_coeffs.data(),
                                          kernel_args_.rhoE_coeffs.data(),
                                          kernel_args_.wave_omegas.data(),
                                          kernel_args_.k_dot_x_p_psi.data(),
                                          rho_.data(), rhoV_.data(),
                                          rhoE_.data());
            }
            else
            {
               throw std::logic_error("Unimplemented kernel type!");
            }
         }
      }
      (), ...);
   }
   (std::index_sequence<1,2,3> {});
}



} // namespace jabber
