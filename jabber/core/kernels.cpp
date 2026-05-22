#include "kernels.hpp"

#include <cmath>

namespace jabber
{

template<std::size_t TDim, bool TGridInnerLoop>
void ComputeKernel(const std::size_t num_pts, const double rho_bar,
                        const double p_bar, const double *U_bar, 
                        const double gamma, const int num_waves,
                        const double t,
                        const double *__restrict__ rho_coeffs,
                        const double *__restrict__ rhoV_coeffs,
                        const double *__restrict__ rhoE_coeffs, 
                        const double *__restrict__ wave_omegas,
                        const double *__restrict__ k_dot_x_p_phi,
                        double *__restrict__ rho,
                        double *__restrict__ rhoV,
                        double *__restrict__ rhoE)
{
   const double rhoE_init = p_bar/(gamma-1.0);

   if constexpr (TGridInnerLoop)
   {
      for (std::size_t i = 0; i < num_pts; i++)
      {
         rho[i] = rho_bar;

         rhoV[i] = U_bar[0];
         if constexpr(TDim > 1)
         {
            rhoV[num_pts + i] = U_bar[1];
         }
         if constexpr(TDim > 2)
         {
            rhoV[2*num_pts + i] = U_bar[2];
         }
         rhoE[i] = rhoE_init;
      }

      // Add contribution of each wave
#ifdef JABBER_WITH_OPENMP
      #pragma omp parallel for reduction(+:rho[0:num_pts],\
                                          rhoV[0:num_pts*TDim],\
                                          rhoE[0:num_pts])
#endif // JABBER_WITH_OPENMP
      for (int w = 0; w < num_waves; w++)
      {
         const double rho_coeff_w = rho_coeffs[w];
         const double rhoV1_coeff_w = rhoV_coeffs[w];
         const double rhoV2_coeff_w = TDim > 1 ? rhoV_coeffs[num_waves + w] : 0;
         const double rhoV3_coeff_w = TDim > 2 ? rhoV_coeffs[2*num_waves + w] : 0;
         const double rhoE_coeff_w = rhoE_coeffs[w];
         const double omt = wave_omegas[w]*t;

         const std::size_t w_offset = w*num_pts;

         for (std::size_t i = 0; i < num_pts; i++)
         {
            const double cos_w = std::cos(k_dot_x_p_phi[w_offset + i] - omt);

            rho[i] += rho_coeff_w*cos_w;
            rhoV[i] += rhoV1_coeff_w*cos_w;
            if constexpr(TDim > 1)
            {
               rhoV[num_pts + i] += rhoV2_coeff_w*cos_w;
            }
            if constexpr(TDim > 2)
            {
               rhoV[2*num_pts + i] += rhoV3_coeff_w*cos_w;
            }
            rhoE[i] += rhoE_coeff_w*cos_w;

         }
      }
   }
   else
   {
#ifdef JABBER_WITH_OPENMP
      #pragma omp parallel for
#endif // JABBER_WITH_OPENMP
      for (std::size_t i = 0; i < num_pts; i++)
      {
         double rho_i = rho_bar;
         double rhoV1_i = U_bar[0];
         double rhoV2_i = TDim > 1 ? U_bar[1] : 0.0;
         double rhoV3_i = TDim > 2 ? U_bar[2] : 0.0;
         double rhoE_i = rhoE_init;

         const std::size_t i_offset = i*num_waves;
         
         for (int w = 0; w < num_waves; w++)
         {
            const double omt = wave_omegas[w]*t;
            const double cos_w = std::cos(k_dot_x_p_phi[i_offset + w] - omt);

            rho_i += rho_coeffs[w]*cos_w;
            rhoV1_i += rhoV_coeffs[w]*cos_w;
            if constexpr (TDim > 1)
            {
               rhoV2_i += rhoV_coeffs[num_waves + w]*cos_w;
            }
            if constexpr (TDim > 2)
            {
               rhoV3_i += rhoV_coeffs[2*num_waves + w]*cos_w;
            }
            rhoE_i += rhoE_coeffs[w]*cos_w;
         }
         rho[i] = rho_i;
         rhoV[i] = rhoV1_i;
         if constexpr (TDim > 1)
         {
            rhoV[num_pts + i] = rhoV2_i;
         }
         if constexpr (TDim > 2)
         {
            rhoV[2*num_pts + i] = rhoV3_i;
         }
         rhoE[i] = rhoE_i;
      }
   }

   for (std::size_t i = 0; i < num_pts; i++)
   {
      double mag_u = 0.0;
      
      const double val0 = rhoV[i];
      mag_u += val0*val0;
      if constexpr (TDim > 1)
      {
         const double val1 = rhoV[num_pts + i];
         mag_u += val1*val1;
      }
      if constexpr (TDim > 2)
      {
         const double val2 = rhoV[2*num_pts + i];
         mag_u += val2*val2;
      }
      rhoE[i] += 0.5*rho[i]*mag_u;
   }
   
   for (std::size_t i = 0; i < num_pts; i++)
   {
      rhoV[i] *= rho[i];
      if constexpr (TDim > 1)
      {
         rhoV[num_pts + i] *= rho[i];
      }
      if constexpr (TDim > 2)
      {
         rhoV[2*num_pts + i] *= rho[i];
      }
   }
}

// Explicit instantiation for Dims 1-3
template void ComputeKernel<1, true>(const std::size_t, const double,
                                 const double, const double *, 
                                 const double, const int, const double,
                                 const double *__restrict__, 
                                 const double *__restrict__,
                                 const double *__restrict__,
                                 const double *__restrict__,
                                 const double *__restrict__,
                                 double *__restrict__,
                                 double *__restrict__,
                                 double *__restrict__);
                                
template void ComputeKernel<2, true>(const std::size_t, const double,
                                 const double, const double *, 
                                 const double, const int, const double,
                                 const double *__restrict__, 
                                 const double *__restrict__,
                                 const double *__restrict__,
                                 const double *__restrict__,
                                 const double *__restrict__,
                                 double *__restrict__,
                                 double *__restrict__,
                                 double *__restrict__);

template void ComputeKernel<3, true>(const std::size_t, const double,
                                 const double, const double *, 
                                 const double, const int, const double,
                                 const double *__restrict__, 
                                 const double *__restrict__,
                                 const double *__restrict__,
                                 const double *__restrict__,
                                 const double *__restrict__,
                                 double *__restrict__,
                                 double *__restrict__,
                                 double *__restrict__);

template void ComputeKernel<1, false>(const std::size_t, const double,
                                 const double, const double *, 
                                 const double, const int, const double,
                                 const double *__restrict__, 
                                 const double *__restrict__,
                                 const double *__restrict__,
                                 const double *__restrict__,
                                 const double *__restrict__,
                                 double *__restrict__,
                                 double *__restrict__,
                                 double *__restrict__);
                                
template void ComputeKernel<2, false>(const std::size_t, const double,
                                 const double, const double *, 
                                 const double, const int, const double,
                                 const double *__restrict__, 
                                 const double *__restrict__,
                                 const double *__restrict__,
                                 const double *__restrict__,
                                 const double *__restrict__,
                                 double *__restrict__,
                                 double *__restrict__,
                                 double *__restrict__);

template void ComputeKernel<3, false>(const std::size_t, const double,
                                 const double, const double *, 
                                 const double, const int, const double,
                                 const double *__restrict__, 
                                 const double *__restrict__,
                                 const double *__restrict__,
                                 const double *__restrict__,
                                 const double *__restrict__,
                                 double *__restrict__,
                                 double *__restrict__,
                                 double *__restrict__);


} // namespace jabber
