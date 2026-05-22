#ifndef JABBER_KERNELS
#define JABBER_KERNELS

#include <cstddef>

namespace jabber
{

/**
 * @defgroup kernels_group Compute Kernels
 * @{
 * 
 */
   

/**
   * @brief Kernel function for evaluating perturbed base flow, with series
   * summation inner loop/vectorization over each gridpoint.
   * 
   * @details This function was designed following Intel guidelines for 
   * auto-vectorizable code.  Because only inner-most loops are candidates for
   * vectorization, this function is templated with \p TDim for "loop
   * unrolling" on the momentum terms. Vectorization is across \p num_pts as
   * that is expected to be the largest value, so data dimensioned by it must
   * be stored in an SoA-format for contiguous memory accesses across hardware
   * threads.
   * 
   * All inner loops have been verified to be vectorized by Intel `icpx` 
   * 2025.3.1 using the flags `-O3 -xhost`. Proper vectorization by Intel
   * compilers can be checked via:
   * 
   * ```
   * icpx -O3 -qopt-report=3 \
   *          -qopt-report-file="report.yaml" \
   *          -xhost \
   *          -c kernels.cpp
   * ```
   * 
   * 
   * @tparam TDim            Physical dimension.
   * @tparam TGridInnerLoop  If true, use grid point axis in series
   *                         summation inner loop. If false, use wave axis.
   *                         This impacts vectorization.
   * 
   * @param num_pts          Number of physical points to evaluate at.
   * @param rho_bar          Base flow density.
   * @param p_bar            Base flow pressure.
   * @param U_bar            Base flow velocity.
   * @param gamma            Specific heat ratio.
   * @param num_waves        Number of acoustic waves to compute.
   * @param t                Time.
   * @param rho_coeffs       \copybrief AcousticField::rho_coeffs Sized
   *                         \p num_waves.
   * @param rhoV_coeffs      \copybrief AcousticField::rhoV_coeffs Sized
   *                         \p TDim x \p num_waves.
   * @param rhoE_coeffs      \copybrief AcousticField::rhoE_coeffs Sized
   *                         \p num_waves.
   * @param wave_omegas      \copybrief AcousticField::wave_omegas Sized 
   *                         \p num_waves.
   * @param k_dot_x_p_phi    \copybrief AcousticField::k_dot_x_p_phi Sized
   *                         \p num_waves x \p num_points with ordering 
   *                         [wave][point] for \p TGridInnerLoop true or
   *                         [point][wave] for \p TGridInnerLoop false.
   * @param rho              Output flow density to compute, sized \p num_pts.
   * @param rhoV             Output flow momentum vector to compute, sized
   *                         \p TDim x \p num_pts with ordering [dim][point].
   * @param rhoE             Output flow energy to compute, sized \p num_pts.
*/
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
                        double *__restrict__ rhoE);

/// @}
// end of kernels_group

} // namespace jabber

#endif // JABBER_KERNELS
