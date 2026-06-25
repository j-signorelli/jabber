#ifndef JABBER_KERNELS
#define JABBER_KERNELS

#include <cstddef>

namespace jabber
{

/**
 * @defgroup kernels_group Compute Kernels
 * @{
 *
 * @brief Kernels for high-performance evaluation of an acoustically-disturbed
 * freestream.
 *
 * @details In order to optimize the computation of an acoustically-perturbed
 * freestream with a "large" number of waves \f$N_w\f$ on a grid every
 * timestep, the kernels provided require that the maximal evaluations that can
 * be pre-computed are to reduce repetitive inner-loop FLOPs. Consider the
 * conserved variable equations again derived in \ref env_group to evaluate:
 * \f[
 * \rho=\rho_\infty + \rho'=
 * \rho_\infty+\frac{1}{c_\infty^2}
 *    \sum_{j=1}^{N_w}p'_j\cos\left(\vec{k}_j\cdot\vec{x} +
 *                                     \psi_j-\omega_jt\right),
 * \f]
 * \f[
 * \rho\vec{u}=(\rho_\infty+\rho')(\vec{U}_\infty+\vec{u}')=
 * (\rho_\infty+\rho')\left(\vec{U}_\infty +
 * \frac{(\pm1)}{\rho_\infty c_\infty}
 * \sum_{j=0}^{N_w}\hat{k}_jp'_j\cos\left(\vec{k}_j\cdot\vec{x} +
 *                                           \psi_j-\omega_jt\right)\right),
 * \text{ and}
 * \f]
 * \f[
 * \rho E=\frac{p}{\gamma-1}+\frac{1}{2}\rho||\vec{u}||^2.
 * \f]
 *
 * The following terms can be pre-computed to reduce
 * FLOPs, as they do not depend on time:
 *
 *    - \f$\vec{k}_j\cdot\vec{x} + \psi_j\f$ (**assuming a fixed grid**),
 *    - \f$\omega_j=2\pi f_j\f$,
 *    - in \f$\rho\f$: \f$\frac{1}{c_\infty^2}p'_j\f$,
 *    - in \f$\rho\vec{u}\f$: \f$\frac{\pm1}{\rho_\infty c_\infty}\hat{k}_j\f$,
 *       and
 *    - in \f$\rho E\f$: \f$\frac{1}{\gamma-1}p'_j\f$.
 *
 */


/**
   * @brief Kernel function for evaluating perturbed freestream.
   *
   * @details This function was written carefully for all (inner) loops to
   * *hopefully* be auto-vectorized, but as it goes,
   * <a href="https://pharr.org/matt/blog/2018/04/18/ispc-origins">
   * "auto-vectorization is not a programming model"</a>. If performance is
   * important, please ensure that the appropriate vector compiler flags are
   * set when building and check to see if all loops were successfully
   * vectorized. All inner loops have been verified to be auto-vectorizable by
   * the Intel 2025.3.1 C++ compiler using `-O3 -xhost`.
   *
   * This function is templated with \p TDim for "loop unrolling" on the
   * momentum terms, keeping loops over either the wave axis or the grid
   * point axis (for \p TGridInnerLoop false or true respectively) as inner
   * loops (candidates for auto-vectorization). Note that the storage pattern
   * of \p k_dot_x_p_psi depends on \p TGridInnerLoop for ensuring
   * contiguous memory accesses across hardware threads.
   *
   * When compiled with `JABBER_ENABLE_OPENMP`, for \p TGridInnerLoop true,
   * the outer loop of the series summation (wave axis) is parallelized across
   * OpenMP threads, and for \p TGridInnerLoop false, the outer loop of the
   * series summation (grid point axis) is parallelized across OpenMP threads.
   *
   *
   * @tparam TDim            Physical dimension.
   * @tparam TGridInnerLoop  If true, use grid point axis in series
   *                         summation inner loop. If false, use wave axis.
   *                         This impacts what axis is a candidate for
   *                         auto-vectorized operations.
   *
   * @param num_pts          Number of physical points to evaluate at.
   * @param rho_infty        Base flow density.
   * @param p_infty          Base flow pressure.
   * @param U_infty          Base flow velocity.
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
   * @param k_dot_x_p_psi    \copybrief AcousticField::k_dot_x_p_psi Sized
   *                         \p num_waves x \p num_points with ordering
   *                         [wave][point] for \p TGridInnerLoop true or
   *                         [point][wave] for \p TGridInnerLoop false.
   * @param rho              Output flow density to compute, sized \p num_pts.
   * @param rhoV             Output flow momentum vector to compute, sized
   *                         \p TDim x \p num_pts with ordering [dim][point].
   * @param rhoE             Output flow energy to compute, sized \p num_pts.
*/
template<std::size_t TDim, bool TGridInnerLoop>
void ComputeKernel(const std::size_t num_pts, const double rho_infty,
                   const double p_infty, const double *U_infty,
                   const double gamma, const int num_waves,
                   const double t,
                   const double *__restrict__ rho_coeffs,
                   const double *__restrict__ rhoV_coeffs,
                   const double *__restrict__ rhoE_coeffs,
                   const double *__restrict__ wave_omegas,
                   const double *__restrict__ k_dot_x_p_psi,
                   double *__restrict__ rho,
                   double *__restrict__ rhoV,
                   double *__restrict__ rhoE);

/// @}
// end of kernels_group

} // namespace jabber

#endif // JABBER_KERNELS
