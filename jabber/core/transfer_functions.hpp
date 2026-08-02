#ifndef JABBER_TRANSFER_FUNCTIONS
#define JABBER_TRANSFER_FUNCTIONS

#include <span>

namespace jabber
{
/**
 * @defgroup tf_group Transfer Functions
 * @{
 * @brief Transfer functions \f$\chi=S_0/S_\infty\f$ to recover the freestream
 * PSD from a pitot-measured PSD.
 *
 *
 *
 * \image html jabber_pitot_tf.png width=400px
 *
 * @details Ideally, one has the pressure PSD of the tunnel freestream, such
 * as from DNS of the tunnel walls. However in general, experimentally-measured
 * pressure PSDs from pitot probes placed in the freestream are more readily
 * available at relevant high-frequencies. A number of different
 * phenomena may alter the measured PSD \f$S_0(f)\f$ from the freestream PSD
 * \f$S_\infty(f)\f$.
 *
 * To this end, a transfer function \f$\chi=S_0(f)/S_\infty(f)\f$ can be used
 * to recover the freestream PSD from the experimental one. More details can
 * be found in \cite duan2019, \cite chaudhry2017, and \cite chaudhry2019.
 *
 *
 * @note Transfer functions in literature (see above) have been defined such
 * that \f$S_0(f)\f$ is the PSD taken of the nondimensional pressure
 * \f$p'_{02}/p_{02}\f$ and \f$S_\infty(f)\f$ is the PSD taken of the
 * nondimensional pressure \f$p'/p_\infty\f$.
 */

/**
 * @brief Compute the transfer function \f$\chi^*\f$, the analytical
 * low-frequency limit for flow-parallel disturbances.
 * 
 * 
 * @details Compute \f$\chi^*\f$ given by Equation 15 in \cite chaudhry2017,
 *
 * \f[
 *
 * \chi^*_\pm=\left(\frac{M_\infty^2\pm 2M_\infty \mp1/M_\infty}
 *          {\gamma M_\infty^2-(\gamma-1)/2}\right)^2,
 * \f]
 *
 * where the subscript \f$\pm\f$ corresponds to the slow and fast
 * acoustic wave transfer function respectively.
 *
 * @param mach_bar      Freestream Mach number.
 * @param gamma         Specific heat ratio.
 * @param speed         Wave speed, 'S' for slow, 'F' for fast.
 *
 * @return \f$\chi^*\f$
 */
double LowFrequencyLimitTF(double mach_bar, double gamma, char speed);

/**

/**
 * @brief Compute approximate \f$\chi(f)\f$ for flow-normal
 * disturbances from a re-dimensionalization of the collapsed transfer 
 * function in \cite chaudhry2017.
 *
 * @details Specifically, this function uses a Bezier curve fit of 
 * Figure 14b of \cite chaudhry2017, the collapsed flow-normal disturbance
 * transfer function. This requires Newton's method to determine
 * the curve's parameter t associated with the input frequency
 *  \p freq, which is then used to determine \f$\chi\f$.
 *
 * @param chi_star        Low frequency limit transfer function,
 *                        \f$\chi^*\f$. See \ref LowFrequencyLimitTF().
 * @param f_s             Shock stand-off frequency,
 *                        \f$f_s=\frac{c_0}{2\Delta}\f$, where \f$c_0\f$ is
 *                        the speed of sound at stagnation conditions and
 *                        \f$\Delta\f$ is the shock standoff distance from
 *                        the pitot probe.
 * @param freq            Frequency to evaluate transfer function at.
 *
 * @return \f$\chi(f)\f$
 */
double FlowNormalFitTF(double chi_star, double f_s, double freq);

/// @}
// end of tf_group

} // namespace jabber

#endif // JABBER_TRANSFER_FUNCTIONS
