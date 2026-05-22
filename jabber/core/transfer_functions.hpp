#ifndef JABBER_TRANSFER_FUNCTIONS
#define JABBER_TRANSFER_FUNCTIONS

#include <span>

namespace jabber
{
/**
 * @defgroup tf_group Transfer Functions
 * @{
 * 
 * @details @todo \f$\chi=\frac{\left(\frac{p'_{02}}{p_{02}}\right)^2}{\left(\frac{p'}{\bar{p}}\right)^2}\f$
 */

/**
 * @brief Compute the transfer function \f$\chi^*\f$, the analytical 
 * low-frequency limit form given in Equation 15 of Chaudhry and
 * Chandler, 2017.
 * 
 * @param mach_bar      Freestream Mach number.
 * @param gamma         Specific heat ratio.
 * @param speed         Wave speed, 'S' for slow, 'F' for fast.
 * 
 * @return \f$\chi^*\f$
 */
double LowFrequencyLimitTF(double mach_bar, double gamma, char speed);

/**
 * @brief Compute the transfer function \f$\chi\f$ obtained via 
 * re-dimensionalization of a fit of the curve in Figure 14b of 
 * Chaudhry and Chandler, 2017.
 * 
 * @param chi_star        Low frequency limit transfer function, 
 *                        \f$\chi^*\f$. See \ref LowFrequencyLimitTF().
 * @param f_s             Shock stand-off frequency, 
 *                        \f$f_s=\frac{c_0}{2\Delta}\f$.
 * @param freq            Frequency to evaluate transfer function at.
 * 
 * @return \f$\chi(f)\f$
 */
double FlowNormalFitTF(double chi_star, double f_s, double freq);

/// @}
// end of tf_group

} // namespace jabber

#endif // JABBER_TRANSFER_FUNCTIONS
