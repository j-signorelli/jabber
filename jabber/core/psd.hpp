#ifndef JABBER_PSD
#define JABBER_PSD

#include "interpolant.hpp"

#include <span>
#include <map>
#include <cstdint>

namespace jabber
{

/**
 * @defgroup psd_group Power Spectral Density (PSD) Discretization
 * @{
 * @brief Tools to compute the energy-conserved powers of a discrete set of
 * frequencies from a PSD.
 *
 * @details
 * ## Theory
 * The power \f$V^2\f$ within generic frequency range \f$[f_1, f_2]\f$ of a
 * continuous, one-sided freestream pressure PSD \f$S_{\infty}(f)\f$ with
 * units \f$(p^\prime/\bar{p})^2/\text{Hz}\f$ is given by
 *
 * \f[
 * V^2_{f_1\text{-}f_2}=\int_{f_1}^{f_2}S(\hat{f})d\hat{f}.
 * \f]
 *
 * As outlined in Appendix B of \cite tam2010, a broadband
 * spectrum of acoustic waves can be formulated by discretizing the PSD
 * into a set of wave frequencies \f$\{f_j\}\f$ as shown below.
 *
 * \image html jabber_psd_disc.png width=400px
 *
 * The cosine wave amplitude is then be computed for each frequency by
 *
 * \f[
 * \frac{1}{2}p^{\prime 2}_j=\int_{f^-_j}^{f^+_j}S_\infty(\hat{f})d\hat{f}.
 * \f]
 *
 * ## Continuous PSD Representation
 * For an arbitrary discretization of frequencies in jabber, the input PSD must
 * be represented in a continuous form, in which case this integral may then be
 * exactly evaluated for each frequency bin. To support this, lightweight
 * classes are provided to formulate a continuous representation of a digitized
 * or discrete PSD in @ref cont_psd_group.
 *
 * ## Interval/Bin Determination
 * After defining a continuous form of a PSD and selecting center frequencies
 * \f$\{f_j\}\f$, the bounds (or interval) for a given center frequency must be
 *  determined. See \ref Interval::Method for the options available.
 */

/**
 * @brief Struct of Δf bin information for PSD discretization.
 */
struct Interval
{
public:
   /**
    * @brief Method for determining interval/bin Δf for a given center
    * frequency in discretization of a PSD.
    */
   enum class Method : std::uint8_t
   {
      /// Compute interval as midpoint between adjacent frequencies.
      Midpoint,

      /**
       * @brief Compute interval as midpoint **on a log10 scale** between
       * adjacent frequencies.
       *
       * @details It is advised to use this if the discretization of
       * frequencies are sampled on a log scale.
       */
      MidpointLog10,

      Size
   };

   /**
    * @brief Compute Interval Δf for given discrete frequency at index \p i
    * in \p freqs using \p method.
    *
    * @details Let \f$N=N_w-1\f$.
    *
    * For \ref Method::Midpoint :
    * \f[
    * \left[f^-_j, f^+_j\right] =
    * \begin{cases}
    * \left[f_{j-1}, f_{j+1}\right] & 0<j<N, \\
    * \left[f_{0}, f_{0} + (f_{1}-f_{0})/2\right]& j=0, \\
    * \left[f_{N-1} + (f_{N}-f_{N-1})/2, f_{N}\right] &j=N.
    * \end{cases}
    * \f]
    *
    * For \ref Method::MidpointLog10 :
    * \f[
    * \left[f^-_j, f^+_j\right] =
    * \begin{cases}
    * \left[\sqrt{f_jf_{j-1}}, \sqrt{f_jf_{j+1}}\right] & 0<j<N, \\
    * \left[f_0, \sqrt{f_0f_1}\right] & j=0, \\
    * \left[\sqrt{f_{N-1}f_N} , f_N \right] & j=N.
    * \end{cases}
    * \f]
    *
    * @param freqs        Span of all discrete frequencies, of size \f$N_w\f$.
    * @param i            Index of frequency to compute interval Δf for.
    * @param method       Method to use for interval computation.
    */
   static Interval ComputeInterval(const std::span<const double> &freqs,
                                   std::size_t i, Interval::Method method);


   /// Left-bound of interval.
   double f_left;

   /// Right-bound of interval.
   double f_right;

   /// Δf
   double DeltaF() const
   {
      return f_right - f_left;
   }
};

/**
 * @defgroup cont_psd_group Continuous PSD Representations
 * @{
 */
/// Base abstract class for a PSD.
class BasePSD
{
public:

   /// Lower frequency bound of PSD.
   virtual double Min() const = 0;

   /// Upper frequency bound of PSD.
   virtual double Max() const = 0;

   /**
    * @brief Compute the power/energy from \p f1 to \p f2
    *
    * @details Note that \p f1 and \p f2 should be within [ \ref Min() ,
    * \ref Max() ].
    */
   virtual double Integrate(double f1, double f2) const = 0;

   /**
    * @brief Compute energy-conserved powers using exact integration.
    *
    * @note Integration is done from \ref Min() to \ref Max().
    *
    * @warning Frequencies \p freqs must be sorted!
    *
    * @param freqs         Input discrete center frequencies in ascending
    *                      order, in range [ \ref Min(), \ref Max() ].
    * @param method        Interval::Method enumerator.
    * @param powers        Output powers.
    */
   void Discretize(const std::span<const double> &freqs, 
                   Interval::Method method,
                   std::span<double> powers) const;

   virtual ~BasePSD() = default;
};

/**
 * @brief Piecewise linear interpolation of discrete PSD data.
 *
 */
class PWLinearPSD : public PWLinear, public BasePSD
{

public:
   /**
    * @brief Construct a new PWLinearPSD object
    *
    * @param freq     Set of discrete frequencies to fit lines in log10 space
    *                 to. PSD bounds \ref Min() and \ref Max() are defined by
    *                 the minimum and maximum discrete frequencies provided.
    * @param psd      PSD associated with each frequency in \p freq.
    */
   PWLinearPSD(const std::span<const double> &freq, 
               std::span<const double> psd)
      : PWLinear(freq, psd) {}

   double Min() const override
   {
      return Map().begin()->first;
   }

   double Max() const override
   {
      return std::prev(Map().end())->first;
   }

   double Integrate(double f1, double f2) const override;

   using BasePSD::Discretize;
};

/**
 * @brief Piecewise log-log interpolation of discrete PSD data.
 *
 * @details If a PSD is digitized from a log-log-scaled plot, this option is
 * advised.
 */
class PWLogLogPSD : public PWLogLog, public BasePSD
{
public:

   /**
    * @brief Construct a new PWLogLogPSD object
    *
    * @param freq     Set of discrete frequencies to fit lines in log space
    *                 to. PSD bounds \ref Min() and \ref Max() are defined by
    *                 the minimum and maximum discrete frequencies provided.
    * @param psd      PSD associated with each frequency in \p freq.
    */
   PWLogLogPSD(const std::span<const double> &freq,
               std::span<const double> psd)
      : PWLogLog(freq, psd) {}

   double Min() const override
   {
      return Map().begin()->first;
   }

   double Max() const override
   {
      return std::prev(Map().end())->first;
   }

   double Integrate(double f1, double f2) const override;

   using BasePSD::Discretize;
};

/// @}
// end of cont_psd_group

/// @}
// end of psd_group

} // namespace jabber

#endif // JABBER_PSD
