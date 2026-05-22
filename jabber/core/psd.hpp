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
 * The power within generic frequency range \f$[f_1, f_2]\f$ of a 
 * continuous, one-sided PSD \f$S(f)\f$ is given by
 * 
 * \f[
 * P_{f_1\text{-}f_2}=\int_{f_1}^{f_2}S(f)df.
 * \f]
 * 
 * As outlined in Appendix B of Tam et al., 2010, "Continuation of the Near 
 * Acoustic Field of a Jet to the Far Field. Part I: Theory", a broadband
 * spectrum of acoustic waves can be formulated by discretizing the PSD
 * into a set of wave frequencies \f${f_k}\f$ and ensuring energy conservation
 * by setting their amplitudes according to the power within an interval 
 * \f$\Delta f_k\f$. Note that for a PSD with units \f$V^2/\text{Hz}\f$, after
 * computing a discrete set of powers and applying any transfer function, 
 * the cosine wave amplitude can then be computed by
 * 
 * 
 * \f[
 * V_k=\sqrt{2P_k}.
 * \f]
 * 
 * ## Continuous PSD Representation
 * For an arbitrary discretization of frequencies in Jabber, the input PSD must
 * be represented in a continuous form, in which case this integral may then be
 * exactly evaluated for each frequency bin. To support this, lightweight
 * classes are provided to formulate a continuous representation of a digitized
 * or discrete PSD.
 * 
 * ## Interval/Bin Determination
 * After defining a continuous form of a PSD and selecting center frequencies 
 * \f${f_k}\f$, the bounds (or interval) for a given center frequency must be
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
       * frequencies are taken on a log scale.
       */
      MidpointLog10,

      Size
   };

   /**
    * @brief Compute Interval Δf for given discrete frequency at index \p i
    * in \p freqs using \p method.
    * 
    * @todo Update docs.
    * 
    * @details
    * For \ref Method::Midpoint :
    * \f[
    * \Delta f_k = 
    * \begin{cases}
    * (f_{k+1}-f_{k-1})/2 & 0<k<N, \\
    * (f_1+f_0)/2 & k=0, \\
    * (f_N-f_{N-1})/2 &k=N.
    * \end{cases}
    * \f]
    * 
    * For \ref Method::MidpointLog10 :
    * \f[
    * \Delta f_k =
    * \begin{cases}
    * \sqrt{f_kf_{k+1}}-\sqrt{f_kf_{k-1}} & 0<k<N, \\
    * \sqrt{f_0f_1}-f_0 & k=0, \\
    * f_N-\sqrt{f_Nf_{N-1}} &k=N.
    * \end{cases}
    * \f]
    * 
    * @param freqs        Span of all discrete frequencies, of size \f$N+1\f$.
    * @param i            Index of frequency to compute interval Δf for.
    * @param method       Method to use for interval computation.
    */
   static Interval ComputeInterval(std::span<const double> freqs,
                                    std::size_t i, Interval::Method method);


   /// Left-bound of interval.
   double f_left;

   /// Right-bound of interval.
   double f_right;

   /// Δf
   double DeltaF() const { return f_right - f_left; }
};

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
   void Discretize(std::span<const double> freqs, Interval::Method method,
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
    * @param freq     Set of discrete frequencies to fit lines in log space
    *                 to. PSD bounds \ref Min() and \ref Max() are defined by
    *                 the minimum and maximum discrete frequencies provided.
    * @param psd      PSD associated with each frequency in \p freq.
    */
   PWLinearPSD(std::span<const double> freq, std::span<const double> psd)
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
   PWLogLogPSD(std::span<const double> freq, std::span<const double> psd)
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
// end of psd_group

} // namespace jabber

#endif // JABBER_PSD
