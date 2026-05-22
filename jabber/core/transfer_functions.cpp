#include "transfer_functions.hpp"

#include <random>
#include <format>

namespace jabber
{

double LowFrequencyLimitTF(double mach_bar, double gamma, char speed)
{
   const double sign = speed == 'S' ? -1 : 1;
   const double num = (mach_bar*mach_bar + sign*2.0*mach_bar - sign/mach_bar);
   const double denom = (gamma*mach_bar*mach_bar - (gamma - 1.0)/2.0);
   const double frac = num/denom;

   return frac*frac;
}

struct FlowNormalBezierFit
{
   static constexpr int kN = 9;

   /// Control points for frequencies for the flow-normal TF Bezier curve fit.
   static constexpr std::array<double, kN> kControlFreqs
   {{0.0, 0.16204, 0.36886, 0.27954, 0.47948, 0.36266, 0.50724, 0.63124, 
      0.80111}};

   /**
    * @brief Control points for transfer functions for the flow-normal TF
    * Bezier curve fit.
    */
   static constexpr std::array<double, kN> kControlChis
   {{1.0, 1.0, 2.1021, 4.6026, 4.6111, 6.9099, 2.5099, 1.6516, 1.3674}};

   static constexpr double kfMin = kControlFreqs.front();
   static constexpr double kfMax = kControlFreqs.back();

   /**
    * @brief Compute the Bezier curve defined by the points within
    * [ \p IMin , \p IMax ], using de Casteljau's algorithm.
    */
   template<std::size_t IMin, std::size_t IMax>
   static double B(const std::array<double, kN> controlPts, double t)
   {
      if constexpr (IMin == IMax)
      {
         return controlPts[IMin];
      }
      else
      {
         return (1-t)*B<IMin,IMax-1>(controlPts, t)
                  + t*B<IMin+1,IMax>(controlPts, t);
      }
   }
   
   /// Evaluate the flow normal TF curve-fit frequency at a given t.
   static double EvalFreq(double t)
   {
      return B<0,kN-1>(kControlFreqs, t);
   }

   /// Evaluate the flow normal TF curve-fit \f$\chi\f$ at a given t.
   static double EvalChi(double t)
   {
      return B<0,kN-1>(kControlChis, t);
   }

   /// Evaluate the derivative of frequency of the curve fit at a given t.
   static double DEvalFreq(double t)
   {
      return kN*(B<1,kN-1>(kControlFreqs, t) - B<0,kN-2>(kControlFreqs, t));
   }
};

double FlowNormalFitTF(double chi_star, double f_s, double freq)
{
   const double &f_norm = freq/f_s;
   if (f_norm < FlowNormalBezierFit::kfMin || 
         f_norm > FlowNormalBezierFit::kfMax)
   {
      std::string err_string = 
         std::format("Normalized frequency {} (non-normalized frequency {})"
                     " must be within range [{},{}].",
                     f_norm, freq, FlowNormalBezierFit::kfMin, 
                     FlowNormalBezierFit::kfMax);

      throw std::invalid_argument(err_string);
   }

   // Compute Bezier curve parameter t for freq using Newton's method
   constexpr std::size_t kNumNewtonIt = 100;
   constexpr double kNewtonTol = 1e-12;
   constexpr std::size_t kMinNewtonIt = 10;

   double tn = 0.5;
   double tnp1;
   for (std::size_t n = 0; n < kNumNewtonIt; n++)
   {
      tnp1 = tn - (FlowNormalBezierFit::EvalFreq(tn) - f_norm)
                     / FlowNormalBezierFit::DEvalFreq(tn);
      if (tnp1 < 0)
      {
         tnp1 = 0.0;
      }
      else if (tnp1 > 1)
      {
         tnp1 = 1.0;
      }

      // If converged, end loop
      if (n >= kMinNewtonIt && std::abs(tnp1-tn) < kNewtonTol)
      {
         n = kNumNewtonIt;
      }
      else if (n + 1 == kNumNewtonIt)
      {
         std::string err_string = 
            std::format("Unable to converge to Bezier parameter t for "
                        "frequency {} (normalized frequency {})", 
                        freq, f_norm);
         throw std::logic_error(err_string);
      }
      else
      {
         tn = tnp1;
      }
   }

   return FlowNormalBezierFit::EvalChi(tnp1)*chi_star;
}

} // namespace jabber
