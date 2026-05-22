#include "psd.hpp"

#include <cmath>
#include <stdexcept>

namespace jabber
{


Interval Interval::ComputeInterval(std::span<const double> freqs,
                                    std::size_t i,
                                    Interval::Method method)
{
   const std::size_t N = freqs.size()-1;

   if (N == 0)
   {
      return Interval{freqs[0], freqs[0]};
   }

   if (method == Interval::Method::Midpoint)
   {
      // Left-boundary
      if (i == 0)
      {
         return Interval{freqs[0],
                         (freqs[0] + freqs[1])/2.0};
      }
      // Right-boundary
      else if (i == N)
      {
         return Interval{(freqs[N-1] + freqs[N])/2.0, 
                         freqs[N]};
      }
      // Interior
      else
      {
         return Interval{(freqs[i] + freqs[i-1])/2.0,
                         (freqs[i] + freqs[i+1])/2.0};
      }
   }
   else if (method == Interval::Method::MidpointLog10)
   {
      // Left-boundary
      if (i == 0)
      {
         return Interval{freqs[0],
                         std::sqrt(freqs[0]*freqs[1])};
      }
      // Right-boundary
      else if (i == N)
      {
         return Interval{std::sqrt(freqs[N]*freqs[N-1]),
                         freqs[N]};
      }
      // Interior
      else
      {
         return Interval{std::sqrt(freqs[i]*freqs[i-1]), 
                         std::sqrt(freqs[i]*freqs[i+1])};
      }
   }
   else
   {
      throw std::logic_error("Invalid/Unimplemented Interval::Method");
   }
}

void BasePSD::Discretize(std::span<const double> freqs,
                           Interval::Method method,
                           std::span<double> powers) const
{
   for (std::size_t i = 0; i < freqs.size(); i++)
   {
      Interval iv = Interval::ComputeInterval(freqs, i, method);
      if (i == 0)
      {
         iv.f_left = Min();
      }
      if (i+1 == freqs.size())
      {
         iv.f_right = Max();
      }
      powers[i] = Integrate(iv.f_left, iv.f_right);
   }
}

double PWLinearPSD::Integrate(double f1, double f2) const
{
   // Determine pieces where f1 and f2 fall within
   auto it = Map().upper_bound(f1);
   if (it == Map().end())
   {
      it = Map().begin();
      it++;
   }

   auto it_last = Map().upper_bound(f2);
   if (it_last == Map().end())
   {
      it_last = std::prev(Map().end());
   }
   auto it_stop = std::next(it_last);

   double integral_val = 0.0;
   double x0 = f1;
   double y0, x1, y1, m;

   // Integrate through
   for(; it != it_stop; it++)
   {
      y0 = operator()(x0);
      m = it->second.m;
      if (it != it_last)
      {
         x1 = it->first;
         y1 = it->second.y;
      }
      else
      {
         x1 = f2;
         y1 = operator()(x1);
      }

      integral_val += 0.5*m*(x1*x1-x0*x0)+(y1-x1*m)*(x1-x0);

      x0 = x1;
   }

   return integral_val;
}

double PWLogLogPSD::Integrate(double f1, double f2) const
{
   // Determine pieces where f1 and f2 fall within
   auto it = Map().upper_bound(f1);
   if (it == Map().end())
   {
      it = Map().begin();
      it++;
   }

   auto it_last = Map().upper_bound(f2);
   if (it_last == Map().end())
   {
      it_last = std::prev(Map().end());
   }
   auto it_stop = std::next(it_last);

   double integral_val = 0.0;
   double x0 = f1;
   double y0, x1, y1, m;

   // Integrate through
   for(; it != it_stop; it++)
   {
      y0 = operator()(x0);
      m = it->second.m;
      if (it != it_last)
      {
         x1 = it->first;
         y1 = it->second.y;
      }
      else
      {
         x1 = f2;
         y1 = operator()(x1);
      }

      // If slope within [-1-1e-8,-1+1e-8], use m=-1 integral
      if (std::abs(m+1) > 1e-8)
      {
         integral_val += y0/(m+1)*(x1*std::pow(x1/x0, m)-x0);
      }
      else
      {
         integral_val += y0*x0*std::log(x1/x0);
      }
      x0 = x1;
   }

   return integral_val;
}

} // namespace jabber
