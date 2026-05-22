#include "interpolant.hpp"

#include <cmath>

namespace jabber
{


PWLinear::PWLinear(std::span<const double> x_k, std::span<const double> y_k)
{
   for (std::size_t i = 1; i < x_k.size(); i++)
   {
      const double m = (y_k[i] - y_k[i-1])/(x_k[i] - x_k[i-1]);
      func_map_.insert({x_k[i], Line{ y_k[i], m }});
   }
   func_map_.insert({x_k[0], Line{y_k[0], func_map_[x_k[1]].m}});
}

double PWLinear::operator() (const double &x) const
{
   auto it = func_map_.upper_bound(x);

   // Extrapolate from first/last interior curve
   //  if f falls below/above min/max frequency
   if (it == func_map_.begin())
   {
      it++;
   }
   if (it == func_map_.end())
   {
      it--;
   }

   const double x_2 = it->first;
   const double y_2 = it->second.y;
   const double m = it->second.m;
   it--;
   const double x_1 = it->first;
   const double y_1 = it->second.y;

   return m*(x-x_1) + y_1;
}

PWLogLog::PWLogLog(std::span<const double> x_k, std::span<const double> y_k)
{
   for (std::size_t i = 1; i < x_k.size(); i++)
   {
      const double m = std::log10(y_k[i]/y_k[i-1])/std::log10(x_k[i]/x_k[i-1]);
      func_map_.insert({x_k[i], Line{ y_k[i], m }});
   }
   func_map_.insert({x_k[0], Line{y_k[0], func_map_[x_k[1]].m}});
}

double PWLogLog::operator() (const double &x) const
{
   auto it = func_map_.upper_bound(x);

   // Extrapolate from first/last interior curve
   //  if f falls below/above min/max frequency
   if (it == func_map_.begin())
   {
      it++;
   }
   if (it == func_map_.end())
   {
      it--;
   }

   const double x_2 = it->first;
   const double y_2 = it->second.y;
   const double m = it->second.m;
   it--;
   const double x_1 = it->first;
   const double y_1 = it->second.y;

   return y_1*std::pow(x/x_1, m);
}

} // namespace jabber
