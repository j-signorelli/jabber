#include "base_flow.hpp"

#include <cmath>

namespace jabber
{

double ComputeNorm(const std::span<const double> &vec)
{
   double mag_sq = 0.0;
   for (int i = 0; i < vec.size(); i++)
   {
      const double &v = vec[i];
      mag_sq += v*v;
   }
   return std::sqrt(mag_sq);
}

BaseFlow::BaseFlow(const double &gamma_, const double &p_, const double &rho_,
                   const std::span<const double> &U_)
   : gamma(gamma_),
     p(p_),
     rho(rho_),
     U(U_.begin(), U_.end()),
     c(std::sqrt(gamma*p/rho)),
     M(ComputeNorm(U)/c)
{

}

} // namespace jabber