#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <jabber/jabber.hpp>

#include <cmath>

using namespace jabber;
using namespace Catch::Matchers;

namespace jabber_test
{

// Discrete x for interpolation.
static constexpr std::array<double, 4> kX={0.1e3, 1e3, 10e3, 100e3};

// Discrete y for interpolation.
static constexpr std::array<double, 4> kY={1e-8, 5e-8, 7e-9, 2e-14};

// Samples to evaluate at.
static constexpr std::array<double, 6> kXSample={0.05e3, 0.7e3, 5e3, 10e3,
                                                 50e3, 150e3};

// Piecewise interval START to compute in of kX/kY of kXSamples.
static constexpr std::array<std::size_t, 6> kXInterval={0, 0, 1, 2, 2, 2};


double ExactLinear(double x0,double y0,
                   double x1, double y1, double x)
{
   const double m = (y1-y0)/(x1-x0);
   return m*(x-x0) + y0;
}

double ExactLogLog(double x0, double y0, 
                     double x1, double y1, double x)
{
   const double m = log10l(y1/y0)/log10l(x1/x0);
   return y0*powl(x/x0, m);
}

TEST_CASE("Evaluation function verification", "[Interpolant]")
{
   SECTION("Linear")
   {
      REQUIRE_THAT(ExactLinear(1e1, 1e-8, 1e3, 1e-12, 1e2), 
                     WithinULP(9.091e-9, 0));
   }

   SECTION("LogLog")
   {
      REQUIRE_THAT(ExactLogLog(1e1, 1e-8, 1e3, 1e-12, 1e2), 
                     WithinULP(1e-10, 0));
   }
}

TEST_CASE("PWLinear", "[Interpolant]")
{
   const std::array<double, kXSample.size()> kYSample=
   []()
   {
      std::array<double, kXSample.size()> exact_y;
      for (std::size_t i = 0; i < kXSample.size(); i++)
      {
         const std::size_t idx = kXInterval[i];
         exact_y[i] = ExactLinear(kX[idx], kY[idx], 
                                  kX[idx+1], kY[idx+1], kXSample[i]);
      }
      return exact_y;
   }();

   PWLinear interp(kX, kY);
   
   for (std::size_t i = 0; i < kXSample.size(); i++)
   {
      CAPTURE(i);
      CHECK_THAT(interp(kXSample[i]), WithinULP(kYSample[i], 1));
   }
}

TEST_CASE("PWLogLog", "[Interpolant]")
{
   const std::array<double, kXSample.size()> kYSample=
   []()
   {
      std::array<double, kXSample.size()> exact_y;
      for (std::size_t i = 0; i < kXSample.size(); i++)
      {
         const std::size_t idx = kXInterval[i];
         exact_y[i] = ExactLogLog(kX[idx], kY[idx], 
                                  kX[idx+1], kY[idx+1], kXSample[i]);
      }
      return exact_y;
   }();

   PWLogLog interp(kX, kY);

   for (std::size_t i = 0; i < kXSample.size(); i++)
   {
      CAPTURE(i);
      CHECK_THAT(interp(kXSample[i]), WithinULP(kYSample[i], 1));
   }
}

} // namespace jabber_test
