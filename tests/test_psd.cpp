#include "test_utils.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_vector.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/generators/catch_generators_range.hpp>

#include <jabber/jabber.hpp>

#include <array>
#include <cmath>
#include <functional>

using namespace jabber;
using namespace Catch::Matchers;
using namespace Catch::Generators;

namespace jabber_test
{


double ExactLinearIntegral(double x0, double y0, double x1, double y1,
                           double a, double b)
{
   const double m = (y1-y0)/(x1-x0);
   // Definite integral of m(x-x0)+y0 computed via Wolfram Alpha
   return 0.5*(b-a)*(m*(a+b-2*x0)+2*y0);
}


double ExactLogLogIntegral(double x0, double y0, double x1, double y1, 
                           double a, double b)
{
   const double m = std::log10(y1/y0)/std::log10(x1/x0);
   if (std::abs(m+1) > std::numeric_limits<double>::epsilon())
   {
      // Definite integral of y0*(x/x0)^m computed via Wolfram Alpha
      return (y0/(m+1))*(b*std::pow(b/x0, m) - a*std::pow(a/x0, m));
   }
   else
   {
      // Definite integral of y0*(x/x0)^(-1) computed via Wolfram Alpha
      return x0*y0*(std::log(b/a));
   }
}

// Verify verification functions above
TEST_CASE("Integration function verification", "[PSD]")
{
   SECTION("Linear integral")
   {
      REQUIRE_THAT(ExactLinearIntegral(1e1,1e-8,1e3,1e-12,1e1,1e3), 
                     WithinULP(4.950495e-6, 0));
   }

   SECTION("Log-log integral")
   {
      SECTION("m != -1")
      {
         REQUIRE_THAT(ExactLogLogIntegral(1e1,1e-8,1e3,1e-12,1e1,1e3),
                      WithinULP(9.9e-8, 0));
      }
      SECTION("m=-1")
      {
         REQUIRE_THAT(ExactLogLogIntegral(1e1, 1e-8, 1e3, 1e-10, 1e1,1e3),
                      WithinULP(1e-7*std::log(100), 0));
      }
   }
}

TEST_CASE("ComputeInterval", "[PSD]")
{

   constexpr std::array<double, 3> kFreqSample={0.1e3, 10e3, 40e3};

   // Use generator to ensure that all interval types are checked!;
   const Interval::Method method = GENERATE(options<Interval::Method>());
   const std::uint8_t m = static_cast<std::uint8_t>(method);
   DYNAMIC_SECTION("Interval Method " << m)
   {
      std::array<Interval, 3> exact;
      if (method == Interval::Method::Midpoint)
      {
         exact[0] = Interval{0.1e3, 5.05e3};
         exact[1] = Interval{5.05e3, 25e3};
         exact[2] = Interval{25e3, 40e3};
      }
      else if (method == Interval::Method::MidpointLog10)
      {
         exact[0] = Interval{0.1e3, 1e3};
         exact[1] = Interval{1e3, 20e3};
         exact[2] = Interval{20e3, 40e3};                   
      }
      else
      {
         const int m = static_cast<int>(method);
         FAIL("No unit test for interval method enumerator " << m);
      }

      for (std::size_t i = 0; i < kFreqSample.size(); i++)
      {
         Interval test_interval = 
                     Interval::ComputeInterval(kFreqSample, i, method);
         CAPTURE(i);
         CHECK(exact[i].f_left == test_interval.f_left);
         CHECK(exact[i].f_right == test_interval.f_right);
      }
   }

}

TEST_CASE("PWLinearPSD integration", "[PSD]")
{
   SECTION("Single segment")
   {
      constexpr std::array<double,2> kX{1.0, 3.0}, 
                                     kY{2.0, 4.0};
      constexpr double kM = (kY[1]-kY[0])/(kX[1]-kX[0]);

      
      const double kA=GENERATE(0.5,1.0,1.5); 
      const double kB=GENERATE(3.5,4.0,4.5);
      if (kA > kB)
      {
         FAIL("Error in test configuration. kA cannot be > kB.");
      }

      // Definite integral of m(x-x0)+y0 computed via Wolfram Alpha
      const double kExact = ExactLinearIntegral(kX[0], kY[0], kX[1], kY[1], 
                                                   kA, kB);

      PWLinearPSD psd(kX,kY);
      CAPTURE(kA, kB);
      REQUIRE(psd.Integrate(kA, kB) == kExact);
   }

   SECTION("Multi-segment")
   {
      constexpr std::array<double,4> kX{1.0, 3.0,  5.0,  7.0}, 
                                     kY{2.0, 4.0, -6.0, -8.0};

      // Generate table
      // First value = a
      // Second value = index of PW interval to compute with at a
      const auto [kA, kAInterval]=GENERATE(table<double, std::size_t>(
                                            {{0.5, 0},
                                             {1.0, 0},
                                             {1.5, 0},
                                             {3.0, 1},
                                             {3.3, 1}}));

      // First value = b
      // Second value = index of PW interval to compute with at b                         
      const auto [kB, kBInterval]=GENERATE(table<double, std::size_t>(
                                            {{3.5, 1},
                                             {5.0, 2},
                                             {5.5, 2},
                                             {7.0, 2},
                                             {8.0, 2}}));
      if (kA > kB)
      {
         FAIL("Error in test configuration. kA cannot be > kB.");
      }

      constexpr std::size_t kNumIntervals = kX.size() -1;

      // Define function pointers of integral for each interval
      const std::array<std::function<double(double,double)>,kNumIntervals>
      kExact =
      [&]()
      {
         std::array<std::function<double(double,double)>,kNumIntervals> 
            exact_funcs;
         for (std::size_t i = 0; i < kNumIntervals; i++)
         {
            exact_funcs[i] = 
            [&,i](double a, double b) -> double
            {
               return ExactLinearIntegral(kX[i], kY[i], kX[i+1], kY[i+1], 
                                             a, b);
            };
         }
         return exact_funcs;
      }();

      double exact = 0.0;
      for (std::size_t i = kAInterval; i <= kBInterval; i++)
      {
         exact += kExact[i](i == kAInterval ? kA : kX[i], 
                            i == kBInterval ? kB : kX[i+1]);
      }
      
      PWLinearPSD psd(kX,kY);
      CAPTURE(kA, kB);
      REQUIRE_THAT(psd.Integrate(kA, kB), WithinULP(exact, 5));
   }
}

TEST_CASE("PWLogLogPSD integration", "[PSD]")
{
   SECTION("Single segment")
   {
      constexpr std::array<double,2> kX{5e3, 20e3}, 
                                     kY{1e-8, 1e-9};

      // Ensure that m != -1
      const double kM = std::log10(kY[1]/kY[0])/std::log10(kX[1]/kX[0]);
      REQUIRE(kM != -1);

      const double kA=GENERATE(2.5e3,5e3,10e3); 
      const double kB=GENERATE(15e3, 20e3, 25e3);
      if (kA > kB)
      {
         FAIL("Error in test configuration. kA cannot be > kB.");
      }

      // Definite integral of y0*(x/x0)^m computed via Wolfram Alpha
      const double kExact = ExactLogLogIntegral(kX[0], kY[0], kX[1], kY[1],
                                                kA, kB);

      PWLogLogPSD psd(kX,kY);
      CAPTURE(kA, kB);
      CAPTURE(psd.Integrate(kA, kB), kExact);
      REQUIRE_THAT(psd.Integrate(kA, kB), WithinULP(kExact,5));
   }

   SECTION("Single segment, m=-1")
   {
      constexpr std::array<double,2> kX{1e3, 10e3}, 
                                     kY{1e-6, 1e-7};
      
      // Ensure that m == -1
      const double kM = std::log10(kY[1]/kY[0])/std::log10(kX[1]/kX[0]);
      REQUIRE(kM == -1);

      const double kA=GENERATE(0.5e3,1e3,2.5e3); 
      const double kB=GENERATE(7.5e3, 10e3, 15e3);

      if (kA > kB)
      {
         FAIL("Error in test configuration. kA cannot be > kB.");
      }

      // Definite integral of y0*(x/x0)^(-1) computed via Wolfram Alpha
      const double kExact = ExactLogLogIntegral(kX[0], kY[0], kX[1], kY[1],
                                                kA, kB);

      PWLogLogPSD psd(kX,kY);
      CAPTURE(kA, kB);
      CAPTURE(psd.Integrate(kA, kB), kExact);
      REQUIRE_THAT(psd.Integrate(kA, kB), WithinULP(kExact,5));
   }

   SECTION("Multi-segment")
   {
      constexpr std::array<double,4> kX{1e3,  10e3, 50e3, 95e3}, 
                                     kY{1e-6, 1e-7, 5e-7, 5e-8};

      // Generate table
      // First value = a
      // Second value = index of PW interval to compute with at a
      const auto [kA, kAInterval]=GENERATE(table<double, std::size_t>(
                                            {{0.5e3, 0},
                                             {1e3, 0},
                                             {5e3, 0},
                                             {10e3, 1},
                                             {15e3, 1}}));

      // First value = b
      // Second value = index of PW interval to compute with at b                         
      const auto [kB, kBInterval]=GENERATE(table<double, std::size_t>(
                                            {{30e3, 1},
                                             {50e3, 2},
                                             {75e3, 2},
                                             {95e3, 2},
                                             {150e3, 2}}));

      if (kA > kB)
      {
         FAIL("Error in test configuration. kA cannot be > kB.");
      }
      
      constexpr std::size_t kNumIntervals = kX.size() -1;

      // Define function pointers of integral for each interval
      const std::array<std::function<double(double,double)>,kNumIntervals>
      kExact =
      [&]()
      {
         std::array<std::function<double(double,double)>,kNumIntervals> 
            exact_funcs;
         for (std::size_t i = 0; i < kNumIntervals; i++)
         {
            exact_funcs[i] = 
            [&,i](double a, double b) -> double
            {
               return ExactLogLogIntegral(kX[i], kY[i], kX[i+1], kY[i+1], 
                                             a, b);
            };
         }
         return exact_funcs;
      }();

      double exact = 0.0;
      for (std::size_t i = kAInterval; i <= kBInterval; i++)
      {
         exact += kExact[i](i == kAInterval ? kA : kX[i], 
                            i == kBInterval ? kB : kX[i+1]);
      }
      
      PWLogLogPSD psd(kX,kY);
      CAPTURE(kA, kB);
      REQUIRE_THAT(psd.Integrate(kA, kB), WithinULP(exact, 5));
   }
}

// BasePSD::Integrate + Interval::ComputeInterval mostly cover PSD disc.
// PWLogLogPSD used for simple test below.
TEST_CASE("PSD Discretization", "[PSD]")
{
   constexpr std::array<double,4> kX{1e3,  10e3, 50e3, 95e3}, 
                                    kY{1e-6, 1e-7, 5e-7, 5e-8};
   
   constexpr std::array<double, 3> kFreqs{7e3, 32e3, 95e3};

   // Manually specify integration bounds for min/max handling.
   // Use Interval::Method::MidpointLog10
   const std::array<Interval, kFreqs.size()> kIntervals
                             {Interval{kX[0], 
                                       std::sqrt(kFreqs[0]*kFreqs[1])},
                              Interval{std::sqrt(kFreqs[0]*kFreqs[1]), 
                                       std::sqrt(kFreqs[1]*kFreqs[2])},
                              Interval{std::sqrt(kFreqs[1]*kFreqs[2]), 
                                       kX[kX.size()-1]}};
   
   PWLogLogPSD psd(kX, kY);

   const std::array<double, kFreqs.size()> kPowersExact =
   [&]()
   {
      std::array<double, kFreqs.size()> powers_exact;
      for (std::size_t i = 0; i < kFreqs.size(); i++)
      {
         // All BasePSD::Integrate types verified in tests above.
         powers_exact[i] = psd.Integrate(kIntervals[i].f_left, 
                                          kIntervals[i].f_right);
      }
      return powers_exact;
   }();

   std::array<double, kFreqs.size()> powers_test;
   psd.Discretize(kFreqs, Interval::Method::MidpointLog10, powers_test);

   for (std::size_t i = 0; i < kFreqs.size(); i++)
   {
      CHECK(kPowersExact[i] == powers_test[i]);
   }

}
} // namespace jabber_test
