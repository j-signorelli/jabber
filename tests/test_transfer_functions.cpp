#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/generators/catch_generators_all.hpp>

#include <jabber/jabber.hpp>

using namespace jabber;
using namespace Catch::Matchers;
namespace jabber_test
{

TEST_CASE("Low-frequency limit", "[TransferFunctions]")
{
    SECTION("Slow")
    {
        const double chi_star = LowFrequencyLimitTF(6.0, 1.4, 'S');
        CHECK_THAT(chi_star, WithinAbs(0.23175337604870483, 1e-14));
    }
    SECTION("Fast")
    {
        const double chi_star = LowFrequencyLimitTF(6.0, 1.4, 'F');
        CHECK_THAT(chi_star, WithinAbs(0.907933119227385, 1e-14));
    }
}

TEST_CASE("Flow-normal fit", "[TransferFunctions]")
{
    const double mach = GENERATE(take(3,random(1.0, 10.0)));
    const double gamma = GENERATE(take(3,random(0.5,1.6)));
    const char speed = GENERATE('S', 'F');

    SECTION("Low-frequency limit")
    {
        const double chi_star = LowFrequencyLimitTF(mach, gamma, speed);
        const double chi = FlowNormalFitTF(chi_star, 1.0, 0.0);

        REQUIRE(chi == chi_star);
    }

    SECTION("Pre-computed from t")
    {
        // Unused, but kept for documentation:
        constexpr std::array<double, 5> t_arr
        {0.46928212, 0.85471849, 0.03781777, 0.09707581, 0.57806404};
        
        constexpr std::array<double, 5> freq_nd_arr
        {0.3824347318693906, 0.6266372340276438, 0.0499083792778991,
         0.1250276007470192, 0.4292429847043615};
        
        constexpr std::array<double, 5> chi_nd_arr
        {4.270193937297611 , 2.2363223163231556, 1.0444846389474434,
         1.28539563880403  , 4.411469866465641};

        const double chi_star = LowFrequencyLimitTF(mach, gamma, speed);

        const double f_s = GENERATE(take(3,random(5e-6,1e-3)));

        for (std::size_t i = 0; i < 5; i++)
        {
            const double chi = FlowNormalFitTF(chi_star, f_s, 
                                                freq_nd_arr[i]*f_s);
            CHECK_THAT(chi/chi_star, WithinAbs(chi_nd_arr[i], 1e-12));
        }
    }

}

} // namespace jabber_test
