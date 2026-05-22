#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_vector.hpp>
#include <catch2/generators/catch_generators_all.hpp>

#include <jabber/jabber.hpp>

#include <cmath>

using namespace jabber;
using namespace Catch::Matchers;
using namespace Catch::Generators;

namespace jabber_test
{

TEST_CASE("Read + write Waves", "[Wave]")
{
   constexpr std::size_t kNumWaves = 100;
   const int kDim = GENERATE(1,2,3);

   // Create set of waves
   std::array<Wave, kNumWaves> waves;
   for (std::size_t i = 0; i < kNumWaves; i++)
   {
      waves[i].amplitude = GENERATE(take(1, random(0.0, 10.0)));
      waves[i].frequency = GENERATE(take(1, random(1e3, 500e3)));
      waves[i].phase = GENERATE(take(1, random(0.0, 2*M_PI)));
      waves[i].speed = GENERATE(take(1, random(0, 1))) ? 'S' : 'F';
      waves[i].k_hat = GENERATE_REF(take(1, chunk(kDim, random(0.0,1.0))));
   }

   // Write waves
   std::ostringstream os;
   WriteWaves(waves, os);

   // Read waves
   std::string out_string = os.str();
   std::istringstream is(out_string);
   std::vector<Wave> parsed_waves;
   ReadWaves(is, parsed_waves);

   REQUIRE(parsed_waves.size() == waves.size());
   for (std::size_t i = 0; i < waves.size(); i++)
   {
      CHECK(waves[i].amplitude == parsed_waves[i].amplitude);
      CHECK(waves[i].frequency == parsed_waves[i].frequency);
      CHECK(waves[i].phase == parsed_waves[i].phase);
      CHECK(waves[i].speed == parsed_waves[i].speed);
      CHECK_THAT(waves[i].k_hat, Equals(parsed_waves[i].k_hat));
   }
}

} // jabber_test
