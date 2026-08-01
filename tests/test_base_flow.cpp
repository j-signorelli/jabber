#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <jabber/jabber.hpp>

using namespace jabber;
using namespace Catch::Matchers;

TEST_CASE("BaseFlow initialization", "[BaseFlow]")
{
   BaseFlow flow(1.4, 1.0/1.4, 4.0, 
                  std::array<double,3>({1.0, 2.0, 2.0}));

   CHECK_THAT(flow.c, WithinULP(0.5, 0));
   CHECK_THAT(flow.M, WithinULP(6.0, 0));
}
