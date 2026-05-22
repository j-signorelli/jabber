#ifdef JABBER_WITH_APP

#include "test_utils.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_vector.hpp>
#include <catch2/generators/catch_generators_all.hpp>

#include <jabber/jabber.hpp>

#include <filesystem>
#include <fstream>
#include <numeric>
#include <vector>

using namespace jabber::app;
using namespace Catch::Matchers;
using namespace Catch::Generators;

namespace jabber_test
{

TEST_CASE("Normalize", "[App][Common]")
{   
   constexpr std::size_t kDim = 3;
   std::vector<double> vec = GENERATE(take(1,chunk(kDim,random(0.0,1.0))));
   std::vector<double> norm_vec(kDim);

   Normalize(vec, norm_vec);

   // Check vec dot norm_vec == ||vec||
   const double vec_mag = std::sqrt(std::inner_product(vec.begin(), 
                                                      vec.end(), 
                                                      vec.begin(), 
                                                      0.0));
   const double dot_product = std::inner_product(vec.begin(), 
                                                   vec.end(), 
                                                   norm_vec.begin(),
                                                   0.0);
   CAPTURE(vec, norm_vec);
   REQUIRE_THAT(dot_product, WithinRel(vec_mag, 1e-12));
}

TEST_CASE("GetRankPartition", "[App][Common]")
{
   const int kVDim = GENERATE(1,2,3);
   const int kN = GENERATE(take(2,random(10,100)));

   const int kSize = GENERATE_REF(range(1, 2*kN));

   CAPTURE(kVDim, kN, kSize);
   
   const std::vector<double> values = take(1,
                                          chunk(kN*kVDim,random(-10.0,10.0)))
                                          .get();

   // Test by appending all subspans together + comparing against original
   std::span<const double> r_sub_values;
   std::vector<double> values_test;
   for (std::size_t r = 0; r < kSize; r++)
   {
      GetRankPartition<double>(values, kVDim, r, kSize, r_sub_values);

      values_test.insert(values_test.end(), 
                           r_sub_values.begin(),
                           r_sub_values.end());
   }

   REQUIRE_THAT(values, Equals(values_test));
}

// TEST_CASE("InputXYVisitor", "[App]")
// {
//    constexpr int kSeed = 0;
//    const InputXY option = GENERATE(options<InputXYOption>());
//    const std::uint8_t idx = static_cast<std::uint8_t>(option);

//    constexpr std::size_t kN = 20;

//    DYNAMIC_SECTION(InputXYNames[idx])
//    {
//       InputXYParamsVariant ipv = GenerateRandomInputXY(option, kSeed);
   
//       std::vector<double> x = GenerateRandomVec<kN>(kSeed, 10, 100e3);
//       std::vector<double> y = GenerateRandomVec<kN>(kSeed, 1e-9, 1e-12);

//       // Pre-process
//       std::visit(
//       overloads
//       {
//       [&](InputXYParams<InputXY::Here> &ip)
//       {
//          ip.x = x;
//          ip.y = y;
//       },
//       [&](InputXYParams<InputXY::FromCSV> &ip)
//       {
//          // Write sample CSV file.
//          std::string file_name = std::filesystem::temp_directory_path() 
//                                  / "jabber_unit_test_XXXXXX";
//          mkstemp(file_name.data());

//          std::ofstream file(file_name);
//          for (std::size_t i = 0; i < kN; i++)
//          {
//             file << std::format("{},{}", x[i], y[i]) << std::endl;
//          }

//          ip.file = file_name;
//       }
//       }, ipv);

//       // Apply the visit
//       std::vector<double> x_test, y_test;
//       std::visit(InputXYVisitor{x_test, y_test}, ipv);

//       CHECK_THAT(x, Equals(x_test));
//       CHECK_THAT(y, Equals(y_test));

//       // Delete any written files...
//       std::visit(
//       [&]<InputXY I>(const InputXYParams<I> &ip)
//       {
//          if constexpr (I == InputXY::FromCSV)
//          {
//             std::filesystem::remove(ip.file);
//          }
//       }, ipv);
//    }
// }

} // jabber_test

#endif // JABBER_WITH_APP
