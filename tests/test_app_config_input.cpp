#ifdef JABBER_WITH_APP

#include "test_app_utils.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/matchers/catch_matchers_vector.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/generators/catch_generators_range.hpp>
#include <jabber/jabber.hpp>

#include <functional>

using namespace jabber;
using namespace jabber::app;
using namespace Catch::Matchers;
using namespace Catch::Generators;
namespace jabber_test
{

TEST_CASE("TOMLConfigInput::ParseBaseFlow", "[App][TOMLConfigInput]")
{
   const double kRho = GENERATE(take(1,random(0.1, 1.0)));
   const double kPBar = GENERATE(take(1,random(0.1, 2000.0)));
   const int kDim = GENERATE(1,2,3);
   const std::vector<double> kUBar = 
                        GENERATE_REF(take(1,chunk(kDim, random(0.0, 600.0))));
   const double kGamma = GENERATE(take(1,random(0.1, 1.5)));

   const std::string base_flow_str = 
   std::format(R"(
               rho={}
               p={}
               U={}
               gamma={}
               )", kRho, kPBar, TOMLWriteValue(kUBar), kGamma);
   
   BaseFlowParams params;
   TOMLConfigInput::ParseBaseFlow(base_flow_str, params);

   CHECK(params.rho == kRho);
   CHECK(params.p == kPBar);
   CHECK_THAT(params.U, Equals(kUBar));
   CHECK(params.gamma == kGamma);
}

TEMPLATE_TEST_CASE_SIG("TOMLConfigInput Parse Options", 
   "[App][TOMLConfigInput]",
   ((typename T,
      void(*Parser)(std::string, typename T::ParamsVariant&)),
      T, Parser),

   (InputXY, TOMLConfigInput::ParseInputXY)
   
   ,(FunctionType, TOMLConfigInput::ParseFunctionType)

   ,(DiscMethod, TOMLConfigInput::ParseDiscMethod)

   ,(Direction, TOMLConfigInput::ParseDirection)

   ,(TransferFunction, TOMLConfigInput::ParseTransferFunction)

   ,(Source, TOMLConfigInput::ParseSource)
   )
{
   using Option = typename T::Option;
   using ParamsVariant = typename T::ParamsVariant;

   const Option kOption = GENERATE(options<Option>());
   const ParamsVariant opv = 
      GENERATE_REF(take(5, random_params<T>(kOption)));

   std::string in_str = TOMLWriteParams<T>(opv);

   ParamsVariant opv_parsed;
   Parser(in_str, opv_parsed);

   TestParamsEqual<T>(opv, opv_parsed);
}

TEST_CASE("TOMLConfigInput::ParseComputation", "[App][TOMLConfigInput]")
{
   const double kT0 = GENERATE(take(1,random(0.0,100.0)));

   const AcousticField::Kernel kKernel = 
                                 GENERATE(options<AcousticField::Kernel>());

   const std::string comp_str = 
      std::format(R"(
                     t0={}
                     Kernel='{}'
                  )", kT0, 
                  KernelType::kNames[static_cast<std::size_t>(kKernel)]);

   CompParams params;
   TOMLConfigInput::ParseComputation(comp_str, params);

   CHECK(params.t0 == kT0);
   CHECK(params.kernel == kKernel);
}

TEST_CASE("TOMLConfigInput::ParsePrecice", "[App][TOMLConfigInput]")
{
   constexpr std::string_view kPartName = "TestParticipant";
   constexpr std::string_view kConfigFile = "TestConfig.xml";
   constexpr std::string_view kFluidMesh = "TestFluidMesh";
   const std::vector<double> kMeshRegion = {-0.01, 1.01, -1.01, 1.01};

   const std::string precice_str = 
      std::format(R"(
                     ParticipantName="{}"
                     ConfigFile="{}"
                     MeshAccessRegion=[{},{},{},{}]
                     FluidMeshName="{}"
                  )", kPartName, kConfigFile, kMeshRegion[0], kMeshRegion[1],
                     kMeshRegion[2], kMeshRegion[3], kFluidMesh);
   
   PreciceParams params;
   TOMLConfigInput::ParsePrecice(precice_str, params);
   CHECK(params.participant_name == kPartName);
   CHECK(params.config_file == kConfigFile);
   CHECK(params.fluid_mesh_name == kFluidMesh);
   CHECK_THAT(params.mesh_access_region, Equals(kMeshRegion));
}

} // jabber_test

#endif // JABBER_WITH_APP
