#ifndef JABBER_TEST_APP_UTILS
#define JABBER_TEST_APP_UTILS

#ifdef JABBER_WITH_APP
#include "test_utils.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_vector.hpp>
#include <catch2/generators/catch_generators_all.hpp>

#include <jabber/jabber.hpp>

#include <optional>

namespace jabber_test
{

/**
 * @brief Get random parameters, provided an option.
 * 
 * @details This function relies heavily on Catch2 generator wrapper creators.
 * The Catch2 seed is updated every call to them, so this function can be
 * called more than once on the same option enumerator and will yield a 
 * different output parameter. 
 * 
 */
template<typename T>
T::ParamsVariant GetRandomParams
   (const typename T::Option &option=random_option<typename T::Option>().get())
{
   using namespace jabber::app;
   using namespace Catch::Generators;

   using ParamsVariant = typename T::ParamsVariant;

   ParamsVariant opv;
   if constexpr (std::same_as<T,InputXY>)
   {
      using enum InputXY::Option;
      if (option == Here)
      {
         const int in_size = random(10, 100).get();
         InputXY::Params<Here> op;
         op.x = chunk(in_size, random(0.0, 100e3)).get();
         op.y = chunk(in_size, random(1e-12, 1e-9)).get();
         opv = op;
      }
      else if (option == FromCSV)
      {
         const int file_suffix = random(0,100).get();
         InputXY::Params<FromCSV> op;
         op.file = "test_file." + std::to_string(file_suffix) 
                  + ".csv";
         opv = op;
      }
   }
   else if constexpr (std::same_as<T,FunctionType>)
   {
      using enum FunctionType::Option;

      if (option == PiecewiseLinear)
      {
         FunctionType::Params<PiecewiseLinear> op;
         op.input_xy = GetRandomParams<InputXY>();
         opv = op;
      }
      else if (option == PiecewiseLogLog)
      {
         FunctionType::Params<PiecewiseLogLog> op;
         op.input_xy = GetRandomParams<InputXY>();
         opv = op;
      }
   }
   else if constexpr (std::same_as<T,DiscMethod>)
   {
      using enum DiscMethod::Option;

      if (option == Uniform)
      {
         opv = DiscMethod::Params<Uniform>{};
      }
      else if (option == UniformLog)
      {
         opv = DiscMethod::Params<UniformLog>{};
      }
      else if (option == Random)
      {
         DiscMethod::Params<Random> op;
         op.seed = random(0,100).get();
         opv = op;
      }
      else if (option == RandomLog)
      {
         DiscMethod::Params<RandomLog> op;
         op.seed = random(0,100).get();
         opv = op;
      }
   }
   else if constexpr (std::same_as<T,Direction>)
   {
      using enum Direction::Option;
      if (option == Single)
      {
         Direction::Params<Single> op;
         const int dim = random(1,3).get();
         op.direction = chunk(dim, random(0.0,1.0)).get();
         opv = op;
      }
      else if (option == RandomXYAngle)
      {
         Direction::Params<RandomXYAngle> op;
         op.min_angle = random(-20.0,20.0).get();
         op.max_angle = random(20.0,60.0).get();
         op.seed = random(0,100).get();
         opv = op;
      }
   }
   else if constexpr (std::same_as<T,TransferFunction>)
   {
      using enum TransferFunction::Option;

      if (option == LowFrequencyLimit)
      {
         opv = TransferFunction::Params<LowFrequencyLimit>{};
      }
      else if (option == Input)
      {
         TransferFunction::Params<Input> op;
         op.input_tf = GetRandomParams<FunctionType>();
         opv = op;
      }
      else if (option == FlowNormalFit)
      {
         TransferFunction::Params<FlowNormalFit> op;
         op.shock_standoff_dist = take(1,random(0.5e-3, 2.5e-3)).get();
         opv = op;
      }
   }
   else if constexpr (std::same_as<T,Source>)
   {
      using enum Source::Option;

      if (option == SingleWave)
      {
         Source::Params<SingleWave> op;
         op.amp = random(0.1,10.0).get();
         op.freq = random(500.0,1500.0).get();
         const int dim = random(1,3).get();
         op.direction = chunk(dim, random(0.0, 1.0)).get();
         op.phase = random(10.0, 180.0).get();
         op.speed = random(0,1).get() ? 'S' : 'F';
         opv = op;
      }
      else if (option == WaveSpectrum)
      {
         Source::Params<WaveSpectrum> op;
         const std::size_t num_waves = random(1,20).get();
         for (std::size_t i = 0; i < num_waves; i++)
         {
            Source::Params<SingleWave> wave;
            wave = std::get<Source::Params<SingleWave>>
                        (GetRandomParams<Source>(SingleWave));
            op.amps.push_back(wave.amp);
            op.freqs.push_back(wave.freq);
            op.directions.push_back(wave.direction);
            op.phases.push_back(wave.phase);
            op.speeds.push_back(wave.speed);
         }
         opv = op;
      }
      else if (option == PSD)
      {
         Source::Params<PSD> op;
         op.input_psd = GetRandomParams<FunctionType>();
         op.scale_fac = random(0,1).get() ? 
                           std::optional<double>{random(1.0, 10.0).get()}
                           : std::nullopt;
         op.phase_seed = random(1,100).get();
         op.speed = random(0,1).get() ? 'S' : 'F';
         op.min_disc_freq = random(1e3,10e3).get();
         op.max_disc_freq = random(100e3, 1000e3).get();
         op.num_waves = random(1,100).get();
         op.int_method = random_option<jabber::Interval::Method>().get();
         op.disc_params = GetRandomParams<DiscMethod>();
         op.dir_params = GetRandomParams<Direction>();
         if (take(1,random(0,1)).get())
         {
            op.tf_params = GetRandomParams<TransferFunction>();
         }
         else
         {
            op.tf_params = std::nullopt;
         }

         opv = op;
      }
      else if (option == WaveCSV)
      {
         Source::Params<WaveCSV> op;
         op.file = "test_waves." + std::to_string(random(0,100).get()) 
                     + ".csv";
         opv = op;
      }
   }

   return opv;
}


/// Generalized Catch 2 random parameters generator.
template<typename T>
class RandomParamsGenerator
   : public Catch::Generators::IGenerator<typename T::ParamsVariant>
{
private:
   using Option = typename T::Option;
   using ParamsVariant = typename T::ParamsVariant;

   const std::optional<Option> option_;
   RandomOptionGenerator<Option> ro_gen_;

   ParamsVariant opv_;

public:
   RandomParamsGenerator(const std::optional<Option> option=std::nullopt)
   : option_(option)
   {
      static_cast<void>(next());
   }

   const ParamsVariant& get() const override
   {
      return opv_;
   }

   bool next() override
   {
      Option use_option;
      if (!option_)
      {
         ro_gen_.next();
         use_option = ro_gen_.get();
      }
      else
      {
         use_option = *option_;
      }

      opv_ = GetRandomParams<T>(use_option);

      return true;
   }
};

/**
 * @brief Generator-wrapper creator for random params of 
 * \ref OptionEnum. Use by \c random_params<E>(E::V) or with varying param
 * type by \c random_params<E>().
 */
template<typename T>
Catch::Generators::GeneratorWrapper<typename T::ParamsVariant> random_params
   (std::optional<typename T::Option> V=std::nullopt)
{

   return Catch::Generators::GeneratorWrapper<typename T::ParamsVariant>(
            Catch::Detail::make_unique<RandomParamsGenerator<T>>(V));
}

/// Simple type trait for std::vector.
template<typename T>
struct is_std_vector : std::false_type {};

template<typename... VArgs>
struct is_std_vector<std::vector<VArgs...>> : std::true_type {};

/// Helper for writing to a TOML value.
template<typename T>
std::string TOMLWriteValue(const T &val)
{
   std::string out_str;

   if constexpr (is_std_vector<T>::value)
   {
      out_str += "[";
      for (int w = 0; w < val.size(); w++)
      {
         out_str += TOMLWriteValue(val[w]);
         out_str += (w+1==val.size() ? "]" : ",");
      }
   }
   else if constexpr (std::same_as<T, char> || std::same_as<T, std::string> || 
                        std::same_as<T, std::string_view>)
   {
      out_str += std::format("'{}'", val);
   }
   else
   {
      out_str += std::format("{}", val);
   }
   return out_str;
}

/**
 * @brief Write the given option parameter variant to TOML-format.
 */
template<typename T>
std::string TOMLWriteParams
   (const typename T::ParamsVariant &opv, bool inline_table=false)
{
   using namespace jabber::app;

   std::map<std::string, std::string> out_params;

   using Option = T::Option;

   std::visit(
   [=,&out_params]<Option V>(const typename T::template Params<V> &op)
   {
      out_params.emplace("Type", 
                  TOMLWriteValue(T::kNames[static_cast<std::size_t>(V)]));

      if constexpr (std::same_as<T,InputXY>)
      {
         using enum InputXY::Option;
         if constexpr (V == Here)
         {
            out_params.emplace("X", TOMLWriteValue(op.x));
            out_params.emplace("Y", TOMLWriteValue(op.y));
         }
         else if constexpr (V == FromCSV)
         {
            out_params.emplace("File",  TOMLWriteValue(op.file));
         }
      }
      else if constexpr (std::same_as<T,FunctionType>)
      {
         using enum FunctionType::Option;

         if constexpr (V == PiecewiseLinear || V == PiecewiseLogLog)
         {
            out_params.emplace
               ("Data", TOMLWriteParams<InputXY>(op.input_xy, true));
         }
      }
      else if constexpr (std::same_as<T,DiscMethod>)
      {
         using enum DiscMethod::Option;
         if constexpr (V == Random || V == RandomLog)
         {
            out_params.emplace("Seed",  TOMLWriteValue(op.seed));
         }
      }
      else if constexpr (std::same_as<T,Direction>)
      {
         using enum Direction::Option;
         if constexpr (V == Single)
         {
            out_params.emplace("Vector", TOMLWriteValue(op.direction));
         }
         else if constexpr (V == RandomXYAngle)
         {
            out_params.emplace("MinAngle", TOMLWriteValue(op.min_angle));
            out_params.emplace("MaxAngle", TOMLWriteValue(op.max_angle));
            out_params.emplace("Seed", TOMLWriteValue(op.seed));
         }
      }
      else if constexpr (std::same_as<T,TransferFunction>)
      {
         using enum TransferFunction::Option;

         if constexpr (V == Input)
         {
            out_params.emplace("InputTF", 
                  TOMLWriteParams<FunctionType>(op.input_tf, true));
         }
         else if constexpr (V == FlowNormalFit)
         {
            out_params.emplace("ShockStandoffDist", 
                                 TOMLWriteValue(op.shock_standoff_dist));
         }
      }
      else if constexpr (std::same_as<T,Source>)
      {
         using enum Source::Option;
         if constexpr (V == SingleWave)
         {
            out_params.emplace("Amplitude", TOMLWriteValue(op.amp));
            out_params.emplace("Frequency", TOMLWriteValue(op.freq));
            out_params.emplace("DirVector", TOMLWriteValue(op.direction));
            out_params.emplace("Phase", TOMLWriteValue(op.phase));
            out_params.emplace("Speed", TOMLWriteValue(op.speed));
         }
         else if constexpr (V == WaveSpectrum)
         {
            out_params.emplace("Amplitudes", TOMLWriteValue(op.amps));
            out_params.emplace("Frequencies", TOMLWriteValue(op.freqs));
            out_params.emplace("DirVectors", TOMLWriteValue(op.directions));
            out_params.emplace("Phases", TOMLWriteValue(op.phases));
            out_params.emplace("Speeds", TOMLWriteValue(op.speeds));
         }
         else if constexpr (V == PSD)
         {
            out_params.emplace("InputPSD", 
                        TOMLWriteParams<FunctionType>(op.input_psd, true));
            if (op.scale_fac)
            {
               out_params.emplace("ScaleFactor", 
                                    TOMLWriteValue(op.scale_fac.value()));
            }
            out_params.emplace("PhaseSeed", TOMLWriteValue(op.phase_seed));
            out_params.emplace("Speed", TOMLWriteValue(op.speed));
            out_params.emplace("Discretization.Min", 
                                 TOMLWriteValue(op.min_disc_freq));
            out_params.emplace("Discretization.Max", 
                                 TOMLWriteValue(op.max_disc_freq));
            out_params.emplace("Discretization.N", 
                                 TOMLWriteValue(op.num_waves));
            out_params.emplace("Discretization.Interval", 
                                 TOMLWriteValue(
                                    IntervalType::kNames[
                                       static_cast<std::size_t>(op.int_method)]));
            out_params.emplace("Discretization.Method", 
                        TOMLWriteParams<DiscMethod>(op.disc_params, true));
            out_params.emplace("Direction", 
                        TOMLWriteParams<Direction>(op.dir_params, true));
            
            if (op.tf_params.has_value())
            {
               out_params.emplace("TransferFunction", 
                        TOMLWriteParams<TransferFunction>(*op.tf_params, true));
            }
         }
         else if constexpr (V == WaveCSV)
         {
            out_params.emplace("File", TOMLWriteValue(op.file));
         }
      }

   }, opv);


   std::string out_str;
   if (inline_table)
   {
      out_str += "{";
   }
   for (auto iter = out_params.begin(); iter != out_params.end(); iter++)
   {
      const auto &kv = *iter;

      out_str += kv.first + "=" + kv.second;
      if (inline_table && std::next(iter) != out_params.end())
      {
         out_str += ",";
      }
      else if (!inline_table)
      {
         out_str += "\n";
      }
   }
   if (inline_table)
   {
      out_str += "}";
   }

   return out_str;
}


template<typename T>
void TestParamsEqual
   (const typename T::ParamsVariant &opv1,
      const typename T::ParamsVariant &opv2)
{
   using namespace jabber::app;
   using namespace Catch::Matchers;

   using Option = typename T::Option;

   std::visit(
   []<Option V1, Option V2>(const typename T::template Params<V1> &op1, 
                              const typename T::template Params<V2> &op2)
   {
      // First make sure these are indeed the same types
      if constexpr (V1 != V2)
      {
         CAPTURE(T::kNames[static_cast<std::size_t>(V1)], 
                  T::kNames[static_cast<std::size_t>(V2)]);
         FAIL("Non-matching types!");
      }
      else // else the types match!
      {
         if constexpr (std::same_as<T,InputXY>)
         {
            using enum InputXY::Option;

            if constexpr (V1 == Here)
            {
               CHECK_THAT(op1.x, Equals(op2.x));
               CHECK_THAT(op1.y, Equals(op2.y));
            }
            else if constexpr (V1 == FromCSV)
            {
               CHECK(op1.file == op2.file);
            }
         }
         else if constexpr (std::same_as<T,FunctionType>)
         {
            using enum FunctionType::Option;
            if constexpr (V1 == PiecewiseLinear || V1 == PiecewiseLogLog)
            {
               TestParamsEqual<InputXY>(op1.input_xy, op2.input_xy);
            }
         }
         else if constexpr (std::same_as<T, DiscMethod>)
         {
            using enum DiscMethod::Option;
            if constexpr (V1 == Random || V1 == RandomLog)
            {
               CHECK(op1.seed == op2.seed);
            }
         }
         else if constexpr (std::same_as<T, Direction>)
         {
            using enum Direction::Option;

            if constexpr (V1 == Single)
            {
               CHECK_THAT(op1.direction, Equals(op2.direction));
            }
            else if constexpr (V1 == RandomXYAngle)
            {
               CHECK(op1.min_angle == op2.min_angle);
               CHECK(op1.max_angle == op2.max_angle);
               CHECK(op1.seed == op2.seed);
            }
         }
         else if constexpr (std::same_as<T,TransferFunction>)
         {
            using enum TransferFunction::Option;
            if constexpr (V1 == Input)
            {
               TestParamsEqual<FunctionType>(op1.input_tf, op2.input_tf);
            }
            else if constexpr (V1 == FlowNormalFit)
            {
               CHECK(op1.shock_standoff_dist == op2.shock_standoff_dist);
            }
         }
         else if constexpr (std::same_as<T,Source>)
         {
            using enum Source::Option;
            if constexpr (V1 == SingleWave)
            {
               CHECK(op1.amp == op2.amp);
               CHECK(op1.freq == op2.freq);
               CHECK_THAT(op1.direction, Equals(op2.direction));
               CHECK(op1.phase == op2.phase);
               CHECK(op1.speed == op2.speed);
            }
            else if constexpr (V1 == WaveSpectrum)
            {
               CHECK_THAT(op1.amps, Equals(op2.amps));
               CHECK_THAT(op1.freqs, Equals(op2.freqs));
               CHECK_THAT(op1.phases, Equals(op2.phases));
               CHECK_THAT(op2.speeds, Equals(op2.speeds));
               REQUIRE(op1.directions.size() == op2.directions.size());
               for (std::size_t i = 0; i < op1.directions.size(); i++)
               {
                  CHECK_THAT(op1.directions[i], Equals(op2.directions[i]));
               }
            }
            else if constexpr (V1 == PSD)
            {
               CHECK(op1.scale_fac == op2.scale_fac);
               CHECK(op1.min_disc_freq == op2.min_disc_freq);
               CHECK(op1.max_disc_freq == op2.max_disc_freq);
               CHECK(op1.num_waves == op2.num_waves);
               CHECK(op1.int_method == op2.int_method);
               CHECK(op1.phase_seed == op2.phase_seed);
               CHECK(op1.speed == op2.speed);
               TestParamsEqual<FunctionType>(op1.input_psd, op2.input_psd);
               TestParamsEqual<DiscMethod>(op1.disc_params, op2.disc_params);
               TestParamsEqual<Direction>(op1.dir_params, op2.dir_params);
                                                
               REQUIRE(op1.tf_params.has_value() == op2.tf_params.has_value());
               if (op1.tf_params.has_value())
               {
                  TestParamsEqual<TransferFunction>
                     (*op1.tf_params, *op2.tf_params);
               }

            }
            else if constexpr (V1 == WaveSpectrum)
            {
               CHECK(op1.file == op2.file);
            }
         }
      }
   }, opv1, opv2);
}

} // namespace jabber_test

#endif // JABBER_WITH_APP

#endif // JABBER_TEST_APP_UTILS