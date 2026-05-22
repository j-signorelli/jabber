#include "config_input.hpp"

#include <toml.hpp>

#include <iomanip>
#include <format>
#include <sstream>
#include <algorithm>

// Helper type for the std::visit
// (https://en.cppreference.com/w/cpp/utility/variant/visit)
template<class... Ts>
struct overloads : Ts... { using Ts::operator()...; };

namespace jabber
{
namespace app
{

/// Tab size used.
static const int kTabSize = 3;

/// Simple type trait for std::vector.
template<typename T>
struct is_std_vector : std::false_type {};

template<typename... VArgs>
struct is_std_vector<std::vector<VArgs...>> : std::true_type {};

/// Get string form of \p val.
template<typename T>
std::string ToString(T val)
{
   std::string out_str;

   // If it's a vector... 
   if constexpr (is_std_vector<T>::value)
   {
      out_str += "[";

      constexpr int kMaxSize = 4;
      static_assert(kMaxSize > 2);

      if (val.size() > kMaxSize)
      {
         // Print first and last elements only if size > kMaxSize
         out_str += ToString(val.front()) + ",...," + ToString(val.back());
      }
      else
      {
         for (int w = 0; w < val.size(); w++)
         {
            out_str += ToString(val[w]);
            out_str += (w+1==val.size() ? "" : ",");
         }
      }
      out_str += "]";
   }
   else
   {
      out_str += std::format("{}", val);
   }

   return out_str;
}

/// Print \p str tabbed at level \p tab_level.
std::string PrintTabbed(const std::string &str, int tab_level)
{
   return std::format("{:<{}}- {}", "", tab_level*kTabSize, str);
}

/// Simple param-value struct for printing.
struct PV
{
   std::string param;
   std::string value;
};

/// Helper for printing parameter information.
std::string PrintParams
   (const std::vector<PV> &params,
      int tab_level=1)
{
   // Get the size of the largest key
   int max_width = 0;
   for (const PV &pv : params)
   {
      const int width = pv.param.size();
      if (width > max_width)
      {
         max_width = width;
      }
   }
   max_width++;

   // Print all params
   std::stringstream out_ss;
   for (const PV &pv : params)
   {  
      out_ss << PrintTabbed(std::format("{:<{}}= {}", pv.param,
                                       max_width, pv.value),  tab_level)
               << std::endl;
   }

   return out_ss.str();
}

void ConfigInput::PrintBaseFlowParams(std::ostream &out) const
{
   out << "Base Flow" << std::endl;
   const std::vector<PV> params
   ({  
      {"rho",     ToString(base_flow_.rho)},
      {"p",       ToString(base_flow_.p)},
      {"U",       ToString(base_flow_.U)},
      {"gamma",   ToString(base_flow_.gamma)}
   });


   out << PrintParams(params) << std::endl;
}

/// Helper name-getter for variant-parameter types.
template<typename T>
std::string GetName(const typename T::Option &option)
{
   return std::string(T::kNames[static_cast<std::size_t>(option)]);
}

/// Print InputXY params visitor.
struct PrintInputXYVisitor
{
   using enum InputXY::Option;

   const int &tab_level;

   std::string operator()(const InputXY::Params<Here> &op)
   {
      const std::vector<PV> params
      ({
         {"Type", GetName<InputXY>(Here)},
         {"X",    ToString(op.x)},
         {"Y",    ToString(op.y)},
      });
      return PrintParams(params, tab_level);
   }

   std::string operator()(const InputXY::Params<FromCSV> &op)
   {
      const std::vector<PV> params
      ({
         {"Type", GetName<InputXY>(FromCSV)},
         {"File", ToString(op.file)},
      });
      return PrintParams(params, tab_level);
   }
};

/// Print FunctionType params visitor.
struct PrintFunctionTypeVisitor
{
   using enum FunctionType::Option;

   const int &tab_level;

   template<FunctionType::Option O>
   std::string operator()(const FunctionType::Params<O> &op)
   {
      std::string out_str;

      const std::vector<PV> params
      ({
         {"Type", GetName<FunctionType>(O)}
      });
      out_str += PrintParams(params, tab_level);
      out_str += PrintTabbed("Input XY\n", tab_level);
      out_str += std::visit(PrintInputXYVisitor{tab_level+1}, 
                                    op.input_xy);
      return out_str;
   }
};

/// Print DiscMethod params visitor.
struct PrintDiscMethodVisitor
{
  using enum DiscMethod::Option;
  
  const int &tab_level;

  template<DiscMethod::Option O>
  std::string operator()(const DiscMethod::Params<O> &op)
  {
      if constexpr (O == Uniform || O == UniformLog)
      {
         const std::vector<PV> params
         ({
            {"Type", GetName<DiscMethod>(O)}
         });
         return PrintParams(params, tab_level);
      }
      else // else Random or RandomLog
      {
         const std::vector<PV> params
         ({
            {"Type", GetName<DiscMethod>(O)},
            {"Seed", ToString(op.seed)}
         });
         return PrintParams(params, tab_level);
      }
  }
};

/// Print Direction params visitor.
struct PrintDirectionVisitor
{
   using enum Direction::Option;

   const int &tab_level;

   std::string operator()(const Direction::Params<Single> &op)
   {
      const std::vector<PV> params
      ({
         {"Type", GetName<Direction>(Single)},
         {"Vector", ToString(op.direction)}
      });
      return PrintParams(params, tab_level);
   }

   std::string operator()(const Direction::Params<RandomXYAngle> &op)
   {
      const std::vector<PV> params
      ({
         {"Type", GetName<Direction>(RandomXYAngle)},
         {"Min Angle", ToString(op.min_angle)},
         {"Max Angle", ToString(op.max_angle)},
         {"Seed", ToString(op.seed)},
      });
      return PrintParams(params, tab_level);
   }
};

/// Print TransferFunction params visitor.
struct PrintTransferFunctionVisitor
{
   using enum TransferFunction::Option;

   const int &tab_level;

   std::string operator()(const TransferFunction::Params<LowFrequencyLimit> &op)
   {
      const std::vector<PV> params
      ({
         {"Type", GetName<TransferFunction>(LowFrequencyLimit)}
      });
      return PrintParams(params, tab_level);
   }

   std::string operator()(const TransferFunction::Params<Input> &op)
   {
      std::string out_str;

      const std::vector<PV> params
      ({
         {"Type", GetName<TransferFunction>(Input)}
      });
      out_str += PrintParams(params, tab_level);
      out_str += PrintTabbed("Input TF\n", tab_level);
      out_str += std::visit(PrintFunctionTypeVisitor{tab_level+1},
                              op.input_tf);
      return out_str;
   }

   std::string operator()(const TransferFunction::Params<FlowNormalFit> &op)
   {
      const std::vector<PV> params
      ({
         {"Type", GetName<TransferFunction>(FlowNormalFit)},
         {"Shock Standoff Distance", ToString(op.shock_standoff_dist)}
      });

      return PrintParams(params, tab_level);
   }
};

/// Print Source params visitor.
struct PrintSourceVisitor
{
   using enum Source::Option;
   
   const int tab_level = 1;

   std::string operator()(const Source::Params<SingleWave> &op)
   {
      const std::vector<PV> params
      ({
         {"Type",       GetName<Source>(SingleWave)},
         {"Amplitude",  ToString(op.amp)},
         {"Frequency",  ToString(op.freq)},
         {"Direction",  ToString(op.direction)},
         {"Phase",      ToString(op.phase)},
         {"Speed",      ToString(op.speed)}
      });
      return PrintParams(params, tab_level);
   }

   std::string operator()(const Source::Params<WaveSpectrum> &op)
   {
      const std::vector<PV> params
      ({
         {"Type",        GetName<Source>(WaveSpectrum)},
         {"Amplitudes",  ToString(op.amps)},
         {"Frequencies", ToString(op.freqs)},
         {"Directions",  ToString(op.directions)},
         {"Phases",      ToString(op.phases)},
         {"Speeds",      ToString(op.speeds)}
      });
      return PrintParams(params, tab_level);
   }

   std::string operator()(const Source::Params<PSD> &op)
   {
      std::string out_str;

      const std::vector<PV> params
      ({
         {"Type",          GetName<Source>(PSD)},
         {"Scale Factor",  ToString(op.scale_fac.value_or(1.0))},
         {"Phase Seed",    ToString(op.phase_seed)},
         {"Speed",         ToString(op.speed)}
      });
      out_str += PrintParams(params, tab_level);
      out_str += PrintTabbed("Discretization\n", tab_level);

      const std::vector<PV> disc_params
      ({
         {"Min",        ToString(op.min_disc_freq)},
         {"Max",        ToString(op.max_disc_freq)},
         {"N",          ToString(op.num_waves)},
         {"Interval",   GetName<IntervalType>(op.int_method)},
      });
      out_str += PrintParams(disc_params, tab_level+1);
      out_str += PrintTabbed("Method\n", tab_level+1);
      out_str += std::visit(PrintDiscMethodVisitor{tab_level+2},
                              op.disc_params);
      out_str += PrintTabbed("Input PSD\n", tab_level);
      out_str += std::visit(PrintFunctionTypeVisitor{tab_level+1},
                              op.input_psd);
      out_str += PrintTabbed("Direction\n", tab_level);
      out_str += std::visit(PrintDirectionVisitor{tab_level+1}, 
                              op.dir_params);
      if (op.tf_params.has_value())
      {
         out_str += PrintTabbed("Transfer Function\n", tab_level);
         out_str += std::visit(PrintTransferFunctionVisitor{tab_level+1},
                                 op.tf_params.value());
      }

      return out_str;
   }

   std::string operator()(const Source::Params<WaveCSV> &op)
   {
      const std::vector<PV> params
      ({
         {"Type",       GetName<Source>(WaveSpectrum)},
         {"File",       op.file}
      });
      
      return PrintParams(params, tab_level);
   }
};

void ConfigInput::PrintSourceParams(std::ostream &out) const
{

   using enum Source::Option;

   out << "Sources" << std::endl;
   for (const Source::ParamsVariant &source : sources_)
   {
      out << std::visit(PrintSourceVisitor{}, source) << std::endl;
   }
}

void ConfigInput::PrintCompParams(std::ostream &out) const
{
   out << "Computation" << std::endl;

   const std::vector<PV> params
   ({  
      {"t0",      ToString(comp_.t0)},
      {"Kernel",  GetName<KernelType>(comp_.kernel)}
   });

   out << PrintParams(params) << std::endl;
}

void ConfigInput::PrintPreciceParams(std::ostream &out) const
{
   if (precice_.has_value())
   {
      out << "preCICE" << std::endl;

      const std::vector<PV> params
      ({  
         {"Participant Name",   precice_->participant_name},
         {"Config File",        precice_->config_file},
         {"Fluid Mesh Name",    precice_->fluid_mesh_name},
         {"Mesh Access Region", ToString(precice_->mesh_access_region)}
      });
      
      out << PrintParams(params) << std::endl;
   }
}

void TOMLConfigInput::ParseBaseFlow
   (std::string toml_string, BaseFlowParams &op)
{
   toml::value in_val = toml::parse_str(toml_string);

   op.rho = in_val.at("rho").as_floating();
   op.p = in_val.at("p").as_floating();
   op.U = toml::get<std::vector<double>>(in_val.at("U"));
   op.gamma = in_val.at("gamma").as_floating();
}

/// Internal helper function here for getting enumerator from string.
template<typename T>
T::Option GetOption(const std::string &name)
{
   using Option = typename T::Option;

   const auto it = std::find(T::kNames.begin(), T::kNames.end(), name);

   Option option;

   if (it == T::kNames.end())
   {
      throw std::invalid_argument(std::format("Invalid input argument: {}",
                                              name));
   }
   else
   {
      option = static_cast<T::Option>(std::distance(T::kNames.begin(), it));
   }
   return option;
}

void TOMLConfigInput::ParseInputXY
   (std::string toml_string, InputXY::ParamsVariant &opv)
{
   using enum InputXY::Option;

   toml::table in_val = toml::parse_str(toml_string).as_table();
   InputXY::Option option = GetOption<InputXY>(in_val.at("Type").as_string());

   if (option == Here)
   {
      InputXY::Params<Here> op;
      op.x = toml::get<std::vector<double>>(in_val.at("X"));
      op.y = toml::get<std::vector<double>>(in_val.at("Y"));
      opv = op;
   }
   else if (option == FromCSV)
   {
      InputXY::Params<FromCSV> op;
      op.file = in_val.at("File").as_string();
      opv = op;
   }
}

void TOMLConfigInput::ParseFunctionType
   (std::string toml_string, FunctionType::ParamsVariant &opv)
{
   using enum FunctionType::Option;

   toml::value in_val = toml::parse_str(toml_string);
   FunctionType::Option option = 
            GetOption<FunctionType>(in_val.at("Type").as_string());

   if (option == PiecewiseLinear)
   {
      FunctionType::Params<PiecewiseLinear> op;
      in_val.at("Data").as_table_fmt().fmt = toml::table_format::multiline;
      ParseInputXY(toml::format(in_val.at("Data")), op.input_xy);
      opv = op;
   }
   else if (option == PiecewiseLogLog)
   {
      FunctionType::Params<PiecewiseLogLog> op;
      in_val.at("Data").as_table_fmt().fmt = toml::table_format::multiline;
      ParseInputXY(toml::format(in_val.at("Data")), op.input_xy);
      opv = op;
   }
}

void TOMLConfigInput::ParseDiscMethod
   (std::string toml_string, DiscMethod::ParamsVariant &opv)
{
   using enum DiscMethod::Option;

   toml::value in_val = toml::parse_str(toml_string);
   DiscMethod::Option option = 
               GetOption<DiscMethod>(in_val.at("Type").as_string());

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
      op.seed = in_val.at("Seed").as_integer();
      opv = op;
   }
   else if (option == RandomLog)
   {
      DiscMethod::Params<RandomLog> op;
      op.seed = in_val.at("Seed").as_integer();
      opv = op;
   }
}

void TOMLConfigInput::ParseDirection
   (std::string toml_string, Direction::ParamsVariant &opv)
{
   using enum Direction::Option;

   toml::value in_val = toml::parse_str(toml_string);
   Direction::Option option = 
                     GetOption<Direction>(in_val.at("Type").as_string());

   if (option == Single)
   {
      Direction::Params<Single> op;
      op.direction = toml::get<std::vector<double>>(in_val.at("Vector"));
      opv = op;
   }
   else if (option == RandomXYAngle)
   {
      Direction::Params<RandomXYAngle> op;
      op.min_angle = in_val.at("MinAngle").as_floating();
      op.max_angle = in_val.at("MaxAngle").as_floating();
      op.seed = in_val.at("Seed").as_integer();
      opv = op;
   }
}

void TOMLConfigInput::ParseTransferFunction
   (std::string toml_string, TransferFunction::ParamsVariant &opv)
{
   using enum TransferFunction::Option;

   toml::value in_val = toml::parse_str(toml_string);

   TransferFunction::Option option = 
            GetOption<TransferFunction>(in_val.at("Type").as_string());

   if (option == LowFrequencyLimit)
   {
      opv = TransferFunction::Params<LowFrequencyLimit>{};
   }
   else if (option == Input)
   {
      TransferFunction::Params<Input> op;
      in_val.at("InputTF").as_table_fmt().fmt = toml::table_format::multiline;
      ParseFunctionType(toml::format(in_val.at("InputTF")), op.input_tf);
      opv = op;
   }
   else if (option == FlowNormalFit)
   {
      TransferFunction::Params<FlowNormalFit> op;
      op.shock_standoff_dist = in_val.at("ShockStandoffDist").as_floating();
      opv = op;
   }
}

void TOMLConfigInput::ParseSource
   (std::string toml_string, Source::ParamsVariant &opv)
{
   using enum Source::Option;

   toml::value in_val = toml::parse_str(toml_string);

   Source::Option option = GetOption<Source>(in_val.at("Type").as_string());

   if (option == SingleWave)
   {
      Source::Params<SingleWave> op;

      op.amp = in_val.at("Amplitude").as_floating();
      op.freq = in_val.at("Frequency").as_floating();
      op.direction = toml::get<std::vector<double>>(
                                       in_val.at("DirVector"));
      op.phase = in_val.at("Phase").as_floating();
      op.speed = *(in_val.at("Speed").as_string().data());

      opv = op;
   }
   else if (option == WaveSpectrum)
   {
      Source::Params<WaveSpectrum> op;
      
      op.amps = toml::get<std::vector<double>>(
                                          in_val.at("Amplitudes"));
      op.freqs = toml::get<std::vector<double>>(
                                          in_val.at("Frequencies"));
      op.directions = toml::get<std::vector<std::vector<double>>>(
                                          in_val.at("DirVectors"));
      op.phases = toml::get<std::vector<double>>(
                                             in_val.at("Phases"));
      std::vector<std::string> speed_strs = 
                  toml::get<std::vector<std::string>>(in_val.at("Speeds"));
      op.speeds.resize(speed_strs.size());
      std::transform(speed_strs.begin(), speed_strs.end(), op.speeds.begin(),
                     [](std::string &s) -> char { return *(s.data()); });
      
      opv = op;
   }
   else if (option == PSD)
   {
      Source::Params<PSD> op;

      if (in_val.contains("ScaleFactor"))
      {
         op.scale_fac = in_val.at("ScaleFactor").as_floating();
      }
      op.phase_seed = in_val.at("PhaseSeed").as_integer();
      op.speed = *(in_val.at("Speed").as_string().data());
      
      in_val.at("InputPSD").as_table_fmt().fmt = toml::table_format::multiline;
      ParseFunctionType(toml::format(in_val.at("InputPSD")), op.input_psd);
      
      toml::value in_disc_val = in_val.at("Discretization");
      op.min_disc_freq = in_disc_val.at("Min").as_floating();
      op.max_disc_freq = in_disc_val.at("Max").as_floating();
      op.num_waves = in_disc_val.at("N").as_integer();
      op.int_method = 
         GetOption<IntervalType>(in_disc_val.at("Interval").as_string());
      
      in_disc_val.at("Method").as_table_fmt().fmt = 
                                                toml::table_format::multiline;
      ParseDiscMethod(toml::format(in_disc_val.at("Method")), op.disc_params);

      in_val.at("Direction").as_table_fmt().fmt = 
                                                toml::table_format::multiline;
      ParseDirection(toml::format(in_val.at("Direction")), op.dir_params);

      if (in_val.contains("TransferFunction"))
      {
         op.tf_params = TransferFunction::ParamsVariant{};
         in_val.at("TransferFunction").as_table_fmt().fmt = 
                                                toml::table_format::multiline;
         ParseTransferFunction(toml::format(in_val.at("TransferFunction")), 
                                 *op.tf_params);
      }
      opv = op;
   }
   else if (option == WaveCSV)
   {
      Source::Params<WaveCSV> op;
      op.file = in_val.at("File").as_string();
      opv = op;
   }
}

void TOMLConfigInput::ParseComputation
   (std::string toml_string, CompParams &op)
{
   toml::value in_val = toml::parse_str(toml_string);
   op.t0 = in_val.at("t0").as_floating();
   op.kernel = GetOption<KernelType>(in_val.at("Kernel").as_string());
}

void TOMLConfigInput::ParsePrecice
   (std::string toml_string, PreciceParams &op)
{
   toml::value in_val = toml::parse_str(toml_string);
   
   op.participant_name = in_val.at("ParticipantName")
                                             .as_string();
   op.config_file = in_val.at("ConfigFile").as_string();
   op.fluid_mesh_name = in_val.at("FluidMeshName").as_string();
   op.mesh_access_region = 
         toml::get<std::vector<double>>(in_val.at("MeshAccessRegion"));
}

TOMLConfigInput::TOMLConfigInput(std::string config_file, std::ostream *out)
{
   toml::value file = toml::parse(config_file);

   // Parse base flow parameters
   toml::value in_base_flow = file.at("BaseFlow");
   ParseBaseFlow(toml::format(in_base_flow), base_flow_);
   if (out)
   {
      PrintBaseFlowParams(*out);
   }

   // Parse the acoustic sources of run
   toml::array in_sources = file.at("Sources").as_array();
   for (const toml::value &in_source : in_sources)
   {
      ParseSource(toml::format(in_source), sources_.emplace_back());
   }
   if (out)
   {
      PrintSourceParams(*out);
   }

   // Parse compute fields of run
   toml::value in_comp = file.at("Computation");
   ParseComputation(toml::format(in_comp), comp_);
   if (out)
   {
      PrintCompParams(*out);
   }

   // Parse preCICE related fields - if exists
   if (file.contains("preCICE"))
   {
      toml::value in_precice = file.at("preCICE");
      precice_ = PreciceParams{};
      ParsePrecice(toml::format(in_precice), *precice_);
      if (out)
      {
         PrintPreciceParams(*out);
      }
   }
}

} // namespace app
} // namespace jabber
