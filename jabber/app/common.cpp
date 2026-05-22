#include "common.hpp"

#include <cmath>
#include <memory>
#include <random>
#include <fstream>
#include <ranges>
#include <algorithm>
#include <type_traits>

// Helper type for the std::visit
// (https://en.cppreference.com/w/cpp/utility/variant/visit)
template<class... Ts>
struct overloads : Ts... { using Ts::operator()...; };


namespace jabber
{
namespace app
{

void PrintBanner(std::ostream &out)
{
   constexpr std::string_view banner = R"(
      /-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\
      |   ________  ____    ______  ______    _____ ______      |
      \  (___  ___)(    )  (_   _ \(_   _ \  / ___/(   __ \     /
      -      ) )   / /\ \    ) (_) ) ) (_) )( (__   ) (__) )    -
      /     ( (   ( (__) )   \   _/  \   _/  ) __) (    __/     \
      |  __  ) )   )    (    /  _ \  /  _ \ ( (     ) \ \  _    |
      \ ( (_/ /   /  /\  \  _) (_) )_) (_) ) \ \___( ( \ \_))   /
      -  \___/   /__(  )__\(______/(______/   \____\)_) \__/    -
      /                                                         \
      |                                                         |
      \-/|\-/|\-/|\-/|\-/|\-/|\-/|\-/|\-/|\-/|\-/|\-/|\-/|\-/|\-/)";

   out << banner << std::endl << std::endl;
}

void Normalize(std::span<const double> vec, std::span<double> norm_vec)
{
   double sum_sq = 0.0;
   for (std::size_t i = 0; i < vec.size(); i++)
   {
      const double val = vec[i];
      sum_sq += val*val;
   }

   const double mag = std::sqrt(sum_sq);
   for (std::size_t i = 0; i < vec.size(); i++)
   {
      norm_vec[i] = vec[i]/mag;
   }
}


void InputXYVisitor::operator()
   (const InputXY::Params<Here> &op)
{
   x = op.x;
   y = op.y;
}

void InputXYVisitor::operator()
   (const InputXY::Params<FromCSV> &op)
{
   std::ifstream is(op.file);
   if (!is.is_open())
   {
      throw std::invalid_argument(std::format("Cannot find XY-data CSV file "
                                              "'{}'.", op.file));
   }
   for (std::string line; std::getline(is, line);)
   {
      auto row_view = std::ranges::views::split(line,std::string_view(","));
      auto row_it = row_view.begin();
      auto val_it = *row_it;
      x.push_back(std::stod(std::string(val_it.begin(), 
                                                val_it.end())));
      val_it = *(++row_it);
      y.push_back(std::stod(std::string(val_it.begin(), 
                                                val_it.end())));
   }
}

void FunctionTypeVisitor::operator()
   (const FunctionType::Params<PiecewiseLinear> &op)
{
   std::vector<double> x, y;
   std::visit(InputXYVisitor{x,y}, op.input_xy);

   std::visit(
   overloads
   {
   [&](std::unique_ptr<Function>* T_ptr_ptr)
   {
      *T_ptr_ptr = std::make_unique<PWLinear>(x,y);
   },
   [&](std::unique_ptr<BasePSD>* T_ptr_ptr)
   {
      *T_ptr_ptr = std::make_unique<PWLinearPSD>(x,y);
   }
   }, T_ptr_ptr_var);
   
}

void FunctionTypeVisitor::operator()
   (const FunctionType::Params<PiecewiseLogLog> &op)
{
   std::vector<double> x, y;
   std::visit(InputXYVisitor{x,y}, op.input_xy);

   std::visit(
   overloads
   {
   [&](std::unique_ptr<Function>* T_ptr_ptr)
   {
      *T_ptr_ptr = std::make_unique<PWLogLog>(x,y);
   },
   [&](std::unique_ptr<BasePSD>* T_ptr_ptr)
   {
      *T_ptr_ptr = std::make_unique<PWLogLogPSD>(x,y);
   }
   }, T_ptr_ptr_var);

}

void DiscMethodVisitor::operator()
   (const DiscMethod::Params<Uniform> &op)
{
   const double df = (max_freq - min_freq)/(freqs.size()+1);
   for (std::size_t i = 0; i < freqs.size(); i++)
   {
      freqs[i] = min_freq + df*i;
   }
}

void DiscMethodVisitor::operator()
   (const DiscMethod::Params<UniformLog> &op)
{
   const double log_df = std::log10(max_freq/min_freq)/(freqs.size()+1);
   const double log_min_freq = std::log10(min_freq);
   for (std::size_t i = 0; i < freqs.size(); i++)
   {
      freqs[i] = std::pow(10.0, std::log10(min_freq) + log_df*i);
   }
}

void DiscMethodVisitor::operator()
   (const DiscMethod::Params<Random> &op)
{
   std::mt19937 gen(op.seed);
   std::uniform_real_distribution<double> real_dist(min_freq, max_freq);
   for (std::size_t i = 0; i < freqs.size(); i++)
   {
      freqs[i] = real_dist(gen);
   }
}

void DiscMethodVisitor::operator()
   (const DiscMethod::Params<RandomLog> &op)
{
   std::mt19937 gen(op.seed);
   std::uniform_real_distribution<double> real_dist(std::log10(min_freq),
                                                      std::log10(max_freq));
   for (std::size_t i = 0; i < freqs.size(); i++)
   {
      freqs[i] = std::pow(10.0, real_dist(gen));
   }
}

void DirectionVisitor::operator()
   (const Direction::Params<Single> &op)
{
   for (std::size_t i = 0; i < k_hats.size(); i++)
   {
      k_hats[i] = op.direction;
   }
}

void DirectionVisitor::operator()
   (const Direction::Params<RandomXYAngle> &op)
{
   std::mt19937 dir_gen(op.seed);
   std::uniform_real_distribution<double> dir_dist(
                                    op.min_angle*M_PI/180.0,
                                    op.max_angle*M_PI/180.0);
   for (std::size_t i = 0; i < k_hats.size(); i++)
   {
      k_hats[i].resize(3, 0.0);
      const double angle = dir_dist(dir_gen);
      k_hats[i][0] = std::cos(angle);
      k_hats[i][1] = std::sin(angle);
   }
}

void TransferFunctionVisitor::operator()
   (const TransferFunction::Params<LowFrequencyLimit> &op)
{
   const double &p_bar = base_flow_params.p;
   const double &rho_bar = base_flow_params.rho;
   const double &gamma = base_flow_params.gamma;
   const double c_bar = std::sqrt(gamma*p_bar/rho_bar);

   const double vel_bar = 
   [&]()
   {
      double mag_sq = 0.0;
      for (int d = 0; d < base_flow_params.U.size(); d++)
      {
         const double &v = base_flow_params.U[d];
         mag_sq += v*v;
      }
      return std::sqrt(mag_sq);
   }();

   const double mach_bar = vel_bar/c_bar;

   const double chi_star = LowFrequencyLimitTF(mach_bar, gamma, speed);

   for (std::size_t i = 0; i < freqs.size(); i++)
   {
      powers[i] /= chi_star;
   }
}

void TransferFunctionVisitor::operator()
   (const TransferFunction::Params<Input> &op)
{
   std::unique_ptr<Function> tf;
   std::visit(FunctionTypeVisitor{&tf}, op.input_tf);

   for (std::size_t i = 0; i < freqs.size(); i++)
   {
      powers[i] /= (*tf)(freqs[i]);
   }
}

void TransferFunctionVisitor::operator()
   (const TransferFunction::Params<FlowNormalFit> &op)
{
   const double &p_bar = base_flow_params.p;
   const double &rho_bar = base_flow_params.rho;
   const double &gamma = base_flow_params.gamma;
   const double c_bar = std::sqrt(gamma*p_bar/rho_bar);

   const double vel_bar = 
   [&]()
   {
      double mag_sq = 0.0;
      for (int d = 0; d < base_flow_params.U.size(); d++)
      {
         const double &v = base_flow_params.U[d];
         mag_sq += v*v;
      }
      return std::sqrt(mag_sq);
   }();

   const double mach_bar = vel_bar/c_bar;

   const double chi_star = LowFrequencyLimitTF(mach_bar, gamma, speed);

   const double f_s = 
   [&]()
   {
      const double RT = p_bar/rho_bar;
      const double RT_0 = RT*(1 + ((gamma-1.0)/2)*mach_bar*mach_bar);

      const double c_0 = std::sqrt(gamma*RT_0);

      return c_0/(2*op.shock_standoff_dist);
   }();

   for (std::size_t i = 0; i < freqs.size(); i++)
   {
      powers[i] /= FlowNormalFitTF(chi_star, f_s, freqs[i]);
   }
}

void SourceVisitor::operator()
   (const Source::Params<SingleWave> &op)
{
   std::vector<double> k_hat(op.direction.size(), 0.0);
   Normalize(op.direction, k_hat);
   waves.emplace_back(op.amp, op.freq, op.phase*M_PI/180.0, op.speed, k_hat);
}

void SourceVisitor::operator()
   (const Source::Params<WaveSpectrum> &op)
{
   for (int i = 0; i < op.amps.size(); i++)
   {
      std::vector<double> k_hat(op.directions[i].size(), 0.0);
      Normalize(op.directions[i], k_hat);
      waves.emplace_back(op.amps[i], op.freqs[i], op.phases[i]*M_PI/180.0, 
                           op.speeds[i], k_hat);
   }
}

void SourceVisitor::operator()
   (const Source::Params<PSD> &op)
{
   // Process input PSD
   std::unique_ptr<BasePSD> psd;
   std::visit(FunctionTypeVisitor{&psd}, op.input_psd);

   // Apply the discretization method for selecting center frequencies
   std::vector<double> freqs(op.num_waves);
   const double min_freq = op.min_disc_freq;
   const double max_freq = op.max_disc_freq;
   std::visit(DiscMethodVisitor{min_freq,max_freq,freqs}, op.disc_params);

   // Sort the frequencies
   std::sort(freqs.begin(), freqs.end());

   // Compute the powers of each wave
   std::vector<double> powers(freqs.size());
   psd->Discretize(freqs, op.int_method, powers);

   // Apply transfer function to the powers
   if (op.tf_params.has_value())
   {
      std::visit(TransferFunctionVisitor{base_flow_params, freqs, op.speed,
                                          powers},
                  *op.tf_params);
   }
   // Compute the amplitudes of each wave + dimensionalize
   std::vector<double> amps(freqs.size());
   for (std::size_t i = 0; i < amps.size(); i++)
   {
      amps[i] = std::sqrt(2*powers[i]*op.scale_fac.value_or(1.0))
                  *base_flow_params.p;
   }

   // Compute the phases of each wave
   std::vector<double> phases(freqs.size());
   std::mt19937 phase_gen(op.phase_seed);
   std::uniform_real_distribution<double> phase_dist(0,2*M_PI);
   for (std::size_t i = 0; i < phases.size(); i++)
   {
      phases[i] = phase_dist(phase_gen);
   }

   // Compute the normalized directions of each wave
   std::vector<std::vector<double>> k_hats(freqs.size());
   std::visit(DirectionVisitor{k_hats}, op.dir_params);

   // Finally, assemble the individual Wave structs
   for (std::size_t i = 0; i < op.num_waves; i++)
   {
      const Wave w{amps[i], freqs[i], phases[i], op.speed, k_hats[i]};
      waves.emplace_back(w);
   }
}

void SourceVisitor::operator()
   (const Source::Params<WaveCSV> &op)
{
   std::ifstream is(op.file);
   if (!is.is_open())
   {
      throw std::invalid_argument("Wave CSV file not found.");
   }
   ReadWaves(is, waves);
}

AcousticField InitializeAcousticField(const ConfigInput &conf, 
                                       std::span<const double> coords,
                                       int dim)
{
   // Get relevant input metadata
   const BaseFlowParams &base_conf = conf.BaseFlow();
   const CompParams &comp_conf = conf.Comp();
   const std::vector<Source::ParamsVariant> &sources_conf = conf.Sources();

   // Initialize acoustic field
   AcousticField field(dim, coords, base_conf.p, base_conf.rho,
                        base_conf.U, base_conf.gamma, comp_conf.kernel);

   // Assemble vector of wave structs based on input source
   for (const Source::ParamsVariant &source : sources_conf)
   {
      std::visit(SourceVisitor{base_conf, field.Waves()}, source);
   }

   // Finalize the acoustic field initialization
   field.Finalize();

   return std::move(field);

}                                        

} // namespace app
} // namespace jabber
