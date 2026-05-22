#ifndef JABBER_APP_PARAMS
#define JABBER_APP_PARAMS

#include "../core/core.hpp"

#include <cstdint>
#include <vector>
#include <string>
#include <variant>
#include <array>
#include <optional>

namespace jabber
{
namespace app
{
/**
 * @defgroup params_group All Parameters/Settings
 * @{
 * @brief All input parameter options and structs.
 * 
 * @details To support settings that accept multiple different types of input,
 * these "option-specific" settings are enclosed in a struct with the following
 * schema:
 * 
 *    1. An enum `Option` of underlying type `std::uint8_t` with its final
 *       enumerator being `Size`,
 *    2. A string array `kNames` of size `Size`, associating a name with
 *       each `Option`,
 *    3. A primary template for a struct called `Params`, accepting an
 *       `Option`, with explicitly-specialized versions for each `Option` as
 *       needed, and
 *    4. A `std::variant` called `ParamsVar` containing a `Params` type for all
 *       `Option`s.
 * 
 * This schema was carefully crafted to reap the following benefits:
 * 
 *    - Enclosing all relevant information for an option-specific setting into a 
 *      single struct,
 *    - Readability, including explicit input variable names and the ability to 
 *      use `using enum ...;` inside each struct, and
 *    - Flexibility to define a `concept` for it if desired in the future.
 * 
 * All other settings that do not have >1 option are defined in standalone structs.
 */

// ----------------------------------------------------------------------------
/// Base flow parameters.
struct BaseFlowParams
{
   /// Density.
   double rho;

   /// Pressure.
   double p;

   /// Velocity vector.
   std::vector<double> U;

   /// Specific heat ratio.
   double gamma;
};

// ----------------------------------------------------------------------------
/// All options and parameters associated with input x,y data.
struct InputXY
{
   enum class Option : std::uint8_t
   {
      /// Provide x,y data directly in config file.
      Here,

      /// Read in x,y data from CSV file.
      FromCSV,

      /// Number of InputXY options.
      Size
   };

   using enum Option;

   static constexpr std::array<std::string_view, 
                     static_cast<std::size_t>(Size)>
   kNames = 
   {
      "Here",       // Here
      "FromCSV",    // FromCSV
   };


   template<Option V>
   struct Params;

   template<>
   struct Params<Here>
   {
      /// Input x's.
      std::vector<double> x;

      /// Input y's.
      std::vector<double> y;
   };

   template<>
   struct Params<FromCSV>
   {
      /// CSV file address. First column are x's, second are y's. No header.
      std::string file;
   };

   using ParamsVariant = std::variant<Params<Here>,Params<FromCSV>>;
};

// ----------------------------------------------------------------------------
/**
 * @brief All options and parameters associated with R->R continuous
 * function definition.
 */
struct FunctionType
{
   enum class Option : std::uint8_t
   {
      /// Piecewise linear fit.
      PiecewiseLinear,

      /// Piecewise lo10-log10 fit (linear on log10-log10 scale).
      PiecewiseLogLog,

      /// Number of FunctionType options.
      Size
   };

   using enum Option;

   static constexpr std::array<std::string_view, 
                        static_cast<std::size_t>(Size)>
   kNames = 
   {
      "PiecewiseLinear",    // PiecewiseLinear
      "PiecewiseLogLog",    // PiecewiseLogLog
   };


   template<Option V>
   struct Params;

   template<>
   struct Params<PiecewiseLinear>
   {
      /// Input x,y data params.
      InputXY::ParamsVariant input_xy;
   };

   template<>
   struct Params<PiecewiseLogLog>
   {
      /// Input x,y data params.
      InputXY::ParamsVariant input_xy;
   };

   using ParamsVariant = std::variant<Params<PiecewiseLinear>,
                                       Params<PiecewiseLogLog>>;

};
// ----------------------------------------------------------------------------
/**
 * @brief All options associated with interval/frequency bin
 * definition.
 */
struct IntervalType
{
   using Option = Interval::Method;

   using enum Option;

   static constexpr std::array<std::string_view, 
                        static_cast<std::size_t>(Size)>
   kNames = 
   {
      "Midpoint",      // Interval::Method::Midpoint
      "MidpointLog",   // Interval::Method::MidpointLog10
   };

};

// ----------------------------------------------------------------------------
/**
 * @brief All options and parameters associated with 1D discretization, such
 * as frequencies.
 */
struct DiscMethod
{
   enum class Option : std::uint8_t
   {
      /**
       * @brief Uniformly sample across range.
       * @warning If discretizing frequencies, harmonic interaction
       * may occur.
       */
      Uniform,

      /**
       * @brief Uniformly sample across range in log10 scaling.
       * @warning If discretizing frequencies, harmonic interaction
       * may occur.
       */
      UniformLog,

      /// Random sampling with a uniform distribution.
      Random,

      /// Random sampling with a uniform distribution on log10 scale.
      RandomLog,

      /// Number of DiscMethod options
      Size
   };
   
   using enum Option;

   static constexpr std::array<std::string_view, 
                        static_cast<std::size_t>(Size)>
   kNames = 
   {
      "Uniform",      // Uniform
      "UniformLog",   // UniformLog
      "Random",       // Random
      "RandomLog",    // RandomLog
   };


   template<Option V>
   struct Params
   {
      // Default struct has zero params.
   };

   template<>
   struct Params<Random>
   {
      /// Seed to use in randomization.
      int seed;
   };

   template<>
   struct Params<RandomLog>
   {
      /// Seed to use in randomization.
      int seed;
   };

   using ParamsVariant = std::variant<Params<Uniform>,Params<UniformLog>,
                                       Params<Random>,Params<RandomLog>>;
};

// ----------------------------------------------------------------------------
/**
 * @brief All options and parameters associated with wavenumber vector 
 * direction specification for a set of waves as in 
 * \p Source::Params<Source::Option::PSD>.
 */
struct Direction
{
   enum class Option : std::uint8_t
   {
      /// Single direction for all waves.
      Single,

      /// Random angle in XY-plane from x-axis for each wave.
      RandomXYAngle,

      /// Number of Direction options.
      Size
   };

   using enum Option;

   static constexpr std::array<std::string_view, 
                        static_cast<std::size_t>(Size)>
   kNames = 
   {
      "Single",           // Single
      "RandomXYAngle",    // RandomXYAngle
   };

   template<Option V>
   struct Params;

   template<>
   struct Params<Single>
   {
      /// Wavenumber vector direction, can be non-normalized.
      std::vector<double> direction;
   };

   template<>
   struct Params<RandomXYAngle>
   {
      /**
       * @brief Min angle of wavenumber vector from x-axis, in XY-plane
       * (CCW+, CW-) (in degrees).
       */
      double min_angle;

      /**
       * @brief Max angle of wavenumber vector from x-axis, in XY-plane 
       * (CCW+, CW-) (in degrees).
       */
      double max_angle;

      /// Seed to use in randomization.
      int seed;
   };

   using ParamsVariant = std::variant<Params<Single>, Params<RandomXYAngle>>;

};

// ----------------------------------------------------------------------------
/**
 * @brief All options and parameters associated with transfer function
 * selection.
 */
struct TransferFunction
{
      /// Transfer function options.
   enum class Option : std::uint8_t
   {
      /**
       * @brief Analytical low-frequency limit transfer function, from Chaudhry
       * & Candler, 2017.
       */
      LowFrequencyLimit,

      /// Provide a transfer function.
      Input,

      /**
       * @brief Extrapolate from an approximate fit of the collapsed/normalized
       * flow-normal transfer function in Chaudhry & Candler, 2017.
       * 
       */
      FlowNormalFit,

      /// Number of TransferFunction options.
      Size
   };

   using enum Option;

   static constexpr std::array<std::string_view, 
                              static_cast<std::size_t>(Size)>
   kNames = 
   {
      "LowFrequencyLimit",  // LowFrequencyLimit
      "Input",              // Input
      "FlowNormalFit"       // FlowNormalFit
   };

   template<Option V>
   struct Params
   {
      // Default struct has zero params
   };

   template<>
   struct Params<Input>
   {
      /// Transfer function representation, (f, V^2).
      FunctionType::ParamsVariant input_tf;
   };

   template<>
   struct Params<FlowNormalFit>
   {
      /// Shock standoff distance from pitot probe.
      double shock_standoff_dist;
   };

   using ParamsVariant = std::variant<Params<LowFrequencyLimit>, Params<Input>,
                                       Params<FlowNormalFit>>;
};

// ----------------------------------------------------------------------------
/**
 * @brief All options and parameters associated with acoustic forcing source
 * definition.
 */
struct Source
{
   enum class Option : std::uint8_t
   {
      /// Single acoustic wave.
      SingleWave,

      /// Spectrum of N acoustic waves.
      WaveSpectrum,

      /// Power spectral density (PSD).
      PSD,

      /// Read in CSV file of Wave data (output from \ref WriteWaves()).
      WaveCSV,

      /// Number of SourceOptions.
      Size
   };
   
   using enum Option;

   static constexpr std::array<std::string_view, 
                              static_cast<std::size_t>(Size)>
   kNames = 
   {
      "SingleWave",      // SingleWave
      "WaveSpectrum",    // WaveSpectrum
      "PSD",             // PSD
      "WaveCSV"          // WaveCSV
   };

   template<Option V>
   struct Params;

   template<>
   struct Params<SingleWave>
   {
      /// Wave amplitude.
      double amp;

      /// Wave frequency (not angular).
      double freq;

      /// Phase, in deg.
      double phase;

      /// Wavenumber vector direction, can be non-normalized.
      std::vector<double> direction;

      /// Wave speed ('S' or 'F').
      char speed;
   };

   template<>
   struct Params<WaveSpectrum>
   {
      /// Wave amplitudes.
      std::vector<double> amps;

      /// Wave frequencies (not angular).
      std::vector<double> freqs;

      /// Phases, in deg.
      std::vector<double> phases;

      /// Wavenumber vector directions, can be non-normalized.
      std::vector<std::vector<double>> directions;

      /// Wave speeds ('S' or 'F').
      std::vector<char> speeds;
   };

   template<>
   struct Params<PSD>
   {

      /// PSD function representation (f, PSD).
      FunctionType::ParamsVariant input_psd;

      /// (For PSD unit V^2/Hz) Scaling factor to multiply power V^2 by.
      std::optional<double> scale_fac;

      /// Minimum wave frequency in discrete frequency selection range.
      double min_disc_freq;

      /// Maximum wave frequency in discrete frequency selection range.
      double max_disc_freq;

      /// Number of waves to discretize PSD to.
      std::size_t num_waves;

      /// Interval method to use for frequency bins.
      Interval::Method int_method;

      /// Frequency discretization method parameters.
      DiscMethod::ParamsVariant disc_params;

      /// Wavenumber direction parameters.
      Direction::ParamsVariant dir_params;

      /// Seed to use for wave phase randomization.
      int phase_seed;

      /// Wave speeds to use.
      char speed;

      /// Transfer function parameters. None if not set.
      std::optional<TransferFunction::ParamsVariant> tf_params;
   };

   template<>
   struct Params<WaveCSV>
   {
      /// Wave CSV file (output from \ref WriteWaves()).
      std::string file;
   };

   using ParamsVariant = std::variant<Params<SingleWave>,Params<WaveSpectrum>,
                                       Params<PSD>,Params<WaveCSV>>;
};


// ----------------------------------------------------------------------------
/**
 * @brief All options associated with kernel type selection.
 */
struct KernelType
{
   using Option = AcousticField::Kernel;

   using enum Option;

   static constexpr std::array<std::string_view, 
                  static_cast<std::size_t>(Size)>
   kNames = 
   {
      "GridPoint",      // AcousticField::Kernel::GridPoint
      "Wave",           // AcousticField::Kernel::Wave
   };

};

// ----------------------------------------------------------------------------
/// Struct for computation parameters.
struct CompParams
{  
   /// Initial time.
   double t0;

   /// Kernel type.
   AcousticField::Kernel kernel;
};

// ----------------------------------------------------------------------------

/// Struct for preCICE parameters.
struct PreciceParams
{

   /// Jabber participant name.
   std::string participant_name;

   /// Address to preCICE config file.
   std::string config_file;

   /// Name of mesh to get coordinates from for computation onto.
   std::string fluid_mesh_name;

   /**
    * @brief Mesh access region, defined according to 
    * precice::Participant::setMeshAccessRegion()
    */
   std::vector<double> mesh_access_region;
};

// ----------------------------------------------------------------------------
/// @}
// end of params_group

} // namespace app

} // namespace jabber

#endif // JABBER_APP_PARAMS
