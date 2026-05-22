#ifndef JABBER_APP_COMMON
#define JABBER_APP_COMMON

#include "config_input.hpp"

#include "../core/core.hpp"
#include <iostream>
#include <memory>

namespace jabber
{
namespace app
{

constexpr static std::string_view LINE = 
"-----------------------------------------------------------------------";

/// Print Jabber banner.
void PrintBanner(std::ostream &out);

/// Normalize the provided vector data
void Normalize(std::span<const double> vec, std::span<double> norm_vec);

/**
 * @defgroup pproc_group Parameter Processing
 * @{
 * 
 */

/**
 * @brief All visitor options for each InputXY::Params, for initializing 
 * \ref x and \ref y vectors.
 * 
 */
struct InputXYVisitor
{  
   using enum InputXY::Option;

   /// Reference to x-data vector to set.
   std::vector<double> &x;

   /// Reference to y-data vector to set.
   std::vector<double> &y;

   void operator() (const InputXY::Params<Here> &op);
   void operator() (const InputXY::Params<FromCSV> &op);
};

/**
 * @brief All visitor options for each FunctionType::Params, for initializing a
 * \ref Function or \ref BasePSD type.
 * 
 */
struct FunctionTypeVisitor
{
   using enum FunctionType::Option;

   /// Function to initialize
   std::variant<std::unique_ptr<Function>*,
                std::unique_ptr<BasePSD>*> T_ptr_ptr_var;

   void operator() (const FunctionType::Params<PiecewiseLinear> &op);
   void operator() (const FunctionType::Params<PiecewiseLogLog> &op);
};

/**
 * @brief All visitor options for each DiscMethod::Params, for initializing
 * a discretized frequency range, \p freqs.
 * 
 */
struct DiscMethodVisitor
{  
   using enum DiscMethod::Option;

   /// Minimum frequency bound.
   const double &min_freq;

   /// Maximum frequency bound.
   const double &max_freq;

   /// **Sized** frequency vector to initialize.
   std::vector<double> &freqs;

   void operator() (const DiscMethod::Params<Uniform> &op);
   void operator() (const DiscMethod::Params<UniformLog> &op);
   void operator() (const DiscMethod::Params<Random> &op);
   void operator() (const DiscMethod::Params<RandomLog> &op);
};

/**
 * @brief All visitor options for each Direction::Params, for initializing
 * wavenumber vector directions in \p k_hats.
 * 
 */
struct DirectionVisitor
{
   using enum Direction::Option;

   /// **Sized** vector of direction vectors for each wave.
   std::vector<std::vector<double>> &k_hats;

   void operator() (const Direction::Params<Single> &op);
   void operator() (const Direction::Params<RandomXYAngle> &op);
};

/**
 * @brief All visitor options for each TransferFunction::Params, for
 * applying a transfer function to given \ref powers.
 * 
 */
struct TransferFunctionVisitor
{
   using enum TransferFunction::Option;

   /// Base flow parameters.
   const BaseFlowParams &base_flow_params;

   /// Array of frequencies to evaluate for.
   const std::vector<double> &freqs;

   /// Wave speed (for all).
   const char &speed;

   /// Array of powers to update; same size as \ref freqs.
   std::vector<double> &powers;

   void operator() (const TransferFunction::Params<LowFrequencyLimit> &op);
   void operator() (const TransferFunction::Params<Input> &op);
   void operator() (const TransferFunction::Params<FlowNormalFit> &op);
};

/**
 * @brief All visitor options for each Source::Params, for initializing 
 * \ref Wave's for each type and appending to \p waves.
 * 
 */
struct SourceVisitor
{
   using enum Source::Option;

   const BaseFlowParams &base_flow_params;

   /// Reference of wave vector to append Wave structs to.
   std::vector<Wave> &waves;
   
   void operator() (const Source::Params<SingleWave> &op);
   void operator() (const Source::Params<WaveSpectrum> &op);
   void operator() (const Source::Params<PSD> &op);
   void operator() (const Source::Params<WaveCSV> &op);
};

/// @}
// end of pproc_group

/**
 * @brief Initialize a \ref AcousticField object from user input and 
 * grid.
 * 
 * @param conf          Input config object.
 * @param coords        Mesh coordinates vector, in XYZ XYZ ordering.
 * @param dim           Spatial dimension of mesh.
 * @return AcousticField    Finalized acoustic field.
 */
AcousticField InitializeAcousticField(const ConfigInput &conf, 
                                                std::span<const double> coords,
                                                int dim);


/**
 * @brief Extremely simple function to get subspan of data from \p global 
 * to \p local for basic MPI-partitioning. Data must be ordered in a AoS 
 * format.
 * 
 * @tparam T         Type.
 * @param global     Global data set to get partition of. Note that all
 *                   ranks must have this defined and equal.
 * @param rank       Rank to get partition of.
 * @param vdim       Vector dimension of each data.
 * @param size       Number of ranks to partition across.
 * @param local      Output rank-local data sub-span to set based on args.
 */
template<typename T>
void GetRankPartition(std::span<const T> global, int vdim, int rank, int size, 
                        std::span<const T> &local)
{
   const std::size_t global_num_dat = global.size()/vdim;
   const std::size_t rank_num_dat = vdim*(global_num_dat/size +
                                          (rank < (global_num_dat % size) 
                                             ? 1 
                                             : 0));

   const std::size_t rank_offset = vdim*(rank*(global_num_dat/size) +
                                          (rank < (global_num_dat % size)
                                             ? rank 
                                             : (global_num_dat % size)));
   
   
   local = global.subspan(rank_offset, rank_num_dat);
}

} // namespace app
} // namespace jabber

#endif // JABBER_APP_COMMON
