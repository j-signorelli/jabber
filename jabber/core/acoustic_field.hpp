#ifndef JABBER_ACOUSTIC_FIELD
#define JABBER_ACOUSTIC_FIELD

#include <vector>
#include <span>
#include <iostream>
#include <cstdint>

namespace jabber
{

/**
 * @defgroup env_group Disturbance Environment Definition & Evaluation
 * @{
 * 
 */

/// Base acoustic wave definition.
struct Wave
{
   /// Wave amplitude, p'.
   double amplitude;

   /// Wave frequency, f.
   double frequency;

   /// Wave phase, φ, in radians.
   double phase;

   /// Wave speed ('S' or 'F').
   char speed;

   /// **Normalized** wavenumber vector direction.
   std::vector<double> k_hat;
};

/**
 * @brief Write span of \ref Wave structs to \p out as a CSV, with columns
 * [Amplitude, Frequency, Phase, Speed, k_hat].
 */
void WriteWaves(std::span<const Wave> waves, std::ostream &out);

/**
 * @brief Read in Wave structs from CSV file, as outputted by 
 * \ref WriteWaves().
 * 
 * @details Parsed waves are appended to \p waves.
 */
void ReadWaves(std::istream &in, std::vector<Wave> &waves);

/**
 * @brief Class for specifying and computing a broadband-spectrum acoustic
 * field onto a provided grid and base flow.
 * 
 * @details All acoustic wave data is stored in a Struct-of-Arrays style in
 * this class.
 * 
 */
class AcousticField
{
public:

   /**
    * @brief Kernel type to use in \ref Compute "Compute()".
    * 
    * @details See file kernels.hpp for more information.
    */
   enum class Kernel : std::uint8_t
   {
      /// Use grid-point axis in series summation inner-loop.
      GridPoint,

      /// Use wave axis in series summation inner-loop.
      Wave,
      
      /// Number of Kernel enumerators.
      Size,
   };

private:

   /// Spatial dimension.
   const int dim_;

   /// Number of points/coordinates of field.
   const std::size_t num_pts_;

   /// Base flow pressure.
   const double p_bar_;

   /// Base flow density.
   const double rho_bar_;

   /// Base flow velocity vector, of size \ref Dim().
   const std::vector<double> U_bar_;

   /// Base flow specific heat ratio, γ.
   const double gamma_;

   /// Base flow speed of sound
   const double c_bar_;

   /// Kernel type to use.
   const Kernel kernel_;
   
   /// SoA coordinates to compute waves on, [dim][node].
   std::vector<std::vector<double>> coords_;

   /// Array of all wave data (AoS).
   std::vector<Wave> waves_;

   /**
    * @brief Struct of kernel-prepped data structures, initialized in 
    * \ref Finalize().
    * 
    * @details Pre-assembling this before calls to \ref Compute() is required.
    */
   struct
   {
      /**
       * @brief Density series coefficients, \f$\frac{1}{\bar{c}^2}p'_j\f$.
       * 
       * @details Size is \ref NumWaves().
       */
      std::vector<double> rho_coeffs;

      /**
       * @brief Momentum series coefficients,
       * \f$\frac{1}{\bar{\rho}\bar{c}}(\pm 1)\hat{k_j}\f$.
       * 
       * @details Size is \ref Dim() x \ref NumWaves(). Ordered as [dim][wave].
       */
      std::vector<double> rhoV_coeffs;

      /**
       * @brief Energy series coefficients, \f$\frac{1}{(\gamma-1)}p'_j\f$.
       * 
       * @details Size is \ref NumWaves().
       */
      std::vector<double> rhoE_coeffs;

      /**
       * @brief Acoustic wave angular frequencies, \f$\omega=2\pi f\f$.
       * 
       * @details Size is \ref NumWaves().
       */
      std::vector<double> wave_omegas;

      /**
       * @brief \f$\vec{k}\cdot x+\phi\f$ term computed for all waves at all
       * points.
       * 
       * @details Size is \ref NumWaves() x \ref NumPoints(). Ordering depends
       * on \ref kernel_.
       */
      std::vector<double> k_dot_x_p_phi;

   } kernel_args_;

   /**
    * @brief Fluid density \f$\rho\f$, computed in \ref Compute().
    * 
    * @details Size is \ref NumPoints().
    */
   std::vector<double> rho_;

   /**
    * @brief Fluid momentum \f$\rho\vec{u}\f$, computed in \ref Compute().
    * 
    * @details Size is \ref NumPoints() * \ref Dim(), with data in XXX YYY 
    * ordering.
    */
   std::vector<double> rhoV_;

   /**
    * @brief Fluid energy \f$\rho E\f$, computed in \ref Compute().
    * 
    * @details Size is \ref NumPoints().
    */
   std::vector<double> rhoE_;

public:
   /**
    * @brief Construct a new AcousticField object.
    * 
    * @param dim        Spatial dimension of mesh.
    * @param coords     Mesh coordinates to compute acoustic forcing on, in
    *                   XYZ XYZ ordering.
    * @param p_bar      Base flow pressure.
    * @param rho_bar    Base flow density.
    * @param U_bar      Base flow velocity vector, of size \p dim.
    * @param gamma      Base flow specific heat ratio, γ.
    * @param kernel     Kernel type to use.
    */
   AcousticField(int dim, std::span<const double> coords,
                  double p_bar, double rho_bar,
                  const std::vector<double> U_bar, double gamma,
                  Kernel kernel=Kernel::GridPoint);

   /// Get the spatial dimension.
   int Dim() const { return dim_; }

   /// Get the number of points/coordinates associated with this field.
   std::size_t NumPoints() const { return num_pts_; }

   /// Get the base flow velocity vector.
   const std::vector<double>& BaseVelocity() const { return U_bar_; }
   
   /// Get the base flow pressure
   double BasePressure() const { return p_bar_; }

   /// Get the base flow density
   double BaseDensity() const { return rho_bar_; }

   /// Get the base flow specific heat ratio, γ.
   double Gamma() const { return gamma_; }

   /// Get the number of waves.
   int NumWaves() const { return waves_.size(); }

   /// Set the number of waves.
   void SetNumWaves(int num_waves) { waves_.resize(num_waves); }

   /**
    * @brief Reserve a number of waves to be included in the field.
    * 
    * @details Calling this prior to any \ref AddWave() can improve
    * initialization efficiency.
    */
   void ReserveNumWaves(int res_waves) { waves_.reserve(res_waves); }

   /// Add a Wave to the acoustic field.
   void AddWave(const Wave &w) { waves_.push_back(w); }

   /// Get reference to Wave vector.
   std::vector<Wave>& Waves() { return waves_; }

   /// Get const reference to Wave vector.
   const std::vector<Wave>& Waves() const { return waves_; }

   /**
    * @brief Finalize the acoustic field, to be called after specifying all
    * waves, before \ref Compute().
    * 
    * @details This function initializes \ref kernel_args_ based on
    * \ref kernel_ and allocates the flowfield solution \ref rho_, 
    * \ref rhoV_, and \ref rhoE_ vectors.
    */
   void Finalize();

   /**
    * @brief Compute the perturbed flowfield at time \p t, **after** calling
    * adding all wave data and calling \ref Finalize()
    * 
    * @warning \ref Finalize() must be called once prior to calls to this,
    * after adding all wave data.
    */
   void Compute(double t);
   
   /**
    * @brief Get span of computed flow densities.
    * 
    * @warning This should only be called after \ref Compute().
    */
   std::span<double> Density() { return rho_; }
   
   /**
    * @brief Get const span of computed flow densities.
    * 
    * @warning This should only be called after \ref Compute().
    */
   std::span<const double> Density() const { return rho_; }

   /**
    * @brief Get span of flow momentum across all components.
    * 
    * @warning This should only be called after \ref Compute().
    */
   std::span<double> Momentum()
   { 
      return std::span<double>(rhoV_);
   }
   
   /**
    * @brief Get const span of flow momentum across all components.
    * 
    * @warning This should only be called after \ref Compute().
    */
   std::span<const double> Momentum() const
   { 
      return std::span<const double>(rhoV_);
   }

   /**
    * @brief Get span of computed flow momentum for component \p comp.
    * 
    * @warning This should only be called after \ref Compute().
    */
   std::span<double> Momentum(int comp)
   { 
      return std::span<double>(rhoV_).subspan(num_pts_*comp, num_pts_);
   }
   
   /**
    * @brief Get const span of computed flow momentum for component \p comp.
    * 
    * @warning This should only be called after \ref Compute().
    */
   std::span<const double> Momentum(int comp) const
   { 
      return std::span<const double>(rhoV_).subspan(num_pts_*comp, num_pts_);
   }

   /**
    * @brief Get span of computed flow energy.
    * 
    * @warning This should only be called after \ref Compute().
    */
   std::span<double> Energy() { return rhoE_; }
   
   /**
    * @brief Get const span of computed flow energy.
    * 
    * @warning This should only be called after \ref Compute().
    */
   std::span<const double> Energy() const { return rhoE_; }

};

/// @}
// end of env_group

} // namespace jabber

#endif // JABBER_ACOUSTIC_FIELD
