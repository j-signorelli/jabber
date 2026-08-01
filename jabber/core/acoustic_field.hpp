#ifndef JABBER_ACOUSTIC_FIELD
#define JABBER_ACOUSTIC_FIELD

#include "base_flow.hpp"

#include <vector>
#include <span>
#include <iostream>
#include <cstdint>

namespace jabber
{

/**
 * @defgroup env_group Disturbance Environment Definition
 * @{
 *
 * @brief Modeling of an acoustically-disturbed freestream.
 *
 * @details
 *
 * ## Mathematical Formulation
 * Acoustic waves within tunnel freestreams are governed by
 * the linearized acoustic wave equation in a moving media,
 *
 * \f[
 * \frac{1}{c_\infty^2}\frac{D^2p'}{Dt^2} - \nabla^2p'
 * =\frac{1}{c_\infty^2}\left(\frac{\partial}{\partial t}
 *     + \vec{U}_\infty\cdot\nabla\right)^2p' - \nabla^2p'=0,
 * \f]
 *
 * where \f$p'\f$ is the pressure perturbation, \f$\vec{U}_\infty\f$
 * is the freestream velocity, and \f$c_\infty\f$ is the freestream
 * speed-of-sound.
 *
 * We assume **planar pressure waves** by the ansatz
 *
 * \f[
 * p'(\vec{x},t)=\sum_{j=1}^{N_w}p'_j\cos\left(\vec{k}_j\cdot\vec{x}
 *                - \omega_jt + \psi_j\right),
 * \f]
 *
 * where \f$p'_j\f$ is the provided pressure wave amplitude,
 * \f$\vec{k}_j\f$ is the wavenumber vector whose direction
 * \f$\hat{k}_j\f$ is provided, \f$\omega_j=2\pi f_j\f$ is
 * the provided angular frequency, \f$\psi_j\f$ is the provided
 * random phase, and \f$N_w\f$ are the number of planar pressure waves.
 * Plugging this into the linearized acoustic wave
 * equation, one can solve for the magnitude of the wavenumber
 * vector that satisfies the equation to be
 *
 * \f[
 * ||\vec{k}_j||=\frac{\omega_j}
 *                {\vec{U}_\infty\cdot\hat{k}_j \pm c_\infty},
 * \f]
 *
 * where \f$\pm\f$ corresponds to either a fast or slove acoustic wave
 * respectively.
 *
 * We can then derive the conserved variables from here. It can be shown that
 * the linearized isentropic equation of state is
 * \f$\rho'=(1/c_\infty^2)p'\f$. This allows us to compute freestream density
 * as
 *
 * \f[
 * \rho(\vec{x},t)=\rho_\infty + \rho'=
 * \rho_\infty+\frac{1}{c_\infty^2}
 *    \sum_{j=1}^{N_w}p'_j\cos\left(\vec{k}_j\cdot\vec{x} +
 *                                     \psi_j-\omega_jt\right).
 * \f]
 *
 * For momentum, we start with the linearized isentropic Euler equation,
 * \f[
 * \rho_\infty\frac{D\vec{u}'}{Dt}=
 * \rho_\infty\left(\frac{\partial \vec{u}}{\partial t}
 *                   +\vec{U}_\infty\cdot\nabla\right)\vec{u}'=-\nabla p'.
 * \f]
 *
 * The same traveling wave ansatz velocity perturbations can be written
 * as for pressure,
 *
 * \f[
 * \vec{u}'(\vec{x},t)=\sum_{j=1}^{N_w}\vec{u}'_j
 *    \cos\left(\vec{k}_j\cdot\vec{x} + \psi_j-\omega_jt\right).
 * \f]
 *
 * Plugging the pressure and velocity perturbation ansatz into the linearized
 * isentropic Euler equation, the velocity perturbation coefficients
 * \f$\vec{u}_j\f$ can be derived and the momentum written as
 * \f[
 * \rho\vec{u}=(\rho_\infty+\rho')(\vec{U}_\infty+\vec{u}')=
 * (\rho_\infty+\rho')\left(\vec{U}_\infty
 *    + \frac{(\pm1)}{\rho_\infty c_\infty}
 * \sum_{j=0}^{N_w}\hat{k}_jp'_j\cos\left(\vec{k}\cdot\vec{x}_j +
 *                                           \psi_j-\omega_jt\right)\right),
 * \f]
 * where the \f$\pm\f$ again corresponds to either a fast or slow acoustic wave
 * respectively. Lastly, energy can simply computed by
 * \f[
 * \rho E=\frac{p}{\gamma-1}+\frac{1}{2}\rho||\vec{u}||^2.
 * \f]
 *
 *
 * ## Wave Definition
 * From the previous section, we see that a single wave is fully defined by:
 *
 *    - \f$p'_j\f$: the planar wave amplitude,
 *    - \f$\hat{k}_j\f$: the normalized wavenumber vector
 *                         orientation/direction,
 *    - \f$f_j\f$: the wave temporal frequency,
 *    - \f$\psi_j\f$: the random phase, and
 *    - \f$\pm\f$: the wave speed (\f$+\f$ for fast, \f$-\f$ for slow).
 *
 * In jabber, this is represented by the \ref Wave struct.
 *
 * ## Acoustic Environment Definition & Evaluation
 * For a large \f$N_w\f$, the computation of the acoustically-perturbed
 * freestream at every gridpoint, every timestep can be computationally
 * expensive. For this, jabber provides high-performance \ref kernels_group
 * to accomplish this. The kernels require the data to be formatted a
 * particular way  (depending on which is used). To simplify the
 * translation from a set of acoustic waves defined by \ref Wave the data
 * required by the kernels structures, the class \ref AcousticField was
 * created.
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
void WriteWaves(const std::span<const Wave> &waves, std::ostream &out);

/**
 * @brief Read in Wave structs from CSV file, as outputted by
 * \ref WriteWaves().
 *
 * @details Parsed waves are appended to \p waves.
 */
void ReadWaves(std::istream &in, std::vector<Wave> &waves);

/**
 * @brief Compute the wavenumber vector magnitude for given freestream and
 * wave properties.
 *
 * @details Evaluates the dispersion relation
 *
 * \f[
 *    ||\vec{k}||=\frac{2\pi f}{\hat{k}\cdot\vec{U}_\infty \pm c_\infty}.
 * \f]
 *
 * @param U_infty     Freestream velocity vector.
 * @param c_infty     Freestream speed-of-sound,
 *                    \f$c_\infty=\sqrt{\gamma p_\infty/\rho}
 *                              =\sqrt{\gamma R T_\infty}\f$.
 * @param wave_freq   Wave frequency.
 * @param wave_k_hat  Wave wavenumber orientation. **Note that
 *                   `wave_k_hat.size() == U_infty.size()`**.
 * @param wave_speed  Wave speed. 'S' for slow, 'F' for fast.
 */
double ComputeWavenumber(const std::span<const double> &U_infty,
                         const double &c_infty,
                         const double &wave_freq,
                         const std::span<const double> &wave_k_hat,
                         const char &wave_speed);

/**
 * @brief Compute the wavenumber vector magnitude for a given
 * \ref Wave. See \ref ComputeWavenumber().
 *
 * @details **Note that `wave.k_hat.size() == U_infty.size()`**.
 */
inline double ComputeWavenumber(std::span<const double> U_infty,
                                const double &c_infty,
                                const Wave &wave)
{
   return ComputeWavenumber(U_infty, c_infty, wave.frequency, wave.k_hat,
                            wave.speed);
}

/**
 * @brief Class for specifying and computing a broadband-spectrum acoustic
 * field onto a provided grid and base flow.
 *
 * @details This class serves to simplify + streamline the preparation
 * of kernel argument data structures from a provided set of \ref Wave
 * structs.
 *
 */
class AcousticField
{
public:

   /**
    * @brief Kernel type to use in \ref Compute "Compute()".
    *
    * @details See \ref kernels_group for more information.
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

   /// Base flow properties.
   const BaseFlow base_flow_;

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
       * @brief Density series coefficients, \f$\frac{1}{c^2_\infty}p'_j\f$.
       *
       * @details Size is \ref NumWaves().
       */
      std::vector<double> rho_coeffs;

      /**
       * @brief Momentum series coefficients,
       * \f$\frac{1}{\rho_\infty c_\infty}(\pm 1)\hat{k_j}p'_j\f$.
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
       * @brief \f$\vec{k}\cdot x+\psi\f$ term computed for all waves at all
       * points.
       *
       * @details Size is \ref NumWaves() x \ref NumPoints(). Ordering depends
       * on \ref kernel_.
       */
      std::vector<double> k_dot_x_p_psi;

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
    * @brief Construct a new Acoustic Field object.
    * 
    * @param dim        @copybrief dim_
    * @param coords     Mesh coordinates to compute acoustic forcing on, in
    *                   XYZ XYZ ordering.
    * @param base_flow  @copybrief base_flow_
    * @param kernel     @copybrief kernel_
    */
   AcousticField(int dim, const std::span<const double> &coords,
                  const BaseFlow &base_flow, Kernel kernel=Kernel::GridPoint);
   /**
    * @brief Construct a new AcousticField object.
    *
    * @param dim        @copybrief dim_
    * @param coords     Mesh coordinates to compute acoustic forcing on, in
    *                   XYZ XYZ ordering.
    * @param p_infty    @copybrief BaseFlow::p
    * @param rho_infty  @copybrief BaseFlow::rho
    * @param u_infty    @copybrief BaseFlow::u
    * @param gamma      @copybrief BaseFlow::gamma
    * @param kernel     @copybrief kernel_
    */
   AcousticField(int dim, const std::span<const double> &coords,
                 double p_infty, double rho_infty,
                 const std::span<const double> &u_infty, double gamma,
                 Kernel kernel=Kernel::GridPoint)
   : AcousticField(dim, coords, BaseFlow(gamma, p_infty, rho_infty, u_infty),
                     kernel)
   {}

   /// Get the spatial dimension.
   int Dim() const
   {
      return dim_;
   }

   /// Get the number of points/coordinates associated with this field.
   std::size_t NumPoints() const
   {
      return num_pts_;
   }

   /// Get the base flow struct.
   const BaseFlow& GetBaseFlow() const
   {
      return base_flow_;
   }

   /// Get the number of waves.
   int NumWaves() const
   {
      return waves_.size();
   }

   /// Set the number of waves.
   void SetNumWaves(int num_waves)
   {
      waves_.resize(num_waves);
   }

   /**
    * @brief Reserve a number of waves to be included in the field.
    *
    * @details Calling this prior to any \ref AddWave() can improve
    * initialization efficiency.
    */
   void ReserveNumWaves(int res_waves)
   {
      waves_.reserve(res_waves);
   }

   /// Add a Wave to the acoustic field.
   void AddWave(const Wave &w)
   {
      waves_.push_back(w);
   }

   /// Get reference to Wave vector.
   std::vector<Wave> &Waves()
   {
      return waves_;
   }

   /// Get const reference to Wave vector.
   const std::vector<Wave>& Waves() const
   {
      return waves_;
   }

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
   std::span<double> Density()
   {
      return rho_;
   }

   /**
    * @brief Get const span of computed flow densities.
    *
    * @warning This should only be called after \ref Compute().
    */
   std::span<const double> Density() const
   {
      return rho_;
   }

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
   std::span<double> Energy()
   {
      return rhoE_;
   }

   /**
    * @brief Get const span of computed flow energy.
    *
    * @warning This should only be called after \ref Compute().
    */
   std::span<const double> Energy() const
   {
      return rhoE_;
   }

};

/// @}
// end of env_group

} // namespace jabber

#endif // JABBER_ACOUSTIC_FIELD
