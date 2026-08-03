#ifndef JABBER_BASE_FLOW
#define JABBER_BASE_FLOW

#include <vector>
#include <span>
namespace jabber
{

struct BaseFlow
{
   /// Specific heat ratio.
   const double gamma;

   /// Freestream pressure.
   const double p;

   /// Freestream density.
   const double rho;

   /// Freestream velocity vector.
   const std::vector<double> U;

   /// Freestream speed-of-sound.
   const double c;

   /// Freestream Mach number.
   const double M;

   /**
    * @brief Construct a new BaseFlow object.
    *
    * @param gamma_     @copybrief gamma
    * @param p_         @copybrief p
    * @param rho_       @copybrief rho
    * @param U_         @copybrief U
    */
   BaseFlow(const double &gamma_, const double &p_, const double &rho_,
            const std::span<const double> &U_);
};

}

#endif // JABBER_BASE_FLOW
