#ifndef P3_NUCLEATE_ICE_LP05_IMPL_HPP
#define P3_NUCLEATE_ICE_LP05_IMPL_HPP

#include "p3_functions.hpp" // for ETI only but harmless for GPU
#include "p3_subgrid_variance_scaling_impl.hpp"
#include "mam4xx/nucleate_ice.hpp"
#include "mam4xx/wv_sat_methods.hpp"
#include "share/util/eamxx_common_physics_functions.hpp"

namespace scream {
namespace p3 {

using PF = scream::PhysicsFunctions<DefaultDevice>;

/*
 * Implementation of p3 contact and immersion freezing droplets function.
 * Clients should NOT #include this file, but include p3_functions.hpp instead.
 */

template <typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>
::nucleate_ice_lp05(
  // Inputs (what we need to call nucleati)
  const Spack& T_atm, const Spack& P_atm,
  const Spack& qv_atm, const Spack& omega_atm,
  const Spack& qc_incld, const Spack& rhoair_atm,
  // Output tendencies
  Spack& nihf, Spack& niimm, Spack& nidep, Spack& nimey,
  const P3Runtime& runtime_options,
  const Smask& context)
{
  constexpr Scalar qsmall   = C::QSMALL;
  constexpr Scalar T_rainfrz = C::T_rainfrz;
  constexpr Scalar T_zerodegc = C::T_zerodegc;
  constexpr Scalar CONS5    = C::CONS5;
  constexpr Scalar CONS6    = C::CONS6;
  const Scalar immersion_freezing_exponent =
      runtime_options.immersion_freezing_exponent;

  auto updraft_w_atm = PF::calculate_vertical_velocity(omega_atm, rhoair_atm);
  for (int p = 0; p < Spack::n; ++p) {
      Real tair = T_atm[p];
      Real pmid = P_atm[p];
      Real qv_scalar = qv_atm[p];
      Real wbar = updraft_w_atm[p];
      Real es = 0;
      Real qs = 0;
      Real cldn = 1;
      Real rhoair = rhoair_atm[p];
      mam4::wv_sat_methods::wv_sat_qsat_water(tair, pmid, es, qs);
      const Real relhum = qv_scalar / qs;
      Real so4_num = 10;
      Real dst3_num = 1e-2;
      Real nuci = 0;
      Real subgrid = 1;
      mam4::nucleate_ice::nucleati(
          wbar, tair, pmid, relhum, cldn, rhoair, so4_num, dst3_num, subgrid,
          nuci, nihf[p], niimm[p], nidep[p], nimey[p]
      );
  }
}

} // namespace p3
} // namespace scream

#endif
