/* Numerical module that computes the long-range interaction integrals */
#ifndef __LONGRANGE__
#define __LONGRANGE__

#include "background.h"
#include "common.h"
#include "input.h"

#define SQR(x)  ((x)*(x))  // square of a number
#define CUB(x)  ((x)*(x)*(x)) // cube of a number

extern const double _eV4_to_rho_class;
extern const double _J4_to_rho_class;
extern const double _Mpc_times_eV;
struct background_parameters_and_redshift {
  /* structures containing fixed input parameters (indices, ...) */
  struct background * pba;

  /* redshift */
  double z;
};

#ifdef __cplusplus
extern "C" {
#endif
  // Initialize the quadrature weights for the long-range interaction integrals
  int background_lrs_init(struct precision *ppr, struct background *pba);

  // Compute various integrals of the fermion distribution
  int background_lrs_momenta(double * qvec, double * wvec, int qsize,
			     double m_over_T0, double factor, double z,
			     double * rho, double * p,
			     double * drho_dm_over_T0, double * pseudo_p,
			     double * I1, double * I2);

  // Returns the effective fermion mass divided by its present day temperature
  double get_mT_over_T0_lrs(struct background * pba, double phi_M);

  // Returns the scalar field times its mass solving the equations of motion
  int get_phi_M_lrs(struct background * pba, double z, double *phi_M);

  // Returns the scale factor at adiabatic instability onset
  int instabilityOnset_lrs(struct background * pba, double * a_rel);
#ifdef __cplusplus
}
#endif

#endif
