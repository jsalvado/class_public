/* Numerical module that computes the long-range interaction integrals */
#ifndef __LONGRANGE__
#define __LONGRANGE__

extern const double _eV4_to_rho_class;

#ifdef __cplusplus
extern "C" {
#endif
double getphi_M(double m_F, double g_over_M, int g_F, double T);
// Returns the scalar field times its mass solving the equations of motion

double getdphiM_dloga(double m_F, double g_over_M, int g_F, double T, double phi_M);
// Returns the derivative of (scalar field times its mass) wrt logarithm of scale factor

double getdPF_dloga(double m_F, double g_over_M, int g_F, double T, double phi_M);
// Returns the derivative of the fermion pressure wrt logarithm of scale factor

double getrhoF(double m_F, double g_over_M, int g_F, double T, double phi_M);
// Returns the fermion energy density

double getPF(double m_F, double g_over_M, int g_F, double T, double phi_M);
// Returns the fermion pressure

#ifdef __cplusplus
}
#endif

#endif
