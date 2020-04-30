/*
 * Numerical module that computes the long-range interaction integrals
 */
#include "longrange.h"
#include "background.h"
#include "common.h"

#define SQR(x)  ((x)*(x))  // square of a number
#define CUB(x)  ((x)*(x)*(x)) // cube of a number

/* Abscisae and weights for the 16-point rule integration */
const double x16[8] = {0.0950125098376374401853193,0.2816035507792589132304605,0.4580167776572273863424194,0.6178762444026437484466718,0.7554044083550030338951012,0.8656312023878317438804679,0.9445750230732325760779884,0.9894009349916499325961542};
const double w16[8] = {0.1894506104550684962853967,0.1826034150449235888667637,0.1691565193950025381893121,0.1495959888165767320815017,0.1246289712555338720524763,0.0951585116824927848099251,0.0622535239386478928628438,0.0271524594117540948517806};

double precision_longrange = 1e-5; // Precision when numerically computing the scalar field
/* rho_class = 8piG/3 * rho_physical / Mpc^-2 = 1374.3 * rho_physical/eV^4
   Alternatively, 
   1374.3 = 8*pi*G/3 / (hbar^3*c^7) * Mpc_over_m^2 * eV^4

   This factor will be stored here
*/
const double _eV4_to_rho_class = 8*M_PI*_G_/3 / (pow(_h_P_/(2*M_PI), 3) * pow(_c_, 7)) * pow(_Mpc_over_m_, 2) * pow(_eV_, 4);

inline double I1(double u, double mtilde, double T){
  /* Integrand for computing the scalar field:
   *  u^2/sqrt(u^2 + mtilde_over_T^2) 1/[exp(u) + 1]
   * u      --- Dummy integration variable
   * mtilde --- m_F + g*phi
   * T      --- Temperature
   */
  return SQR(u)/sqrt(SQR(u) + SQR(mtilde/T)) * 1/(exp(u)+1);
}

inline double I2(double u, double mtilde, double T){
  /* Integrand for computing the scalar field derivative
   *  u^2 (2u^2 + 3 mtilde^2) / (u^2 + mtilde_over_T^2)^(3/2) 1/[exp(u) + 1]
   * u      --- Dummy integration variable
   * mtilde --- m_F + g*phi
   * T      --- Temperature
   */
  return SQR(u) * (2*SQR(u) + 3*SQR(mtilde/T)) / CUB(sqrt(SQR(u) + SQR(mtilde/T))) * 
    1/(exp(u)+1);
}

inline double I3(double u, double mtilde, double T){
  /* Integrand for computing the scalar field derivative
   *  u^4 / (u^2 + mtilde_over_T^2)^(3/2) 1/[exp(u) + 1]
   * u      --- Dummy integration variable
   * mtilde --- m_F + g*phi
   * T      --- Temperature
   */
  return SQR(SQR(u)) / CUB(sqrt(SQR(u) + SQR(mtilde/T))) * 
    1/(exp(u)+1);
}

inline double I4(double u, double mtilde, double T){
  /* Integrand for computing the fermion pressure derivative
   *  u^4 / sqrt(u^2 + mtilde_over_T^2) 1/[exp(u) + 1]
   * u      --- Dummy integration variable
   * mtilde --- m_F + g*phi
   * T      --- Temperature
   */
  return SQR(SQR(u)) / sqrt(SQR(u) + SQR(mtilde/T)) *
    1/(exp(u)+1);
}

inline double I5(double u, double mtilde, double T){
  /* Integrand for computing the fermion energy density
   * u^2 sqrt(u^2 + mtilde_over_T^2) 1/[exp(u) + 1]
   * u      --- Dummy integration variable
   * mtilde --- m_F + g*phi
   * T      --- Temperature
   */
  return SQR(u) * sqrt(SQR(u) + SQR(mtilde/T)) *
    1/(exp(u)+1);
}

inline double I6(double u, double mtilde, double T){
  /* Integrand for computing the fermion pressure
   * u^4 / sqrt(u^2 + mtilde_over_T^2) 1/[exp(u) + 1]
   * u      --- Dummy integration variable
   * mtilde --- m_F + g*phi
   * T      --- Temperature
   */
  return SQR(SQR(u)) / sqrt(SQR(u) + SQR(mtilde/T)) *
    1/(exp(u)+1);
}  

double potentialPrime(double phi_M, double m_F, double g_over_M, int g_F, double T){
  /* Returns the function that, when equals to zero, gives the scalar field
   * (it's the derivative of the potential with respect to the field)
   *
   * phi_M    --- Scalar field times its mass 
   * m_F      --- Mass of the fermion sourcing the scalar field
   * g_over_M --- Coupling divided by the scalar mass
   * g_F      --- Number of fermionic degrees of freedom
   * T        --- Temperature
   */
  double mtilde = m_F + g_over_M*phi_M;
  double res = phi_M;

  /* Carry out the integral. We will integrate from 0 to 20 */
  double A = (20 - 0)/2.;
  double B = (20 + 0)/2.;

  double integral = 0;
  for(int i = 0; i<8; ++i)
    integral += A * w16[i] * (I1(A*x16[i] + B, mtilde, T) + 
			      I1(-A*x16[i] + B, mtilde, T));
  res += g_over_M * g_F * SQR(T)/(2*SQR(M_PI)) * mtilde * integral;

  return res;
}

double getphi_M(double m_F, double g_over_M, int g_F, double T){
  /*
   * Returns the scalar field times its mass solving the equations of motion
   *
   * m_F      --- Mass of the fermion sourcing the scalar field
   * g_over_M --- Coupling divided by the scalar mass
   * g_F      --- Number of fermionic degrees of freedom
   * T        --- Temperature
   */
  // Use Ridders' method for root finding
  double x2 = 0, x1=-m_F/g_over_M; // We know that -m <= g*phi <= 0

  ErrorMsg errmsg = "longrange.c:";
  int j,MAXIT=10000;
  double ans,fh,fl,fm,fnew,s,xh,xl,xm,xnew, xtol=1e-5;
  fl = potentialPrime(x1, m_F, g_over_M, g_F, T);
  fh = potentialPrime(x2, m_F, g_over_M, g_F, T);

  if ((fl > 0.0 && fh < 0.0) || (fl < 0.0 && fh > 0.0)) {
    xl=x1;
    xh=x2;
    ans=-1.11e11;
    for (j=1;j<=MAXIT;j++) {
      xm=0.5*(xl+xh);
      fm = potentialPrime(xm, m_F, g_over_M, g_F, T);
      s=sqrt(fm*fm-fl*fh);
      if (s == 0.0){
	return ans;
      }
      xnew=xm+(xm-xl)*((fl >= fh ? 1.0 : -1.0)*fm/s);
      if (fabs(xnew-ans) <= xtol*ans){
	return ans;
      }
      ans=xnew;
      fnew = potentialPrime(ans, m_F, g_over_M, g_F, T);
      if (fnew == 0.0){
	return ans;
      }
      if (NRSIGN(fm,fnew) != fm) {
        xl=xm;
        fl=fm;
        xh=ans;
        fh=fnew;
      } else if (NRSIGN(fl,fnew) != fl) {
        xh=ans;
        fh=fnew;
      } else if (NRSIGN(fh,fnew) != fh) {
        xl=ans;
        fl=fnew;
      } else{
	class_stop(errmsg,"zriddr failed");
      }
      if (fabs(xh-xl) <= fabs(xtol*xh)) {
	return ans;
      }
    }
    class_stop(errmsg,"zriddr exceed maximum iterations");
  }
  else {
    if (fl == 0.0) return x1;
    if (fh == 0.0) return x2;
    class_stop(errmsg,"root must be bracketed in zriddr.");
  }
  class_stop(errmsg,"Failure in int.");
}

double getdphiM_dloga(double m_F, double g_over_M, int g_F, double T, double phi_M){
  /*
   * Returns the derivative of (scalar field times its mass) wrt logarithm of scale factor
   *
   * m_F      --- Mass of the fermion sourcing the scalar field
   * g_over_M --- Coupling divided by the scalar mass
   * g_F      --- Number of fermionic degrees of freedom
   * T        --- Temperature
   * phi_M    --- Scalar field times mass
   */
  double mtilde = m_F + g_over_M*phi_M;
  
  /* Carry out the integral. We will integrate from 0 to 20 */
  double A = (20 - 0)/2.;
  double B = (20 + 0)/2.;

  double integral1 = 0;
  for(int i = 0; i<8; ++i)
    integral1 += A * w16[i] * (I2(A*x16[i] + B, mtilde, T) + 
			       I2(-A*x16[i] + B, mtilde, T));
  double integral2 = 0;
  for(int i = 0; i<8; ++i)
    integral2 += A * w16[i] * (I3(A*x16[i] + B, mtilde, T) + 
			       I3(-A*x16[i] + B, mtilde, T));
  
  double res = g_over_M * SQR(T)/(2*SQR(M_PI)) * g_F * mtilde * integral1;
  res /= 1 + SQR(g_over_M) * SQR(T)/(2*SQR(M_PI)) * g_F * mtilde * integral2;

  return res;    
}

double getdPF_dloga(double m_F, double g_over_M, int g_F, double T, double phi_M){
  /*
   * Returns the derivative of the fermion pressure wrt logarithm of scale factor
   * m_F      --- Mass of the fermion sourcing the scalar field
   * g_over_M --- Coupling divided by the scalar mass
   * g_F      --- Number of fermionic degrees of freedom
   * T        --- Temperature
   * phi_M    --- Scalar field times mass
   */
  double mtilde = m_F + g_over_M*phi_M;

  /* Carry out the integral. We will integrate from 0 to 20 */
  double A = (20 - 0)/2.;
  double B = (20 + 0)/2.;

  double integral1 = 0;
  for(int i = 0; i<8; ++i)
    integral1 += A * w16[i] * (I3(A*x16[i] + B, mtilde, T) + 
			       I3(-A*x16[i] + B, mtilde, T));
  integral1 /= 3.0;

  double integral2 = 0;
  for(int i = 0; i<8; ++i)
    integral2 += A * w16[i] * (I4(A*x16[i] + B, mtilde, T) + 
			       I4(-A*x16[i] + B, mtilde, T));
  integral2 *= 4.0/3.0;

  double res = - getdphiM_dloga(m_F, g_over_M, g_F, T, phi_M) *
    g_over_M * SQR(T)/(2*SQR(M_PI)) * g_F * mtilde * integral1;

  res += -SQR(SQR(T))/(2*SQR(M_PI)) * g_F * integral2 -
    SQR(T)/(2*SQR(M_PI)) * SQR(mtilde) * g_F * integral1;

  return res;
}

double getrhoF(double m_F, double g_over_M, int g_F, double T, double phi_M){
  /*
   * Returns the fermion energy density
   * m_F      --- Mass of the fermion sourcing the scalar field
   * g_over_M --- Coupling divided by the scalar mass
   * g_F      --- Number of fermionic degrees of freedom
   * T        --- Temperature
   * phi_M    --- Scalar field times mass
   */
  double mtilde = m_F + g_over_M*phi_M;

  /* Carry out the integral. We will integrate from 0 to 20 */
  double A = (20 - 0)/2.;
  double B = (20 + 0)/2.;

  double integral = 0;
  for(int i = 0; i<8; ++i)
    integral += A * w16[i] * (I5(A*x16[i] + B, mtilde, T) + 
			      I5(-A*x16[i] + B, mtilde, T));
  return SQR(SQR(T))/(2*SQR(M_PI)) * g_F * integral;
}

double getPF(double m_F, double g_over_M, int g_F, double T, double phi_M){
  /*
   * Returns the fermion pressure
   * m_F      --- Mass of the fermion sourcing the scalar field
   * g_over_M --- Coupling divided by the scalar mass
   * g_F      --- Number of fermionic degrees of freedom
   * T        --- Temperature
   * phi_M    --- Scalar field times mass
   */
  double mtilde = m_F + g_over_M*phi_M;

  /* Carry out the integral. We will integrate from 0 to 20 */
  double A = (20 - 0)/2.;
  double B = (20 + 0)/2.;

  double integral = 0;
  for(int i = 0; i<8; ++i)
    integral += A * w16[i] * (I6(A*x16[i] + B, mtilde, T) + 
			      I6(-A*x16[i] + B, mtilde, T));
  return SQR(SQR(T))/(2*SQR(M_PI)) * g_F * integral / 3.0;
}
