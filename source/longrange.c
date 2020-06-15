/*
 * Numerical module that computes the long-range interaction integrals
 */
#include "longrange.h"

/* rho_class = 8piG/3 * rho_physical / Mpc^-2 = 1374.3 * rho_physical/eV^4
   Alternatively, 
   1374.3 = 8*pi*G/3 / (hbar^3*c^7) * Mpc_over_m^2 * eV^4

   This factor will be stored here
*/
const double _eV4_to_rho_class = 8*M_PI*_G_/3 / (pow(_h_P_/(2*_PI_), 3) * pow(_c_, 7)) * pow(_Mpc_over_m_, 2) * pow(_eV_, 4);
const double _J4_to_rho_class = 8*M_PI*_G_/3 / (pow(_h_P_/(2*_PI_), 3) * pow(_c_, 7)) * pow(_Mpc_over_m_, 2);
const double _Mpc_times_eV = _eV_/(_h_P_/(2*_PI_) * _c_) * _Mpc_over_m_; // Multiply by this to convert eV to Mpc^-1

/**
 * This is the routine where the distribution function f0(q) of the
 * fermion undergoing long-range interactions is specified
 *
 * @param pbadist Input:  structure containing all parameters defining f0(q)
 * @param q       Input:  momentum
 * @param f0      Output: phase-space distribution
 */
int background_lrs_distribution(
                                 void * pbadist,
                                 double q,
                                 double * f0
                                 ) {
  struct background * pba;
  struct background_parameters_for_distributions * pbadist_local;

  /** - extract from the input structure pbadist all the relevant information */
  pbadist_local = pbadist;          /* restore actual format of pbadist */
  pba = pbadist_local->pba;         /* extract the background structure from it */

  *f0 = pba->lrs_g_F * 1.0/CUB(2*_PI_) * 1./(exp(q)+1.);
  // Fermi-Dirac distribution with g_F degrees of freedom

  return _SUCCESS_;
}

/**
 * This is the routine where the derivative of the
 * distribution function f0(q) is specified
 *
 * @param pbadist         Input:  structure containing all parameters defining f0(q)
 * @param q               Input:  momentum
 * @param dlnf0_dlnq      Output: derivative of log[phase-space distribution] wrt log[q]
 */
int background_lrs_derivative(
                                 void * pbadist,
                                 double q,
                                 double * dlnf0_dlnq
                                 ) {
  struct background * pba;
  struct background_parameters_for_distributions * pbadist_local;

  /** - extract from the input structure pbadist all the relevant information */
  pbadist_local = pbadist;          /* restore actual format of pbadist */
  pba = pbadist_local->pba;         /* extract the background structure from it */

  double f0 = 1./(exp(q)+1.);
  // Fermi-Dirac distribution with g_F degrees of freedom

  //Avoid underflow in extreme tail:
  if (fabs(f0)==0.)
    *dlnf0_dlnq = -q;
  else
    *dlnf0_dlnq = -q * exp(q)/(exp(q)+1.);

  return _SUCCESS_;
}

/**
 * This function is only used for the purpose of finding optimal
 * quadrature weights. The logic is: if we can accurately convolve
 * f0(q) with this function, then we can convolve it accurately with
 * any other relevant function.
 *
 * @param pbadist Input:  structure containing all background parameters
 * @param q       Input:  momentum
 * @param test    Output: value of the test function test(q)
 */

int background_lrs_test_function(
                                  void * pbadist,
                                  double q,
                                  double * test
                                  ) {

  double c = 2.0/(3.0*_zeta3_);
  double d = 120.0/(7.0*SQR(SQR(_PI_)));
  double e = 2.0/(45.0*_zeta5_);

  /** Using a + bq creates problems for otherwise acceptable distributions
      which diverges as \f$ 1/r \f$ or \f$ 1/r^2 \f$ for \f$ r\to 0 \f$*/
  *test = CUB(2.0*_PI_)/6.0*(c*q*q-d*q*q*q-e*q*q*q*q);

  return _SUCCESS_;
}


/**
 * This function finds optimal quadrature weights for the
 * long-range interaction integrals
 *
 * @param ppr Input: precision structure
 * @param pba Input/Output: background structure
 */

int background_lrs_init(
                         struct precision *ppr,
                         struct background *pba
                         ) {

  int index_q, tolexp,row,status,filenum;
  double q, dlnf0_dlnq;
  struct background_parameters_for_distributions pbadist;
  FILE *psdfile;

  pbadist.pba = pba;

  pbadist.q = NULL;
  pbadist.tablesize = 0;

  /* Handle perturbation qsampling: */
  if (pba->lrs_quadrature_strategy==qm_auto){
    /** Automatic q-sampling for this species */
    class_alloc(pba->q_lrs,_QUADRATURE_MAX_*sizeof(double),pba->error_message);
    class_alloc(pba->w_lrs,_QUADRATURE_MAX_*sizeof(double),pba->error_message);

    class_call(get_qsampling(pba->q_lrs,
			     pba->w_lrs,
			     &(pba->q_size_lrs),
			     _QUADRATURE_MAX_,
			     ppr->tol_lrs,
			     pbadist.q,
			     pbadist.tablesize,
			     background_lrs_test_function,
			     background_lrs_distribution,
			     &pbadist,
			     pba->error_message),
	       pba->error_message,
	       pba->error_message);
    pba->q_lrs=realloc(pba->q_lrs,pba->q_size_lrs*sizeof(double));
    pba->w_lrs=realloc(pba->w_lrs,pba->q_size_lrs*sizeof(double));


    if (pba->background_verbose > 0)
      printf("lrs species sampled with %d points for purpose of perturbation integration\n",
	     pba->q_size_lrs);

    /* Handle background q_sampling: */
    class_alloc(pba->q_lrs_bg,_QUADRATURE_MAX_BG_*sizeof(double),pba->error_message);
    class_alloc(pba->w_lrs_bg,_QUADRATURE_MAX_BG_*sizeof(double),pba->error_message);

    class_call(get_qsampling(pba->q_lrs_bg,
			     pba->w_lrs_bg,
			     &(pba->q_size_lrs_bg),
			     _QUADRATURE_MAX_BG_,
			     ppr->tol_lrs_bg,
			     pbadist.q,
			     pbadist.tablesize,
			     background_lrs_test_function,
			     background_lrs_distribution,
			     &pbadist,
			     pba->error_message),
	       pba->error_message,
	       pba->error_message);


    pba->q_lrs_bg=realloc(pba->q_lrs_bg,pba->q_size_lrs_bg*sizeof(double));
    pba->w_lrs_bg=realloc(pba->w_lrs_bg,pba->q_size_lrs_bg*sizeof(double));

    /** - in verbose mode, inform user of number of sampled momenta
	for background quantities */
    if (pba->background_verbose > 0)
      printf("lrs species sampled with %d points for purpose of background integration\n",
	     pba->q_size_lrs_bg);
  }
  else{
    /** Manual q-sampling for this species. Same sampling used for both perturbation and background sampling, since this will usually be a high precision setting anyway */
    pba->q_size_lrs_bg = pba->lrs_input_q_size;
    pba->q_size_lrs = pba->lrs_input_q_size;
    class_alloc(pba->q_lrs_bg,pba->q_size_lrs_bg*sizeof(double),pba->error_message);
    class_alloc(pba->w_lrs_bg,pba->q_size_lrs_bg*sizeof(double),pba->error_message);
    class_alloc(pba->q_lrs,pba->q_size_lrs*sizeof(double),pba->error_message);
    class_alloc(pba->w_lrs,pba->q_size_lrs*sizeof(double),pba->error_message);
    class_call(get_qsampling_manual(pba->q_lrs,
				    pba->w_lrs,
				    pba->q_size_lrs,
				    pba->lrs_qmax,
				    pba->lrs_quadrature_strategy,
				    pbadist.q,
				    pbadist.tablesize,
				    background_lrs_distribution,
				    &pbadist,
				    pba->error_message),
	       pba->error_message,
	       pba->error_message);
    for (index_q=0; index_q<pba->q_size_lrs; index_q++) {
      pba->q_lrs_bg[index_q] = pba->q_lrs[index_q];
      pba->w_lrs_bg[index_q] = pba->w_lrs[index_q];
    }
    /** - in verbose mode, inform user of number of sampled momenta
	for background quantities */
    if (pba->background_verbose > 0)
      printf("lrs species sampled with %d points for purpose of background andperturbation integration using the manual method\n",
	     pba->q_size_lrs);
  }

  class_alloc(pba->dlnf0_dlnq_lrs,
	      pba->q_size_lrs*sizeof(double),
	      pba->error_message);


  for (index_q=0; index_q<pba->q_size_lrs; index_q++) {
    q = pba->q_lrs[index_q];
    class_call(background_lrs_derivative(&pbadist,q,&dlnf0_dlnq),
	       pba->error_message,pba->error_message);
    pba->dlnf0_dlnq_lrs[index_q] = dlnf0_dlnq;
  }

  /* rho_class = 8piG/3 * rho_physical / Mpc^-2
     The 4pi (kB T0)^4 will appear in front of integrals*/
  pba->factor_lrs=4*_PI_*SQR(SQR(pba->T_cmb*pba->lrs_T_F*_k_B_)) * _J4_to_rho_class;

  return _SUCCESS_;
}


/**
 * For a given long-range interacting species: given the quadrature weights, the mass
 * and the redshift, find background FERMION quantities by a quick weighted
 * sum over.  Input parameters passed as NULL pointers are not
 * evaluated for speed-up
 *
 * @param qvec             Input: sampled momenta
 * @param wvec             Input: quadrature weights
 * @param qsize            Input: number of momenta/weights
 * @param m_over_T0        Input: mass divided by *present day* fermion temperature
 * @param factor           Input: normalization factor for the p.s.d.
 * @param z                Input: redshift
 * @param rho              Output: energy density
 * @param p                Output: pressure
 * @param drho_dm_over_T0  Output: derivative used in next function
 * @param pseudo_p         Output: pseudo-pressure used in perturbation module for fluid approx
 * @param I1               Output: \integ d^3u m/T (2 eps^2 + (m/T)^2)/eps^3 f(u)
 *                         where eps^2 = u^2 + (m/T)^2
 * @param I2               Output: \integ d^3u u^2/eps^3 f(u)
 *                         where eps^2 = u^2 + (m/T)^2
 */

int background_lrs_momenta(
                            /* Only calculate for non-NULL pointers: */
                            double * qvec,
                            double * wvec,
                            int qsize,
                            double m_over_T0,
                            double factor,
                            double z,
                            double * rho, // density
                            double * p,   // pressure
                            double * drho_dm_over_T0,  // d rho / d(m_over_T0) used in next function
                            double * pseudo_p,  // pseudo-p used in ncdm fluid approx
			    double * I1,
			    double * I2
                            ) {

  int index_q;
  double epsilon;
  double q2;
  double factor2;
  /** Summary: */

  /** - rescale normalization at given redshift */
  factor2 = factor*SQR(SQR(1+z));

  /** - initialize quantities */
  if (rho!=NULL) *rho = 0.;
  if (p!=NULL) *p = 0.;
  if (drho_dm_over_T0!=NULL) *drho_dm_over_T0 = 0.;
  if (pseudo_p!=NULL) *pseudo_p = 0.;
  if (I1!=NULL) *I1 = 0.;
  if (I2!=NULL) *I2 = 0.;

  /** - loop over momenta */
  for (index_q=0; index_q<qsize; index_q++) {

    /* squared momentum */
    q2 = qvec[index_q]*qvec[index_q];

    /* energy */
    epsilon = sqrt(q2+m_over_T0*m_over_T0/(1.+z)/(1.+z));

    /* integrand of the various quantities */
    if (rho!=NULL) *rho += q2 * epsilon*wvec[index_q];
    if (p!=NULL) *p += q2 * q2/3./epsilon*wvec[index_q];
    if (drho_dm_over_T0!=NULL) *drho_dm_over_T0 += q2 * m_over_T0/SQR(1.+z)/epsilon*wvec[index_q];
    if (pseudo_p!=NULL) *pseudo_p += q2 * SQR(q2)/CUB(epsilon)/3.0*wvec[index_q];
    if (I1!=NULL) *I1 += q2 * m_over_T0/(1.+z) * (2*SQR(epsilon) + SQR(m_over_T0/(1.+z))) / CUB(epsilon) * wvec[index_q];
    if (I2!=NULL) *I2 += q2 * q2/CUB(epsilon) * wvec[index_q];
  }

  /** - adjust normalization */
  if (rho!=NULL) *rho *= factor2;
  if (p!=NULL) *p *= factor2;
  if (drho_dm_over_T0!=NULL) *drho_dm_over_T0 *= factor2;
  if (pseudo_p!=NULL) *pseudo_p *=factor2;
  if (I1!=NULL) *I1 *= 4*_PI_;
  if (I2!=NULL) *I2 *= 4*_PI_;

  return _SUCCESS_;
}

inline double get_mT_over_T0(struct background * pba, double phi_M){
  // Convert phi_M to J and divide by kB*T0
  double phi_M_over_T0 = phi_M * _eV_ / (_k_B_*pba->lrs_T_F*pba->T_cmb);

  return pba->lrs_m_F_over_T0 + phi_M_over_T0 * pba->lrs_g_over_M;
}

/* Function that, when equals to zero, gives the scalar field
 * (it's the derivative of the potential with respect to the field)
 *
 * phi_M    --- Scalar field times its mass [eV^2]
 * param    --- pointer to background_parameters_and_redshift struct containing background parameters and redshift
 * y        --- Output
 */
int potentialPrime(double phi_M, void *param, double *y, ErrorMsg error_message){
  struct background * pba;
  struct background_parameters_and_redshift * pbaz;
  pbaz = param;
  pba = pbaz->pba;
  
  double drhodM; // drhodM = T_0 \integ m/E
  class_call(background_lrs_momenta(pba->q_lrs_bg,
				    pba->w_lrs_bg,
				    pba->q_size_lrs_bg,
				    get_mT_over_T0(pba, phi_M),
				    pba->factor_lrs,
				    pbaz->z,
				    NULL,
				    NULL,
				    &drhodM,
				    NULL,
				    NULL,
				    NULL),
             error_message,
             error_message);
  drhodM /= _J4_to_rho_class; // Convert back to J^4
  drhodM /= pba->T_cmb*pba->lrs_T_F*_k_B_; // Divide by T0 (in J)
  drhodM /= CUB(_eV_); // Convert to eV^3
  
  *y = phi_M +  pba->lrs_g_over_M * drhodM;

  return _SUCCESS_;
}

/**
 * Returns the scalar field times its mass solving the equations of motion.
 * Hubble friction is neglected (i.e., the scalar mass is assumed to be much
 * heavier than the Hubble scale), and the resulting trascendental equation
 * is solved with the Ridder method
 *
 * @param pba    Input: pointer to background structure
 * @param z      Input: redshift
 * @param phi_M  Output: scalar field times massnumber of momenta/weights
 */
int getPhi_M(struct background * pba, double z, double * phi_M){
  double x2 = 0, x1=-pba->lrs_m_F/pba->lrs_g_over_M; // We know that -m <= g*phi <= 0. This ensures that mTilde will always be positive
  struct background_parameters_and_redshift bpaz;
  bpaz.pba = pba;
  bpaz.z = z;

  *phi_M = 0.;
  int fevals = 0;
  class_call(class_fzero_ridder(potentialPrime,
                                x1,
                                x2,
                                1e-5*MAX(fabs(x1),fabs(x2)),
                                &bpaz,
                                NULL,
                                NULL,
                                phi_M,
                                &fevals,
                                pba->error_message),
             pba->error_message,pba->error_message);

  if(pba->lrs_m_F + pba->lrs_g_over_M * (*phi_M) < 0) // mTilde < 0 because of round-off uncertainties
    *phi_M = -pba->lrs_m_F/pba->lrs_g_over_M * (1-1e-7);
      
  class_test(pba->lrs_m_F + pba->lrs_g_over_M * (*phi_M) < 0,
	     pba->error_message,
	     "Negative effective fermion mass (%.3e eV) for m_F/T_F=%.3e, g m_F/M =%.3e (phi M = %.3e eV). Check numerical precision",
	     pba->lrs_m_F + pba->lrs_g_over_M * (*phi_M),
	     pba->lrs_m_F_over_T0/(1+z),
	     pba->lrs_g_over_M * pba->lrs_m_F,
	     *phi_M);
  
  return _SUCCESS_;
}
