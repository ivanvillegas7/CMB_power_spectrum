#include"Perturbations.h"

//====================================================
// Constructors
//====================================================

Perturbations::Perturbations(
    BackgroundCosmology *cosmo, 
    RecombinationHistory *rec) : 
  cosmo(cosmo), 
  rec(rec)
{}

//====================================================
// Do all the solving
//====================================================

void Perturbations::solve(){

  // Integrate all the perturbation equation and spline the result
  integrate_perturbations();

  // Compute source functions and spline the result
  compute_source_functions();
}

//====================================================
// The main work: integrate all the perturbations
// and spline the results
//====================================================

void Perturbations::integrate_perturbations(){
  Utils::StartTiming("Integrate perturbation");

  //===================================================================
  // TODO: Set up the k-array for the k's we are going to integrate over
  // Start at k_min end at k_max with n_k points with either a
  // quadratic or a logarithmic spacing
  //===================================================================
  Vector k_array = Utils::linspace(k_min, k_max, n_k);
  Vector x_array = Utils::linspace(x_start, x_end, n_x);

  // Declare vectors storing values we are going to calculate
  Vector Phi_array(n_x*n_k);
  Vector Psi_array(n_x*n_k);
  Vector delta_CDM_array(n_x*n_k);
  Vector delta_B_array(n_x*n_k);
  Vector v_CDM_array(n_x*n_k);
  Vector v_B_array(n_x*n_k);
  std::vector<Vector> Theta_array(Constants.n_ell_theta, Vector(n_x*n_k));
  std::vector<Vector> Theta_p_array(Constants.n_ell_thetap, Vector(n_x*n_k));
  std::vector<Vector> Nu_array(Constants.n_ell_neutrinos, Vector(n_x*n_k));
  Vector Pi_array(n_x*n_k);

  double H0 = cosmo->get_H0();

  // Loop over all wavenumbers
  for(int ik = 0; ik < n_k; ik++){

    // Progress bar...
    if( (10*ik) / n_k != (10*ik+10) / n_k ) {
      std::cout << "Progress: " << (100*ik+100)/n_k << "%\n" << std::flush;
      if(ik == n_k-1) std::cout << std::endl;
    }

    // Current value of k
    double k = k_array[ik];

    // Find value to integrate to
    double x_end_tight = get_tight_coupling_time(k);

    //===================================================================
    // TODO: Tight coupling integration
    // Remember to implement the routines:
    // set_ic : The IC at the start
    // rhs_tight_coupling_ode : The dydx for our coupled ODE system
    //===================================================================

    // Set up initial conditions in the tight coupling regime
    auto y_tight_coupling_ini = set_ic(x_start, k);

    // The tight coupling ODE system
    ODEFunction dydx_tight_coupling = [&](double x, const double *y, double *dydx){
      return rhs_tight_coupling_ode(x, k, y, dydx);
    };

    // Integrate from x_start -> x_end_tight
    int index_end_tc = int((x_end_tight-x_start)/((x_end-x_start)/n_x));
    Vector x_array_tc = Utils::linspace(x_start, x_end_tight, index_end_tc);

    // Solve the ODE in the tight coupling regime
    ODESolver ode;
    ode.solve(dydx_tight_coupling, x_array_tc, y_tight_coupling_ini);
    auto solution_tc     = ode.get_data();
    auto y_tc_end = ode.get_final_data();

    //====================================================================
    // TODO: Full equation integration
    // Remember to implement the routines:
    // set_ic_after_tight_coupling : The IC after tight coupling ends
    // rhs_full_ode : The dydx for our coupled ODE system
    //===================================================================

    // Set up initial conditions (y_tight_coupling is the solution at the end of tight coupling)
    auto y_full_ini = set_ic_after_tight_coupling(y_tc_end, x_end_tight, k);

    // The full ODE system
    ODEFunction dydx_full = [&](double x, const double *y, double *dydx){
      return rhs_full_ode(x, k, y, dydx);
    };

    // Integrate from x_end_tight -> x_end
    Vector x_array_full = Utils::linspace(x_end_tight, x_end, n_x-index_end_tc);

    // Solve the ODE in the tight coupling regime
    ode.solve(dydx_full, x_array_full, y_full_ini);
    auto solution_full = ode.get_data();

    //===================================================================
    // TODO: remember to store the data found from integrating so we can
    // spline it below.
    //
    // To compute a 2D spline of a function f(x,k) the data must be given 
    // to the spline routine as a 1D array f_array with the points f(ix, ik) 
    // stored as f_array[ix + n_x * ik]
    //
    // Example:
    // Vector x_array(n_x);
    // Vector k_array(n_k);
    // Vector f(n_x * n_k);
    // Spline2D y_spline;
    // f_spline.create(x_array, k_array, f_array);
    // We can now use the spline as f_spline(x, k)
    //
    // NB: If you use Theta_spline then you have to allocate it first,
    // before using it e.g.
    // Theta_spline = std::vector<Spline2D>(n_ell_theta);
    //===================================================================

    // Start filling arrays from the tight coupling regime 
    for (int ix=0; ix < index_end_tc; ix++){
      int index = ix+n_x*ik;
      auto y_tc = solution_tc[ix];
      double x  = x_array_tc[ix];

      // Get quantities needed for computation
      double Hp         = cosmo->Hp_of_x(x);
      double Omega_R    = cosmo->get_OmegaR();
      double Omega_Nu   = cosmo->get_OmegaNu();
      double dtaudx     = rec->dtaudx_of_x(x);
      double ck_over_Hp = Constants.c*k/Hp ;

      // Copy from below (line "References to the tight coupling quantities")
      double &delta_CDM =  y_tc[Constants.ind_deltacdm_tc];
      double &delta_B   =  y_tc[Constants.ind_deltab_tc];
      double &v_CDM     =  y_tc[Constants.ind_vcdm_tc];
      double &v_B       =  y_tc[Constants.ind_vb_tc];
      double &Phi       =  y_tc[Constants.ind_Phi_tc];
      double *Theta     = &y_tc[Constants.ind_start_theta_tc];

      // Calculate quantities
      Phi_array[index]       = Phi;
      delta_CDM_array[index] = delta_CDM;
      delta_B_array[index]   = delta_B;
      v_CDM_array[index]     = v_CDM;
      v_B_array[index]       = v_B;
      Theta_array[0][index]  = Theta[0];
      Theta_array[1][index]  = Theta[1];
      if (Constants.polarization) Theta_array[2][index] = -8.*ck_over_Hp/(15.*dtaudx)*Theta[1];
      else Theta_array[2][index] = -20.*ck_over_Hp/(45.*dtaudx)*Theta[1];
      // Solve for rest of the l-values
      for(int l=3; l < Constants.n_ell_theta; l++){
        Theta_array[l][index] = -l/(2.*l+1.)*ck_over_Hp/dtaudx*Theta_array[l-1][index];
      }
      if(Constants.polarization){
        Theta_p_array[0][index]   = 5.*Theta[2]/4.;
        Theta_p_array[1][index]   = -ck_over_Hp*Theta[2]/(4.*dtaudx);
        Theta_p_array[2][index]   = Theta[2]/4.;
        for (int l = 3; l < Constants.n_ell_thetap; l++)
        {
          Theta_p_array[l][index] = -l*ck_over_Hp*Theta_p_array[l-1][index]/((2.*l+1.)*dtaudx);
        }
      }
      double N2 = 0.0;
      if(Constants.neutrinos){
        double *Nu           = &y_tc[Constants.ind_start_nu_tc];
        Nu_array[0][index]   = -Psi_array[index]/2.;
        Nu_array[1][index]   = ck_over_Hp*Psi_array[index]/6.;
        Nu_array[2][index]   = -(Psi_array[index]+Phi)*pow(Constants.c*k*exp(x)/H0, 2)/(12.*Omega_Nu);
        N2                   = Nu_array[2][index];
        for (int l = 3; l < Constants.n_ell_neutrinos_tc; l++){
          Nu_array[l][index] = ck_over_Hp*Nu_array[l-1][index]/(2.*l+1.);
        }
      }
      Psi_array[index]       = -Phi-12.*pow(H0/(Constants.c*k*exp(x)), 2)*(Omega_R*Theta[2]+Omega_Nu*N2);
    }
    // Now fill rest of arrays from the full regime 
    for (int ix = index_end_tc; ix < n_x; ix++){
      int index   = ix+n_x*ik;
      auto y_full = solution_full[ix-index_end_tc];
      double x    = x_array_full[ix-index_end_tc];

      // Get quantities needed for computation
      double Hp         = cosmo->Hp_of_x(x);
      double OmegaR     = cosmo->get_OmegaR();
      double OmegaNu    = cosmo->get_OmegaNu();
      double dtau       = rec->dtaudx_of_x(x);

      double &delta_CDM =  y_full[Constants.ind_deltacdm];
      double &delta_B   =  y_full[Constants.ind_deltab];
      double &v_CDM     =  y_full[Constants.ind_vcdm];
      double &v_B       =  y_full[Constants.ind_vb];
      double &Phi       =  y_full[Constants.ind_Phi];
      double *Theta     = &y_full[Constants.ind_start_theta];
      double N2         = 0.0;
      if (Constants.neutrinos)
      {
        double *Nu      = &y_full[Constants.ind_start_nu];
        N2              = Nu[2];
      }
      

      // Calculate quantities
      Phi_array[index]            = Phi;
      Psi_array[index]            = -Phi-12.*pow(H0/(Constants.c*k*exp(x)), 2)*OmegaR*Theta[2]+OmegaNu*N2;
      delta_CDM_array[index]      = delta_CDM;
      delta_B_array[index]        = delta_B;
      v_CDM_array[index]          = v_CDM;
      v_B_array[index]            = v_B;
      // Fill in all solved Theta-values
      for(int l=0; l < Constants.n_ell_theta; l++){
        Theta_array[l][index]     = Theta[l];
      }

      Pi_array[index]             = Theta_array[2][index];

      // Fill in all solved Theta_p-values
      if (Constants.polarization)
      {
        double *Theta_p           = &y_full[Constants.ind_start_thetap];
        for(int l=0; l < Constants.n_ell_thetap; l++)
        {
          Theta_p_array[l][index] = Theta_p[l];
        }
        Pi_array[index]            += Theta_p_array[0][index]+Theta_p_array[2][index];
      }
      // Fill in all solved Nu-values
      if (Constants.neutrinos){
        double *Nu                = &y_full[Constants.ind_start_nu];
        for(int l=0; l < Constants.n_ell_neutrinos; l++)
        {
          Nu_array[l][index]      = Nu[l];
        }
      }
    }
  }
  Utils::EndTiming("Integrate perturbation");

  //=============================================================================
  // TODO: Make all splines needed: Θ_0,Θ_1,Θ_2,Φ,Ψ,...
  //=============================================================================
  Phi_spline.create(x_array, k_array, Phi_array, "Φ spline");
  Psi_spline.create(x_array, k_array, Psi_array, "Ψ spline");
  delta_CDM_spline.create(x_array, k_array, delta_CDM_array, "δ_CDM spline");
  delta_B_spline.create(x_array, k_array, delta_B_array, "δ_B spline");
  v_CDM_spline.create(x_array, k_array, v_CDM_array, "v_CDM spline");
  v_B_spline.create(x_array, k_array, v_B_array, "v_B spline");
  Pi_spline.create(x_array, k_array, Pi_array, "Π spline");
  Theta_spline = std::vector<Spline2D>(Constants.n_ell_theta);
  for(int l=0; l < Constants.n_ell_theta; l++){
    Theta_spline[l].create(x_array, k_array, Theta_array[l]);
  }
  if (Constants.polarization)
  {
    Theta_p_spline = std::vector<Spline2D>(Constants.n_ell_thetap);
    for(int l=0; l < Constants.n_ell_theta; l++){
      Theta_p_spline[l].create(x_array, k_array, Theta_p_array[l]);
    }
  }
  if (Constants.neutrinos)
  {
    Nu_spline = std::vector<Spline2D>(Constants.n_ell_neutrinos);
    for(int l=0; l < Constants.n_ell_neutrinos; l++){
      Nu_spline[l].create(x_array, k_array, Nu_array[l]);
    }
  }
}

//====================================================
// Set IC at the start of the run (this is in the
// tight coupling regime)
//====================================================
Vector Perturbations::set_ic(const double x, const double k) const{

  // The vector we are going to fill
  Vector y_tc(Constants.n_ell_tot_tc);

  //=============================================================================
  // Compute where in the y_tc array each component belongs
  // This is just an example of how to do it to make it easier
  // Feel free to organize the component any way you like
  //=============================================================================
  
  // For integration of perturbations in tight coupling regime (Only 2 photon multipoles + neutrinos needed)
  const int n_ell_theta_tc      = Constants.n_ell_theta_tc;
  const int n_ell_neutrinos_tc  = Constants.n_ell_neutrinos_tc;
  const int n_ell_tot_tc        = Constants.n_ell_tot_tc;
  const bool polarization       = Constants.polarization;
  const bool neutrinos          = Constants.neutrinos;

  // References to the tight coupling quantities
  double &delta_CDM    =  y_tc[Constants.ind_deltacdm_tc];
  double &delta_B      =  y_tc[Constants.ind_deltab_tc];
  double &v_CDM        =  y_tc[Constants.ind_vcdm_tc];
  double &v_B          =  y_tc[Constants.ind_vb_tc];
  double &Phi          =  y_tc[Constants.ind_Phi_tc];
  double *Theta        = &y_tc[Constants.ind_start_theta_tc];

  //=============================================================================
  // TODO: Set the initial conditions in the tight coupling regime
  //=============================================================================

  // Set constants
  double c                   = Constants.c;

  // Cosmological parameters and variables
  double H0                  = cosmo->get_H0();
  double Hp                  = cosmo->Hp_of_x(x);
  double Omega_R             = cosmo->get_OmegaR();
  double Omega_Nu            = cosmo->get_OmegaNu();

  //Set recombination parameters
  double R                   = rec->R_of_x(x);
  double dtaudx              = rec->dtaudx_of_x(x);

  // Save flops
  double ck_over_Hp          = c*k/Hp;
  double f_Nu                = 0.0;
  if (neutrinos){
    f_Nu                     = Omega_Nu/(Omega_R+Omega_Nu);
  }

  // SET: Scalar quantities (Gravitational potental, baryons and CDM)
  double Psi                 = -1./(3./2.+2.*f_Nu/5.);
  delta_CDM                  = -3.*Psi/2.;
  delta_B                    = -3.*Psi/2.;
  v_CDM                      = -ck_over_Hp*Psi/2.;
  v_B                        = -ck_over_Hp*Psi/2.;
  Phi                        = -(1.+2.*f_Nu/5.)*Psi;

  // SET: Photon temperature perturbations (Theta_ell)
  Theta[0]                   = -Psi/2.;
  Theta[1]                   = ck_over_Hp*Psi/6.;
  
  // SET: Neutrino perturbations (N_ell)
  if(neutrinos){
    double *Nu               = &y_tc[Constants.ind_start_nu_tc];
    Nu[0]                    = -Psi/2.;
    Nu[1]                    = ck_over_Hp*Psi/6.;
    Nu[2]                    = -(Psi+Phi)*pow(c*k*exp(x)/H0, 2)/(12.*Omega_Nu);
    for (int l = 3; l < n_ell_neutrinos_tc; l++)
    {
      Nu[l]                  = ck_over_Hp*Nu[l-1]/(2.*l+1.);
    }
  }

  return y_tc;
}

//====================================================
// Set IC for the full ODE system after tight coupling 
// regime ends
//====================================================

Vector Perturbations::set_ic_after_tight_coupling(
    const Vector &y_tc, 
    const double x, 
    const double k) const{

  // Make the vector we are going to fill
  Vector y(Constants.n_ell_tot_full);
  
  //=============================================================================
  // Compute where in the y array each component belongs and where corresponding
  // components are located in the y_tc array
  // This is just an example of how to do it to make it easier
  // Feel free to organize the component any way you like
  //=============================================================================

  // Number of multipoles we have in the full regime
  const int n_ell_theta         = Constants.n_ell_theta;
  const int n_ell_thetap        = Constants.n_ell_thetap;
  const int n_ell_neutrinos     = Constants.n_ell_neutrinos;
  const bool polarization       = Constants.polarization;
  const bool neutrinos          = Constants.neutrinos;

  // Number of multipoles we have in the tight coupling regime
  const int n_ell_theta_tc      = Constants.n_ell_theta_tc;
  const int n_ell_neutrinos_tc  = Constants.n_ell_neutrinos_tc;

  // References to the tight coupling quantities
  const double &delta_CDM_tc    =  y_tc[Constants.ind_deltacdm_tc];
  const double &delta_B_tc      =  y_tc[Constants.ind_deltab_tc];
  const double &v_CDM_tc        =  y_tc[Constants.ind_vcdm_tc];
  const double &v_B_tc          =  y_tc[Constants.ind_vb_tc];
  const double &Phi_tc          =  y_tc[Constants.ind_Phi_tc];
  const double *Theta_tc        = &y_tc[Constants.ind_start_theta_tc];

  // References to the quantities we are going to set
  double &delta_CDM       =  y[Constants.ind_deltacdm_tc];
  double &delta_B         =  y[Constants.ind_deltab_tc];
  double &v_CDM           =  y[Constants.ind_vcdm_tc];
  double &v_B             =  y[Constants.ind_vb_tc];
  double &Phi             =  y[Constants.ind_Phi_tc];
  double *Theta           = &y[Constants.ind_start_theta_tc];

  //=============================================================================
  // TODO: fill in the initial conditions for the full equation system below
  // NB: remember that we have different number of multipoles in the two
  // regimes so be careful when assigning from the tc array
  //=============================================================================

  // Set constants
  double c                    = Constants.c;

  // Cosmological parameters and variables
  double H0                   = cosmo->get_H0();
  double Hp                   = cosmo->Hp_of_x(x);
  double Omega_R              = cosmo->get_OmegaR();
  double Omega_Nu             = cosmo->get_OmegaNu();

  //Set recombination parameters
  double R                    = rec->R_of_x(x);
  double dtaudx               = rec->dtaudx_of_x(x);

  // Save flops
  double ck_over_Hp           = c*k/Hp;
  double f_Nu                 = 0.0;
  if (neutrinos){
    f_Nu                      = Omega_Nu/(Omega_R+Omega_Nu);
  }

  // SET: Scalar quantities (Gravitational potental, baryons and CDM)
  delta_CDM                   = delta_CDM_tc;
  delta_B                     = delta_B_tc;
  v_CDM                       = v_CDM_tc;
  v_B                         = v_B_tc;
  Phi                         = Phi_tc;

  // SET: Photon temperature perturbations (Theta_ell)
  Theta[0]                   = Theta_tc[0];
  Theta[1]                   = Theta_tc[1];
  if (polarization) Theta[2] = -(8./15.)*ck_over_Hp*Theta[1]/dtaudx;
  else Theta[2]              = -(20./45.)*ck_over_Hp*Theta[1]/dtaudx;
  for (int l = 3; l < n_ell_theta; l++)
  {
    Theta[l]                 = -l*ck_over_Hp*Theta[l-1]/((2.*l+1.)*dtaudx);
  }

  // SET: Photon polarization perturbations (Theta_p_ell)
  if(polarization){
    double *Theta_p          = &y[Constants.ind_start_thetap_tc];
    Theta_p[0]               = 5.*Theta[2]/4.;
    Theta_p[1]               = -ck_over_Hp*Theta[2]/(4.*dtaudx);
    Theta_p[2]               = Theta[2]/4.;
    for (int l = 3; l < n_ell_thetap; l++)
    {
      Theta_p[l]             = -l*ck_over_Hp*Theta_p[l-1]/((2.*l+1.)*dtaudx);
    }
  }

  // SET: Neutrino perturbations (N_ell)
  if(neutrinos){
    const double *Nu_tc      = &y_tc[Constants.ind_start_nu_tc];
    double *Nu               = &y[Constants.ind_start_nu_tc];
    for (int l = 0; l < n_ell_neutrinos; l++)
    {
      Nu[l]                  = Nu_tc[l];
    }
  }
  
  return y;
}

//====================================================
// The time when tight coupling end
//====================================================

double Perturbations::get_tight_coupling_time(const double k) const{
  double x_tight_coupling_end = 0.0;

  //=============================================================================
  // TODO: compute and return x for when tight coupling ends
  // Remember all the three conditions in Callin
  //=============================================================================
  double dtaudx;
  Vector x_array = Utils::linspace(x_start, x_end, n_x);

  for (int i = 0; i < n_x; i++)
  {
    dtaudx = -rec->dtaudx_of_x(x_array[i]);

    if (dtaudx<10.0 || dtaudx<10.0*k*Constants.c/cosmo->Hp_of_x(x_array[i]) || rec->Xe_of_x(x_array[i])<0.99)
    {
      return x_array[i-1];
    }
  }

  return x_tight_coupling_end;
}

//====================================================
// After integrsating the perturbation compute the
// source function(s)
//====================================================
void Perturbations::compute_source_functions(){
  Utils::StartTiming("Source");

  //=============================================================================
  // TODO: Make the x and k arrays to evaluate over and use to make the splines
  //=============================================================================
  Vector k_array = Utils::linspace(k_min, k_max, n_k);
  Vector x_array = Utils::linspace(x_start, x_end, n_x);

  // Make storage for the source functions (in 1D array to be able to pass it to the spline)
  Vector ST_array(k_array.size()*x_array.size());
  Vector SE_array(k_array.size()*x_array.size());

  // Compute source functions
  for(auto ix = 0; ix < x_array.size(); ix++){
    const double x = x_array[ix];
    for(auto ik = 0; ik < k_array.size(); ik++){
      const double k = k_array[ik];

      // NB: This is the format the data needs to be stored 
      // in a 1D array for the 2D spline routine source(ix,ik) -> S_array[ix + nx * ik]
      const int index = ix+n_x*ik;

      //=============================================================================
      // TODO: Compute the source functions
      //=============================================================================

      // Fetch functions from BackgroundCosmology
      double Hp      = cosmo->Hp_of_x(x);
      double dHpdx   = cosmo->dHpdx_of_x(x);
      double ddHpddx = cosmo->ddHpddx_of_x(x);

      // Fetch functions from Recombination
      double tau      = rec->tau_of_x(x);
      double dtaudx   = rec->dtaudx_of_x(x);
      double ddtauddx = rec->ddtauddx_of_x(x);
      double g        = rec->g_tilde_of_x(x);
      double dgdx     = rec->dgdx_tilde_of_x(x);
      double ddgddx   = rec->ddgddx_tilde_of_x(x);

      // Get functions from Perturbation
      double Psi   = get_Psi(x,k);
      double v_B   = get_v_B(x,k);
      double T0    = get_Theta(x,k,0);
      double T1    = get_Theta(x,k,1);
      double T2    = get_Theta(x,k,2);
      double T3    = get_Theta(x,k,3);
      double P1    = 0.0;
      double P2    = 0.0;
      double P3    = 0.0;
      double Pi    = get_Pi(x, k);
      double dPidx = get_dPidx(x, k);

      double dPsidx = get_dPsidx(x,k);
      double dPhidx = get_dPhidx(x,k);
      double dv_Bdx = get_dv_Bdx(x,k);
      double dT1dx  = get_dThetadx(x,k,1);
      double dT2dx  = get_dThetadx(x,k,2);
      double dT3dx  = get_dThetadx(x,k,3);
      double dP1dx  = 0.0;
      double dP2dx  = 0.0;
      double dP3dx  = 0.0;

      if (Constants.polarization)
      {
        P1    = get_Theta_p(x,k,1);
        P2    = get_Theta_p(x,k,2);
        P3    = get_Theta_p(x,k,3);
        dP1dx  = get_dTheta_pdx(x,k,1);
        dP2dx  = get_dTheta_pdx(x,k,2);
        dP3dx  = get_dTheta_pdx(x,k,3);
      }
      
      // Save FLOPs
      double c           = Constants.c;
      double ck_over_Hp  = c*k/Hp;
      double dHp_over_Hp = dHpdx/Hp;

      //double ddT2ddx = ck_over_Hp/5.*(2.*dT1dx-3.*dT3dx)-ck_over_Hp*dHpdx/(5.*Hp)*(2*T1-3*T3)+9./10.*(ddtauddx*T2+dtaudx*dT2dx);

      double ddPiddx_1 = ((2.*k*c)/(5.*Hp))*(-dHp_over_Hp*T1+dT1dx);
      double ddPiddx_2 = (3./10.)*(ddtauddx*Pi+dtaudx*dPidx);
      double ddPiddx_3 = ((3.*k*c)/(5.*Hp))*(-dHp_over_Hp*(T3+P1+P3)+dT3dx+dP1dx+dP3dx);
      double ddPiddx   = ddPiddx_1+ddPiddx_2-ddPiddx_3;

      
      // Temperature source 
      double term1     = g*(T0+Psi+Pi/4.);
      double term2     = exp(-tau)*(dPsidx-dPhidx);
      double term3     = (dHpdx*g*v_B+Hp*dgdx*v_B+Hp*g*dv_Bdx)/(k*c);
      double term4     = (pow(dHpdx, 2)+Hp*ddHpddx)*g*Pi+3*Hp*dHpdx*(dgdx*Pi+g*dPidx)+pow(Hp, 2)*(ddgddx*Pi+2*dgdx*dPidx+g*ddPiddx);

      // For SW term
      //ST_array[index]  = term1;
      // For ISW term
      //ST_array[index]  = term2;
      // For DOPPLER term
      //ST_array[index]  = -term3;
      // For POLARIZATION term (if polarization is included)
      //ST_array[index]  = term4;
      // All terms
      ST_array[index]  = term1+term2-term3+term4;

      // Polarization source
      if(Constants.polarization){
        SE_array[index] = 3.*g*Pi/(4.*pow(k*(cosmo->eta_of_x(0.0)-cosmo->eta_of_x(x)), 2));
      }
      else SE_array[index] = 0.0;
    }
  }

  // Spline the source functions
  ST_spline.create (x_array, k_array, ST_array, "Source_Temp_x_k");
  SE_spline.create (x_array, k_array, SE_array, "Source_Pol_x_k");

  Utils::EndTiming("Source");
}

//====================================================
// The right hand side of the perturbations ODE
// in the tight coupling regime
//====================================================

// Derivatives in the tight coupling regime
int Perturbations::rhs_tight_coupling_ode(double x, double k, const double *y, double *dydx){

  //=============================================================================
  // Compute where in the y / dydx array each component belongs
  // This is just an example of how to do it to make it easier
  // Feel free to organize the component any way you like
  //=============================================================================

  // Set constants
  double c                      = Constants.c;

  // Set cosmological parameters and variables
  double H0                     = cosmo->get_H0();
  double Hp                     = cosmo->Hp_of_x(x);
  double dHpdx                  = cosmo->dHpdx_of_x(x);
  double Omega_CDM              = cosmo->get_OmegaCDM();
  double Omega_B                = cosmo->get_OmegaB();
  double Omega_R                = cosmo->get_OmegaR();
  double Omega_Nu               = cosmo->get_OmegaNu();
  double eta                    = cosmo->eta_of_x(x);

  //Set recombination parameters
  double R                      = rec->R_of_x(x);
  double dtaudx                 = rec->dtaudx_of_x(x);
  double ddtauddx               = rec->ddtauddx_of_x(x);

  // For integration of perturbations in tight coupling regime (Only 2 photon multipoles + neutrinos needed)
  const int n_ell_theta_tc      = Constants.n_ell_theta_tc;
  const int n_ell_neutrinos_tc  = Constants.n_ell_neutrinos_tc;
  const bool neutrinos          = Constants.neutrinos;

  // The different quantities in the y array
  const double &delta_CDM       =  y[Constants.ind_deltacdm_tc];
  const double &delta_B         =  y[Constants.ind_deltab_tc];
  const double &v_CDM           =  y[Constants.ind_vcdm_tc];
  const double &v_B             =  y[Constants.ind_vb_tc];
  const double &Phi             =  y[Constants.ind_Phi_tc];
  const double *Theta           = &y[Constants.ind_start_theta_tc];

  // References to the quantities we are going to set in the dydx array
  double &ddelta_CDMdx          =  dydx[Constants.ind_deltacdm_tc];
  double &ddelta_Bdx            =  dydx[Constants.ind_deltab_tc];
  double &dv_CDMdx              =  dydx[Constants.ind_vcdm_tc];
  double &dv_Bdx                =  dydx[Constants.ind_vb_tc];
  double &dPhidx                =  dydx[Constants.ind_Phi_tc];
  double *dThetadx              = &dydx[Constants.ind_start_theta_tc];

  // Save flops
  double ck_over_Hp             = c*k/Hp;
  double dHp_over_Hp            = dHpdx/Hp;

  double Nu0                    = 0.0;
  double Nu2                    = 0.0;
  if (neutrinos){
    const double *Nu            = &y[Constants.ind_start_nu_tc];
    Nu0                         = Nu[0];
    Nu2                         = Nu[2];
  }
  double Theta2                 = -ck_over_Hp*Theta[1]/dtaudx;
  if (Constants.polarization)
  {
    Theta2                     *= 8./15.;
  }
  else Theta2                  *= 20./45.;
  
  //=============================================================================
  // TODO: fill in the expressions for all the derivatives
  //=============================================================================

  // SET: Scalar quantities (Φ, δ, v, ...)
  double Psi                  = -Phi-12.*pow(H0/(c*k*exp(x)), 2)*(Omega_R*Theta2+Omega_Nu*Nu2);
  
  dPhidx                      = Psi-pow(ck_over_Hp, 2)*Phi/3.+pow(H0/Hp, 2)*exp(-x)*(Omega_CDM*delta_CDM+Omega_B*delta_B+4.*Omega_R*exp(-x)*Theta[0]+4.*Omega_Nu*exp(-x)*Nu0)/2.;
  ddelta_Bdx                  = ck_over_Hp*v_B-3.*dPhidx;
  ddelta_CDMdx                = ck_over_Hp*v_CDM-3.*dPhidx;
  dv_CDMdx                    = -v_CDM-ck_over_Hp*Psi;

  // SET: Photon multipoles (Θ_ell)
  dThetadx[0]                 = -ck_over_Hp*Theta[1]-dPhidx;

  double num                  = ((1.-R)*dtaudx+(1+R)*ddtauddx)*(3.*Theta[1]+v_B)-ck_over_Hp*Psi+(1.-dHp_over_Hp)*ck_over_Hp*(-Theta[0]+2.*Theta2)-ck_over_Hp*dThetadx[0];
  double den                  = (1.+R)*dtaudx+dHp_over_Hp-1.;
  double q                    = -num/den;

  dv_Bdx                      = (-v_B-ck_over_Hp*Psi+R*(q+ck_over_Hp*(-Theta[0]+2.*Theta2)-ck_over_Hp*Psi))/(1.+R);
  dThetadx[1]                 = (q-dv_Bdx)/3.;

  // SET: Neutrino mutlipoles (Nu_ell)
  if(neutrinos){
    const double *Nu            = &y[Constants.ind_start_nu_tc];
    double *dNudx               = &dydx[Constants.ind_start_nu_tc];
    dNudx[0]                    = -ck_over_Hp*Nu[1]-dPhidx;
    dNudx[1]                    = ck_over_Hp*Nu[0]/3.-(2./3.)*ck_over_Hp*Nu[2]+ck_over_Hp*Psi/3.;
    for (int l = 2; l < n_ell_neutrinos_tc-1; l++){
      dNudx[l]                  = ck_over_Hp/(2.*l+1.)*(l*Nu[l-1]-(l+1.)*Nu[l+1]);
    }
    dNudx[n_ell_neutrinos_tc-1] = ck_over_Hp*Nu[n_ell_neutrinos_tc-2]-c*n_ell_neutrinos_tc*Nu[n_ell_neutrinos_tc-1]/(Hp*eta);
  }
  
  return GSL_SUCCESS;
}

//====================================================
// The right hand side of the full ODE
//====================================================

int Perturbations::rhs_full_ode(double x, double k, const double *y, double *dydx){

  //=============================================================================
  // Compute where in the y / dydx array each component belongs
  // This is just an example of how to do it to make it easier
  // Feel free to organize the component any way you like
  //=============================================================================

  // Index and number of the different quantities
  const int n_ell_theta     = Constants.n_ell_theta;
  const int n_ell_thetap    = Constants.n_ell_thetap;
  const int n_ell_neutrinos = Constants.n_ell_neutrinos;
  const bool polarization   = Constants.polarization;
  const bool neutrinos      = Constants.neutrinos;

  // The different quantities in the y array
  const double &delta_CDM   =  y[Constants.ind_deltacdm];
  const double &delta_B     =  y[Constants.ind_deltab];
  const double &v_CDM       =  y[Constants.ind_vcdm];
  const double &v_B         =  y[Constants.ind_vb];
  const double &Phi         =  y[Constants.ind_Phi];
  const double *Theta       = &y[Constants.ind_start_theta];

  // References to the quantities we are going to set in the dydx array
  double &ddelta_CDMdx      =  dydx[Constants.ind_deltacdm];
  double &ddelta_Bdx        =  dydx[Constants.ind_deltab];
  double &dv_CDMdx          =  dydx[Constants.ind_vcdm];
  double &dv_Bdx            =  dydx[Constants.ind_vb];
  double &dPhidx            =  dydx[Constants.ind_Phi];
  double *dThetadx          = &dydx[Constants.ind_start_theta];

  // Check if polarization and neutrinos are included
  double ThetaP0            = 0.0;
  double ThetaP2            = 0.0;
  if (polarization)
  {
    const double *Theta_p   = &y[Constants.ind_start_thetap];
    ThetaP0                 = Theta_p[0];
    ThetaP2                 = Theta_p[2];
  }
  double Nu0                = 0.0;
  double Nu2                = 0.0;
  if (neutrinos)
  {
    const double *Nu        = &y[Constants.ind_start_thetap];
    Nu0                     = Nu[0];
    Nu2                     = Nu[2];
  }

  // Set constants
  double c                  = Constants.c;

  // Cosmological parameters and variables
  double H0                 = cosmo->get_H0();
  double Hp                 = cosmo->Hp_of_x(x);
  double Omega_CDM          = cosmo->get_OmegaCDM();
  double Omega_B            = cosmo->get_OmegaB();
  double Omega_R            = cosmo->get_OmegaR();
  double Omega_Nu           = cosmo->get_OmegaNu();
  double eta                = cosmo->eta_of_x(x);

  //Set recombination parameters
  double R                  = rec->R_of_x(x);
  double dtaudx             = rec->dtaudx_of_x(x);

  // Save flops
  double ck_over_Hp         = c*k/Hp;
  double Pi                 = Theta[2]+ThetaP0+ThetaP2;

  //=============================================================================
  // TODO: fill in the expressions for all the derivatives
  //=============================================================================

  // SET: Scalar quantities (Φ, δ, v, ...)
  double Psi               = -Phi-12.*pow(H0/(c*k*exp(x)), 2)*(Omega_R*Theta[2]+Omega_Nu*Nu2);
  dPhidx                   = Psi-pow(ck_over_Hp, 2)*Phi/3.+pow(H0/Hp, 2)*exp(-x)*(Omega_CDM*delta_CDM+Omega_B*delta_B+4*Omega_R*exp(-x)*Theta[0]+4.*Omega_Nu*exp(-x)*Nu0)/2.;
  dv_Bdx                   = -v_B-ck_over_Hp*Psi+dtaudx*R*(3.*Theta[1]+v_B);
  ddelta_Bdx               = ck_over_Hp*v_B-3.*dPhidx;
  dv_CDMdx                 = -v_CDM-ck_over_Hp*Psi;
  ddelta_CDMdx             = ck_over_Hp*v_CDM-3.*dPhidx;

  // SET: Photon multipoles (Θ_ell)
  dThetadx[0]              = -ck_over_Hp*Theta[1]-dPhidx;
  dThetadx[1]              = ck_over_Hp*Theta[0]/3.-(2./3.)*ck_over_Hp*Theta[2]+ck_over_Hp*Psi/3.+dtaudx*(Theta[1]+v_B/3.);
  for (int l = 2; l < n_ell_theta-1; l++)
  {
    if (l==2) dThetadx[l]  = ck_over_Hp/(2.*l+1.)*(l*Theta[l-1]-(l+1.)*Theta[l+1])+dtaudx*(Theta[l]-Pi/10.);
    else dThetadx[l]       = ck_over_Hp/(2.*l+1.)*(l*Theta[l-1]-(l+1.)*Theta[l+1])+dtaudx*Theta[l];   
  }
  dThetadx[n_ell_theta-1]  = ck_over_Hp*Theta[n_ell_theta-2]+(dtaudx-c*(n_ell_theta/(Hp*eta)))*Theta[n_ell_theta-1];

  // SET: Photon polarization multipoles (Θ_p_ell)
  if(polarization){
    const double *Theta_p  = &y[Constants.ind_start_thetap];
    double *dTheta_pdx     = &dydx[Constants.ind_start_thetap];
    dTheta_pdx[0]          = -ck_over_Hp*Theta_p[1]+dtaudx*(Theta_p[0]-Pi/2.);
    for (int l = 2; l < n_ell_thetap-1; l++)
  {
    if (l==2){
      dTheta_pdx[l]          = ck_over_Hp/(2.*l+1.)*(l*Theta_p[l-1]-(l+1.)*Theta_p[l+1])+dtaudx*(Theta_p[l]-Pi/10.);
    }
    else{
      dTheta_pdx[l]          = ck_over_Hp/(2.*l+1.)*(l*Theta_p[l-1]-(l+1.)*Theta_p[l+1])+dtaudx*Theta_p[l];
    }
  }
  dTheta_pdx[n_ell_thetap-1] = ck_over_Hp*Theta_p[n_ell_thetap-2]+(dtaudx-c*n_ell_thetap/(Hp*eta))*Theta_p[n_ell_thetap-1];
  }

  // SET: Neutrino mutlipoles (Nu_ell)
  if(neutrinos){
    const double *Nu         = &y[Constants.ind_start_nu];
    double *dNudx            = &dydx[Constants.ind_start_nu];
    dNudx[0]                 = -ck_over_Hp*Nu[1]-dPhidx;
    dNudx[1]                 = ck_over_Hp*Nu[0]/3.-(2./3.)*ck_over_Hp*Nu[2]+ck_over_Hp*Psi/3.;
    for (int l = 2; l < n_ell_neutrinos-1; l++){
      dNudx[l]               = ck_over_Hp/(2.*l+1.)*(l*Nu[l-1]-(l+1.)*Nu[l+1]);
    }
    dNudx[n_ell_neutrinos-1] = ck_over_Hp*Nu[n_ell_neutrinos-2]-c*n_ell_neutrinos*Nu[n_ell_neutrinos-1]/(Hp*eta);
  }
  
  return GSL_SUCCESS;
}

//====================================================
// Get methods
//====================================================

double Perturbations::get_delta_CDM(const double x, const double k) const{
  return delta_CDM_spline(x,k);
}
double Perturbations::get_delta_B(const double x, const double k) const{
  return delta_B_spline(x,k);
}
double Perturbations::get_v_CDM(const double x, const double k) const{
  return v_CDM_spline(x,k);
}
double Perturbations::get_v_B(const double x, const double k) const{
  return v_B_spline(x,k);
}
double Perturbations::get_dv_Bdx(const double x, const double k) const{
  return v_B_spline.deriv_x(x,k);
}
double Perturbations::get_Phi(const double x, const double k) const{
  return Phi_spline(x,k);
}
double Perturbations::get_dPhidx(const double x, const double k) const{
  return Phi_spline.deriv_x(x,k);
}
double Perturbations::get_Psi(const double x, const double k) const{
  return Psi_spline(x,k);
}
double Perturbations::get_dPsidx(const double x, const double k) const{
  return Psi_spline.deriv_x(x,k);
}
double Perturbations::get_Pi(const double x, const double k) const{
  return Pi_spline(x,k);
}
double Perturbations::get_dPidx(const double x, const double k) const{
  return Pi_spline.deriv_x(x,k);
}
double Perturbations::get_Source_T(const double x, const double k) const{
  return ST_spline(x,k);
}
double Perturbations::get_Source_E(const double x, const double k) const{
  return SE_spline(x,k);
}
double Perturbations::get_Theta(const double x, const double k, const int ell) const{
  return Theta_spline[ell](x,k);
}
double Perturbations::get_dThetadx(const double x, const double k, const int ell) const{
  return Theta_spline[ell].deriv_x(x,k);
}
double Perturbations::get_Theta_p(const double x, const double k, const int ell) const{
  return Theta_p_spline[ell](x,k);
}
double Perturbations::get_dTheta_pdx(const double x, const double k, const int ell) const{
  return Theta_p_spline[ell].deriv_x(x,k);
}
double Perturbations::get_Nu(const double x, const double k, const int ell) const{
  return Nu_spline[ell](x,k);
}

//====================================================
// Print some useful info about the class
//====================================================

void Perturbations::info() const{
  std::cout << "\n";
  std::cout << "Info about perturbations class:\n";
  std::cout << "x_start:       " << x_start                << "\n";
  std::cout << "x_end:         " << x_end                  << "\n";
  std::cout << "n_x:           " << n_x                    << "\n";
  std::cout << "k_min (1/Mpc): " << k_min * Constants.Mpc  << "\n";
  std::cout << "k_max (1/Mpc): " << k_max * Constants.Mpc  << "\n";
  std::cout << "n_k:           " << n_k                    << "\n";
  if(Constants.polarization)
    std::cout << "Polarization has been included.\n";
  else
    std::cout << "Polarization has not been included.\n";
  if(Constants.neutrinos)
    std::cout << "Neutrinos have been included.\n";
  else
    std::cout << "Neutrinos have not been included.\n";

  std::cout << "Information about the perturbation system:\n";
  std::cout << "ind_delta_CDM:      " << Constants.ind_deltacdm         << "\n";
  std::cout << "ind_delta_B:        " << Constants.ind_deltab           << "\n";
  std::cout << "ind_v_CDM:          " << Constants.ind_vcdm             << "\n";
  std::cout << "ind_v_B:            " << Constants.ind_vb               << "\n";
  std::cout << "ind_Phi:            " << Constants.ind_Phi              << "\n";
  std::cout << "ind_start_theta:    " << Constants.ind_start_theta      << "\n";
  std::cout << "n_ell_theta:        " << Constants.n_ell_theta          << "\n";
  if(Constants.polarization){
    std::cout << "ind_start_thetap:   " << Constants.ind_start_thetap   << "\n";
    std::cout << "n_ell_thetap:       " << Constants.n_ell_thetap       << "\n";
  }
  if(Constants.neutrinos){
    std::cout << "ind_start_nu:       " << Constants.ind_start_nu       << "\n";
    std::cout << "n_ell_neutrinos     " << Constants.n_ell_neutrinos    << "\n";
  }
  std::cout << "n_ell_tot_full:     " << Constants.n_ell_tot_full       << "\n";

  std::cout << "Information about the perturbation system in tight coupling:\n";
  std::cout << "ind_deltacdm:       " << Constants.ind_deltacdm_tc      << "\n";
  std::cout << "ind_deltab:         " << Constants.ind_deltab_tc        << "\n";
  std::cout << "ind_v_CDM:          " << Constants.ind_vcdm_tc          << "\n";
  std::cout << "ind_v_B:            " << Constants.ind_vb_tc            << "\n";
  std::cout << "ind_Phi:            " << Constants.ind_Phi_tc           << "\n";
  std::cout << "ind_start_theta:    " << Constants.ind_start_theta_tc   << "\n";
  std::cout << "n_ell_theta:        " << Constants.n_ell_theta_tc       << "\n";
  if(Constants.neutrinos){
    std::cout << "ind_start_nu:       " << Constants.ind_start_nu_tc    << "\n";
    std::cout << "n_ell_neutrinos     " << Constants.n_ell_neutrinos_tc << "\n";
  }
  std::cout << "n_ell_tot_tc:       " << Constants.n_ell_tot_tc         << "\n";
  std::cout << std::endl;
}

//====================================================
// Output some results to file for a given value of k
//====================================================

void Perturbations::output(const double k, const std::string filename) const{
  std::ofstream fp(filename.c_str());
  const int npts = 5000;
  auto x_array = Utils::linspace(x_start, x_end, npts);
  auto print_data = [&] (const double x) {
    double arg = k*(cosmo->eta_of_x(0.0)-cosmo->eta_of_x(x));
    fp << x                                            << " ";
    fp << get_Theta(x,k,0)                             << " ";
    fp << get_Theta(x,k,1)                             << " ";
    fp << get_Theta(x,k,2)                             << " ";
    fp << get_Phi(x,k)                                 << " ";
    fp << get_Psi(x,k)                                 << " ";
    if (Constants.polarization)
    {
      fp << get_Theta_p(x,k,0)                         << " ";
      fp << get_Theta_p(x,k,1)                         << " ";
      fp << get_Theta_p(x,k,2)                         << " ";
    }
    else{    
      fp << 0.0                                        << " ";
      fp << 0.0                                        << " ";
      fp << 0.0                                        << " ";
    }
    if (Constants.neutrinos)
    {
      fp << get_Nu(x,k,0)                              << " ";
      fp << get_Nu(x,k,1)                              << " ";
      fp << get_Nu(x,k,2)                              << " ";
    }
    else{    
      fp << 0.0                                        << " ";
      fp << 0.0                                        << " ";
      fp << 0.0                                        << " ";
    }
    fp << get_delta_B(x, k)                            << " ";
    fp << get_delta_CDM(x, k)                          << " ";
    fp << get_v_B(x, k)                                << " ";
    fp << get_v_CDM(x, k)                              << " ";
    fp << get_Pi(x,k)                                  << " ";
    fp << get_Source_T(x,k)                            << " ";
    fp << get_Source_T(x,k) * Utils::j_ell(5,   arg)   << " ";
    fp << get_Source_T(x,k) * Utils::j_ell(50,  arg)   << " ";
    fp << get_Source_T(x,k) * Utils::j_ell(500, arg)   << " ";
    if (Constants.polarization)
    {
      fp << get_Source_E(x,k)  << " ";
      fp << get_Source_E(x,k) * Utils::j_ell(5,   arg) << " ";
      fp << get_Source_E(x,k) * Utils::j_ell(50,  arg) << " ";
      fp << get_Source_E(x,k) * Utils::j_ell(500, arg) << " ";
    }
    else{    
      fp << 0.0     << " ";
      fp << 0.0     << " ";
      fp << 0.0     << " ";
      fp << 0.0     << " ";
    }
    fp << "\n";
  };
  std::for_each(x_array.begin(), x_array.end(), print_data);
}
