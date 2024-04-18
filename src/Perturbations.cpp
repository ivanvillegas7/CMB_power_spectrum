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

  // Loop over all wavenumbers
  for(int ik = 0; ik < n_k; ik++){

    // Progress bar...
    if( (10*ik) / n_k != (10*ik+10) / n_k ) {
      std::cout << (100*ik+100)/n_k << "% " << std::flush;
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
    // ...
    // ...
    // ...
    // ...
    // ...

    //====i===============================================================
    // TODO: Full equation integration
    // Remember to implement the routines:
    // set_ic_after_tight_coupling : The IC after tight coupling ends
    // rhs_full_ode : The dydx for our coupled ODE system
    //===================================================================

    // Set up initial conditions (y_tight_coupling is the solution at the end of tight coupling)
    // auto y_full_ini = set_ic_after_tight_coupling(y_tight_coupling, x_end_tight, k);

    // The full ODE system
    ODEFunction dydx_full = [&](double x, const double *y, double *dydx){
      return rhs_full_ode(x, k, y, dydx);
    };

    // Integrate from x_end_tight -> x_end
    // ...
    // ...
    // ...
    // ...
    // ...

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
    //
    //===================================================================
    //...
    //...

  }
  Utils::EndTiming("Integrate perturbation");

  //=============================================================================
  // TODO: Make all splines needed: Θ_0,Θ_1,Θ_2,Φ,Ψ,...
  //=============================================================================
  // ...
  // ...
  // ...
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
  double *Nu           = &y_tc[Constants.ind_start_nu_tc];

  //=============================================================================
  // TODO: Set the initial conditions in the tight coupling regime
  //=============================================================================

  // Set constants
  double c                    = Constants.c;

  // Cosmological parameters and variables
  double H0                   = cosmo->get_H0();
  double Hp                   = cosmo->Hp_of_x(x);
  double Omega_CDM            = cosmo->get_OmegaCDM(x);//();
  double Omega_B              = cosmo->get_OmegaB(x);//();
  double Omega_R              = cosmo->get_OmegaR(x);//();
  double Omega_Nu             = cosmo->get_OmegaNu(x);//();

  //Set recombination parameters
  double R                    = rec->R_of_x(x);//();
  double dtaudx               = rec->dtaudx_of_x(x);

  // Save flops
  double ck_over_Hp          = c*k/Hp;
  double f_Nu                 = Omega_Nu/(Omega_R+Omega_Nu);

  // SET: Scalar quantities (Gravitational potental, baryons and CDM)
  double Psi                  = -1/(3/2+2*f_Nu/5);
  delta_CDM                   = -3*Psi/2;
  delta_B                     = -3*Psi/2;
  v_CDM                       = -ck_over_Hp*Psi/2;
  v_B                         = -ck_over_Hp*Psi/2;
  Phi                         = -(1+2*f_Nu/5);

  // SET: Photon temperature perturbations (Theta_ell)
  Theta[0]                   = -Psi/2;
  Theta[1]                   = ck_over_Hp*Psi/6;
  if (polarization) Theta[2] = -(8/15)*ck_over_Hp*Theta[1]/dtaudx;
  else Theta[2]              = -(8/15)*ck_over_Hp*Theta[1]/dtaudx;
  for (int l = 3; l <= n_ell_theta_tc; l++)
  {
    Theta[l]                 = -l*ck_over_Hp*Theta[l-1]/((2*l+1)*dtaudx);
  }

  // SET: Neutrino perturbations (N_ell)
  if(neutrinos){
    Nu[0]                    = -Psi/2;
    Nu[1]                    = ck_over_Hp*Psi/6;
    Nu[2]                    = (Psi+Phi)*pow(c*k*exp(x)/H0, 2)/(12*Omega_Nu);
    for (int l = 3; l <= n_ell_neutrinos_tc; l++)
    {
      Nu[l]                  = ck_over_Hp*Nu[l-1]/(2*l+1);
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
  const double *Nu_tc           = &y_tc[Constants.ind_start_nu_tc];

  // References to the quantities we are going to set
  double &delta_CDM       =  y[Constants.ind_deltacdm_tc];
  double &delta_B         =  y[Constants.ind_deltab_tc];
  double &v_CDM           =  y[Constants.ind_vcdm_tc];
  double &v_B             =  y[Constants.ind_vb_tc];
  double &Phi             =  y[Constants.ind_Phi_tc];
  double *Theta           = &y[Constants.ind_start_theta_tc];
  double *Theta_p         = &y[Constants.ind_start_thetap_tc];
  double *Nu              = &y[Constants.ind_start_nu_tc];

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
  double Omega_CDM            = cosmo->get_OmegaCDM(x);//();
  double Omega_B              = cosmo->get_OmegaB(x);//();
  double Omega_R              = cosmo->get_OmegaR(x);//();
  double Omega_Nu             = cosmo->get_OmegaNu(x);//();

  //Set recombination parameters
  double R                    = rec->R_of_x(x);//();
  double dtaudx               = rec->dtaudx_of_x(x);

  // Save flops
  double ck_over_Hp          = c*k/Hp;
  double f_Nu                 = Omega_Nu/(Omega_R+Omega_Nu);

  // SET: Scalar quantities (Gravitational potental, baryons and CDM)
  double Psi                  = -1/(3/2+2*f_Nu/5);
  delta_CDM                   = -3*Psi/2;//delta_CMB_tc;
  delta_B                     = -3*Psi/2;//delta_B_tc;
  v_CDM                       = -ck_over_Hp*Psi/2;//v_CDM_tc;
  v_B                         = -ck_over_Hp*Psi/2;//v_B_tc;
  Phi                         = -(1+2*f_Nu/5);//Phi_tc;

  // SET: Photon temperature perturbations (Theta_ell)
  Theta[0]                   = -Psi/2;//Theta_tc[0];
  Theta[1]                   = ck_over_Hp*Psi/6;//Theta_tc[1];
  if (polarization) Theta[2] = -(8/15)*ck_over_Hp*Theta[1]/dtaudx;//Theta_tc[2];
  else Theta[2]              = -(8/15)*ck_over_Hp*Theta[1]/dtaudx;//Theta_tc[2];
  for (int l = 3; l <= n_ell_theta; l++)
  {
    Theta[l]                 = -l*ck_over_Hp*Theta[l-1]/((2*l+1)*dtaudx);//Theta_tc[l];
  }

  // SET: Photon polarization perturbations (Theta_p_ell)
  if(polarization){
    Theta_p[0]               = 5*Theta[2]/4;
    Theta_p[1]               = ck_over_Hp*Theta[2]/(4*dtaudx);
    Theta_p[2]               = Theta[2]/4;
    for (int l = 3; l <= n_ell_thetap; l++)
    {
      Theta_p[l]             = -l*ck_over_Hp*Theta_p[l-1]/((2*l+1)*dtaudx);
    }
  }

  // SET: Neutrino perturbations (N_ell)
  if(neutrinos){
    Nu[0]                    = -Psi/2;//Nu_tc[0];
    Nu[1]                    = ck_over_Hp*Psi/6;//Nu_tc[1];
    Nu[2]                    = (Psi+Phi)*pow(c*k*exp(x)/H0, 2)/(12*Omega_Nu);//Nu_tc[2];
    for (int l = 3; l <= n_ell_neutrinos; l++)
    {
      Nu[l]                  = ck_over_Hp*Nu[l-1]/(2*l+1);//Nu_tc[l];
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
  // ...
  // ...

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
  Vector ST_array(k_array.size() * x_array.size());
  Vector SE_array(k_array.size() * x_array.size());

  // Compute source functions
  for(auto ix = 0; ix < x_array.size(); ix++){
    const double x = x_array[ix];
    for(auto ik = 0; ik < k_array.size(); ik++){
      const double k = k_array[ik];

      // NB: This is the format the data needs to be stored 
      // in a 1D array for the 2D spline routine source(ix,ik) -> S_array[ix + nx * ik]
      const int index = ix + n_x * ik;

      //=============================================================================
      // TODO: Compute the source functions
      //=============================================================================
      // Fetch all the things we need...
      // const double Hp       = cosmo->Hp_of_x(x);
      // const double tau      = rec->tau_of_x(x);
      // ...
      // ...

      // Temperature source
      ST_array[index] = 0.0;

      // Polarization source
      if(Constants.polarization){
        SE_array[index] = 0.0;
      }
    }
  }

  // Spline the source functions
  ST_spline.create (x_array, k_array, ST_array, "Source_Temp_x_k");
  if(Constants.polarization){
    SE_spline.create (x_array, k_array, SE_array, "Source_Pol_x_k");
  }

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
  double Omega_CDM              = cosmo->get_OmegaCDM(x);//();
  double Omega_B                = cosmo->get_OmegaB(x);//();
  double Omega_R                = cosmo->get_OmegaR(x);//();
  double Omega_Nu               = cosmo->get_OmegaNu(x);//();

  //Set recombination parameters
  double R                      = rec->R_of_x(x);//();
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
  const double *Nu              = &y[Constants.ind_start_nu_tc];

  // References to the quantities we are going to set in the dydx array
  double &ddelta_CDMdx          =  dydx[Constants.ind_deltacdm_tc];
  double &ddelta_Bdx            =  dydx[Constants.ind_deltab_tc];
  double &dv_CDMdx              =  dydx[Constants.ind_vcdm_tc];
  double &dv_Bdx                =  dydx[Constants.ind_vb_tc];
  double &dPhidx                =  dydx[Constants.ind_Phi_tc];
  double *dThetadx              = &dydx[Constants.ind_start_theta_tc];
  double *dNudx                 = &dydx[Constants.ind_start_nu_tc];

  // Save flops
  double ck_over_Hp             = c*k/Hp;
  double dHp_over_Hp            = dHpdx/Hp;

  //=============================================================================
  // TODO: fill in the expressions for all the derivatives
  //=============================================================================

  // SET: Scalar quantities (Φ, δ, v, ...)
  double Psi   = -Phi-12*pow(H0/(c*k*exp(x)), 2)*(Omega_R*Theta[2]+Omega_Nu*Nu[2]);
  double q     = -(((1-R)*dtaudx+(1+R)*ddtauddx)*(3*Theta[1]+v_B)-ck_over_Hp*Psi+(1-dHp_over_Hp)*ck_over_Hp*(-Theta[0]+2+Theta[2])-ck_over_Hp*dThetadx[0])/((1+R)*dtaudx+dHp_over_Hp-1);
  
  ddelta_Bdx   = ck_over_Hp*v_B-3*dPhidx;
  dv_Bdx       = (-v_B-ck_over_Hp*Psi+R*(q+ck_over_Hp*(-Theta[0]+2*Theta[2])-ck_over_Hp*Psi))/(1+R);
  ddelta_CDMdx = ck_over_Hp*v_CDM-3*dPhidx;
  dv_CDMdx     = -v_CDM-ck_over_Hp*Psi;
  dPhidx       = Psi-pow(ck_over_Hp, 2)/3+pow(H0/Hp, 2)*exp(-x)*(Omega_CDM*delta_CDM+Omega_B*delta_B*4*Omega_R*exp(-x)*Theta[0]+4*Omega_Nu*exp(-x)*Nu[0]);

  // SET: Photon multipoles (Θ_ell)
  dThetadx[0] = -ck_over_Hp*Theta[1]-dPhidx;
  dThetadx[1] = (q-dv_Bdx)/3;

  // SET: Neutrino mutlipoles (Nu_ell)
  if(neutrinos){
    dNudx[0]  = -ck_over_Hp*Nu[1]-dPhidx;
    dNudx[1]  = -ck_over_Hp*Nu[0]/3-(2/3)*ck_over_Hp*Nu[2]+ck_over_Hp*Psi/3;
    for (int l = 2; l < n_ell_neutrinos_tc; l++) dNudx[l] = l*ck_over_Hp*Nu[l-1]/(2*l-1)-(l+1)*ck_over_Hp*Nu[l+1]/(2l+1);
    dNudx[n_ell_neutrinos_tc] = ck_over_Hp*Nu[n_ell_neutrinos_tc-1]-c*(n_ell_neutrinos_tc+1)*Nu[n_ell_neutrinos_tc]/(Hp*cosmo->eta_of_x(x));
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
  const double *Theta_p     = &y[Constants.ind_start_thetap];
  const double *Nu          = &y[Constants.ind_start_nu];

  // References to the quantities we are going to set in the dydx array
  double &ddelta_CDMdx      =  dydx[Constants.ind_deltacdm];
  double &ddelta_Bdx        =  dydx[Constants.ind_deltab];
  double &dv_CDMdx          =  dydx[Constants.ind_vcdm];
  double &dv_Bdx            =  dydx[Constants.ind_vb];
  double &dPhidx            =  dydx[Constants.ind_Phi];
  double *dThetadx          = &dydx[Constants.ind_start_theta];
  double *dTheta_pdx        = &dydx[Constants.ind_start_thetap];
  double *dNudx             = &dydx[Constants.ind_start_nu];

  // Set constants
  double c                  = Constants.c;

  // Cosmological parameters and variables
  double H0                 = cosmo->get_H0();
  double Hp                 = cosmo->Hp_of_x(x);
  double Omega_CDM          = cosmo->get_OmegaCDM(x);//();
  double Omega_B            = cosmo->get_OmegaB(x);//();
  double Omega_R            = cosmo->get_OmegaR(x);//();
  double Omega_Nu           = cosmo->get_OmegaNu(x);//();
  double eta                = cosmo->eta_of_x(x);

  //Set recombination parameters
  double R                  = rec->R_of_x(x);//();
  double dtaudx             = rec->dtaudx_of_x(x);

  // Save flops
  double ck_over_Hp         = c*k/Hp;
  double pi                 = Theta[2]+Theta_p[0]+Theta_p[2];

  //=============================================================================
  // TODO: fill in the expressions for all the derivatives
  //=============================================================================

  // SET: Scalar quantities (Φ, δ, v, ...)
  double Psi              = -Phi-12*pow(H0/(c*k*exp(x)), 2)*(Omega_R*Theta[2]+Omega_Nu*Nu[2]);
  dPhidx                  = Psi-pow(ck_over_Hp, 2)/3+pow(H0/Hp, 2)*exp(-x)*(Omega_CDM*delta_CDM+Omega_B*delta_B*4*Omega_R*exp(-x)*Theta[0]+4*Omega_Nu*exp(-x)*Nu[0]);
  dv_Bdx                  = -v_B-ck_over_Hp*Psi+dtaudx*R*(3*Theta[1]+v_B);
  ddelta_Bdx              = ck_over_Hp*v_B-3*dPhidx;
  dv_CDMdx                = -v_CDM-ck_over_Hp*Psi;
  ddelta_CDMdx            = ck_over_Hp*v_CDM-3*dPhidx;

  // SET: Photon multipoles (Θ_ell)
  dThetadx[0]             = -ck_over_Hp*Theta[1]-dPhidx;
  dThetadx[1]             = ck_over_Hp+Theta[0]/3-(2/3)*ck_over_Hp*Theta[2]+ck_over_Hp*Psi/3+dtaudx*(Theta[1]+v_B/3);
  for (int l = 2; l < n_ell_theta; l++)
  {
    if (l==2) dThetadx[l] = l*ck_over_Hp*Theta[l-1]/(2*l+1)-(l+1)*ck_over_Hp*Theta[l+1]/(2*l+1)+dtaudx*(Theta[l]-pi/10);
    else dThetadx[l]      = l*ck_over_Hp*Theta[l-1]/(2*l+1)-(l+1)*ck_over_Hp*Theta[l+1]/(2*l+1)+dtaudx*Theta[l];    
  }
  dThetadx[n_ell_theta+1] = ck_over_Hp*Theta[n_ell_theta]+(dtaudx-c*(n_ell_theta+1)/(Hp*eta))*Theta[n_ell_theta];

  // SET: Photon polarization multipoles (Θ_p_ell)
  if(polarization){
    dTheta_pdx[0]          = -ck_over_Hp*Theta_p[1]+dtaudx*(Theta_p[0]-pi/2);
    for (int l = 2; l < n_ell_thetap; l++)
  {
    if (l==2){
      dTheta_pdx[l]        = l*ck_over_Hp*Theta_p[l-1]/(2*l+1)-(l+1)*ck_over_Hp*Theta_p[l+1]/(2*l+1)+dtaudx*(Theta_p[l]-pi/10);
    }
    else{
      dTheta_pdx[l]        = l*ck_over_Hp*Theta_p[l-1]/(2*l+1)-(l+1)*ck_over_Hp*Theta_p[l+1]/(2*l+1)+dtaudx*Theta_p[l];
    }
  }
  dThetadx[n_ell_thetap+1] = ck_over_Hp*Theta_p[n_ell_thetap]+(dtaudx-c*(n_ell_thetap+1)/(Hp*eta))*Theta_p[n_ell_thetap];
  }

  // SET: Neutrino mutlipoles (Nu_ell)
  if(neutrinos){
    dNudx[0]               = -ck_over_Hp*Nu[1]-dPhidx;
    dNudx[1]               = -ck_over_Hp*Nu[0]/3-(2/3)*ck_over_Hp*Nu[2]+ck_over_Hp*Psi/3;
    for (int l = 2; l < n_ell_neutrinos; l++){
      dNudx[l]             = l*ck_over_Hp*Nu[l-1]/(2*l-1)-(l+1)*ck_over_Hp*Nu[l+1]/(2l+1);
    }
    dNudx[n_ell_neutrinos] = ck_over_Hp*Nu[n_ell_neutrinos-1]-c*(n_ell_neutrinos+1)*Nu[n_ell_neutrinos]/(Hp*cosmo->eta_of_x(x));
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
double Perturbations::get_Phi(const double x, const double k) const{
  return Phi_spline(x,k);
}
double Perturbations::get_Psi(const double x, const double k) const{
  return Psi_spline(x,k);
}
double Perturbations::get_Pi(const double x, const double k) const{
  return Pi_spline(x,k);
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
double Perturbations::get_Theta_p(const double x, const double k, const int ell) const{
  return Theta_p_spline[ell](x,k);
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
  std::cout << "n_x:     " << n_x              << "\n";
  std::cout << "k_min (1/Mpc): " << k_min * Constants.Mpc  << "\n";
  std::cout << "k_max (1/Mpc): " << k_max * Constants.Mpc  << "\n";
  std::cout << "n_k:     " << n_k              << "\n";
  if(Constants.polarization)
    std::cout << "We include polarization.\n";
  else
    std::cout << "We do not include polarization.\n";
  if(Constants.neutrinos)
    std::cout << "We include neutrinos.\n";
  else
    std::cout << "We do not include neutrinos.\n";

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
    double arg = k * (cosmo->eta_of_x(0.0) - cosmo->eta_of_x(x));
    fp << x                  << " ";
    fp << get_Theta(x,k,0)   << " ";
    fp << get_Theta(x,k,1)   << " ";
    fp << get_Theta(x,k,2)   << " ";
    fp << get_Phi(x,k)       << " ";
    fp << get_Psi(x,k)       << " ";
    fp << get_Pi(x,k)        << " ";
    fp << get_Source_T(x,k)  << " ";
    fp << get_Source_T(x,k) * Utils::j_ell(5,   arg)           << " ";
    fp << get_Source_T(x,k) * Utils::j_ell(50,  arg)           << " ";
    fp << get_Source_T(x,k) * Utils::j_ell(500, arg)           << " ";
    fp << "\n";
  };
  std::for_each(x_array.begin(), x_array.end(), print_data);
}

