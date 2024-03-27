#include"RecombinationHistory.h"

//====================================================
// Constructors
//====================================================
   
RecombinationHistory::RecombinationHistory(
    BackgroundCosmology *cosmo, 
    double Yp) :
  cosmo(cosmo),
  Yp(Yp)
{}

//====================================================
// Do all the solving we need to do
//====================================================

void RecombinationHistory::solve(){
    
  // Compute and spline Xe, ne
  solve_number_density_electrons();
   
  // Compute and spline τ, dτdx, ddτddx, g, dgdx, ddgddx, ...
  solve_for_optical_depth_tau();

  //Compute and spline s
  solve_for_sound_horizon();
}

//====================================================
// Solve for X_e and n_e using Saha and Peebles and spline the result
//====================================================

void RecombinationHistory::solve_number_density_electrons(){
  Utils::StartTiming("Electron fraction (Xe)");
  
  //=============================================================================
  // TODO: Set up x-array and make arrays to store X_e(x) and n_e(x) on
  //=============================================================================
  Vector x_array = Utils::linspace(x_start, x_end, npts_rec_arrays);
  Vector Xe_arr(npts_rec_arrays);
  
  Vector Xe_arr_only_Saha = Xe_arr;
  Vector ne_arr           = Xe_arr;

  // Save index of when leaving Saha regime and start using Peebles
  double Xe_Peebles_transition = 0.0;
  int count = 0;

  // Calculate recombination history
  bool saha_regime = true;

  for(int i = 0; i < npts_rec_arrays; i++){

    //==============================================================
    // TODO: Get X_e from solving the Saha equation so
    // implement the function electron_fraction_from_saha_equation
    //==============================================================
    auto Xe_ne_data = electron_fraction_from_saha_equation(x_array[i]);

    // Electron fraction and number density
    const double Xe_current = Xe_ne_data.first;
    const double ne_current = Xe_ne_data.second;

    // Are we still in the Saha regime?
    if(Xe_current < Xe_saha_limit){

      saha_regime = false;

      if (count==0)
      {
        idx_Peebles_transition = i;
        Xe_Peebles_transition = Xe_current;
        x_Saha_to_Peebles = x_array[i];
        count++;
      }
      
    }

    if(saha_regime){
      
      //=============================================================================
      // TODO: Store the result we got from the Saha equation
      //=============================================================================
      Xe_arr[i]           = Xe_current;
      ne_arr[i]           = ne_current;
      Xe_arr_only_Saha[i] = Xe_current; // Keep the solutions obtained by Saha eq only
    } 
    
    else {

      //==============================================================
      // TODO: Compute X_e from current time until today by solving 
      // the Peebles equation (NB: if you solve all in one go remember
      // to exit the for-loop!)
      // Implement rhs_peebles_ode
      //==============================================================

      // Integrate equation from x(i-1) to x(i) and fetch Xe(i)
      Vector x_array_current{x_array[i-1],x_array[i]};
      double Xe_current = Xe_arr[i-1];

      // The Peebles ODE equation
      ODESolver peebles_Xe_ode;
      ODEFunction dXedx = [&](double x, const double *Xe, double *dXedx){
        return rhs_peebles_ode(x, Xe, dXedx);
      };
      
      //=============================================================================
      // TODO: Set up IC, solve the ODE and fetch the result 
      //=============================================================================

      Vector peebles_ini{Xe_current};
      peebles_Xe_ode.solve(dXedx, x_array_current, peebles_ini);
      auto solution = peebles_Xe_ode.get_data_by_component(0);
      double Xe_now = solution.back();
      Xe_arr[i]=Xe_now;
      ne_arr[i]=Xe_now*get_number_density_H(x_array[i]);
    }
  }

  //=============================================================================
  // TODO: Spline the result. Implement and make sure the Xe_of_x, ne_of_x 
  // functions are working
  //=============================================================================

  // Fill rest of Xe_arr_saha with Saha solution with 0
  for (int i=idx_Peebles_transition; i < npts_rec_arrays; i++){
    auto Xe_ne_data = electron_fraction_from_saha_equation(x_array[i]);
    double const Xe_current = Xe_ne_data.first;
    // Check for nan and negative values and set to zero. 
    double Xe_current_non_neg_zero_nan = Xe_current < 1e-9 || std::isnan(Xe_current) ? 0.0 : Xe_current;
    Xe_arr_only_Saha[i] = Xe_current_non_neg_zero_nan;
    //std::cout<<x_array[i]<<" Only Saha: "<<Xe_arr_only_Saha[i]<<" "<<log(Xe_arr_only_Saha[i])<<" "<<exp(log(Xe_arr_only_Saha[i]))<<"\n";
  }

  // Spline the result in logarithmic form. Used in get Xe_of_x and ne_of_x methods
  Vector log_Xe_arr           = log(Xe_arr);
  Vector log_ne_arr           = log(ne_arr);
  
  log_Xe_of_x_spline.create(x_array,log_Xe_arr,"log Xe");
  Xe_of_x_spline_only_Saha.create(x_array,Xe_arr_only_Saha,"Xe Saha");
  log_ne_of_x_spline.create(x_array,log_ne_arr,"log ne");

  Utils::EndTiming("Electron fraction (Xe)");
}

//====================================================
// Solve the Saha equation to get ne and Xe
//====================================================
std::pair<double,double> RecombinationHistory::electron_fraction_from_saha_equation(double x) const{
  const double a = exp(x);
 
  // Physical constants
  const double c         = Constants.c;
  const double k_b       = Constants.k_b;
  const double G         = Constants.G;
  const double m_e       = Constants.m_e;
  const double hbar      = Constants.hbar;
  const double m_H       = Constants.m_H;
  const double epsilon_0 = Constants.epsilon_0;
  const double H0_over_h = Constants.H0_over_h;

  // Fetch cosmological parameters
  const double T_B = cosmo->get_TCMB()/a;  // Baryon temperature approximation
  const double nH  = get_number_density_H(x);

  //Save FLOPS
  const double E_TB       = T_B*k_b;
  const double eps_over_T = epsilon_0/E_TB;

  // Electron fraction and number density
  double Xe = 0.0;
  
  //=============================================================================
  // TODO: Compute Xe and ne from the Saha equation
  //=============================================================================

  // Right hand side of Saha equation
  const double brakets = m_e*pow(c, 2.)*E_TB/(2.*M_PI); 
  const double rhs_Saha = pow(1./(c*hbar), 3)*brakets*sqrt(brakets)*exp(-eps_over_T);

  // Calculate Xe

  // If near endpoint, take care of instability and set solution to basically zero
  if (rhs_Saha<1e-20) Xe = 1e-20;

  // Determine if we have to use the Taylor approximation in the second order equation
  else if (rhs_Saha>1e+9) Xe = 1.0;

  else {
    if (4.0/rhs_Saha<1e-9) Xe = 1.0;
    
    else Xe = rhs_Saha*(-1.+sqrt(1.0+4.0/rhs_Saha))/2.0;
  }

  // Return electron fraction and number density
  return std::pair<double,double>(Xe, Xe*nH);
}

//====================================================
// The right hand side of the dXedx Peebles ODE
//====================================================
int RecombinationHistory::rhs_peebles_ode(double x, const double *Xe, double *dXedx){

  // Current value of a and X_e
  const double X_e         = Xe[0];
  const double a           = exp(x);

  // Physical constants in SI units
  const double k_b         = Constants.k_b;
  const double c           = Constants.c;
  const double m_e         = Constants.m_e;
  const double hbar        = Constants.hbar;
  const double sigma_T     = Constants.sigma_T;
  const double lambda_2s1s = Constants.lambda_2s1s;
  const double epsilon_0   = Constants.epsilon_0;

  // Cosmological parameters
  const double T_B         = cosmo->get_TCMB()/a;  // Baryon temperature approximation
  const double nH          = get_number_density_H(x);
  const double H           = cosmo->H_of_x(x);

  //=============================================================================
  // TODO: Write the expression for dXedx
  //=============================================================================

  // Some constants to save FLOPS
  const double E_TB       = k_b*T_B;
  const double eps_over_T = epsilon_0/E_TB;

  // Expression needed in dXedx
  const double phi_2        = 0.448*log(eps_over_T);
  const double alpha_2      = 64.*M_PI*(3*sigma_T/(8.*M_PI))*phi_2*sqrt(eps_over_T/(27.*M_PI))*c;
  const double beta_        = alpha_2*pow(m_e*pow(c, 2)*E_TB/(2*M_PI), 3./2.)/pow(c*hbar, 3.);
  const double beta_2       = beta_*exp(-1.*eps_over_T/4.);
  const double n_1s         = (1.0-X_e)*nH;
  const double lambda_alpha = H*pow(3.*epsilon_0/(hbar*c), 3.)/(pow(8.*M_PI, 2.)*n_1s);
  const double C_r          = (lambda_2s1s+lambda_alpha)/(lambda_2s1s+lambda_alpha+beta_2);
  
  dXedx[0] = C_r*(beta_*exp(-eps_over_T)*(1.-X_e)-nH*alpha_2*pow(X_e, 2.))/H;

  return GSL_SUCCESS;
}

//====================================================
// Solve for the optical depth (τ), compute the 
// visibility function and spline the result
//====================================================

void RecombinationHistory::solve_for_optical_depth_tau(){
  Utils::StartTiming("Optical depth (τ)");

  // Set up x-arrays to integrate over. We split into three regions as we need extra points in reionisation
  const int npts = npts_rec_arrays;
  Vector x_array = Utils::linspace(x_start, x_end, npts);

  // The ODE system dτ/dx, dτ_noreion/dx and dτ_baryon/dx
  ODEFunction dtaudx = [&](double x, const double *tau, double *dtaudx){

    // Set the derivative for photon optical depth
    dtaudx[0] = Constants.c*ne_of_x(-x)*Constants.sigma_T/cosmo->H_of_x(-x);

    return GSL_SUCCESS;
  };

  //=============================================================================
  // TODO: Set up and solve the ODE and make τ splines
  //=============================================================================
  
  Vector tau_ini{0.0};

  Vector x_arr_reversed = Utils::linspace(-x_end, -x_start, npts_rec_arrays);

  ODESolver ode;
  ode.solve(dtaudx, x_arr_reversed, tau_ini);

  auto tau_array_reversed = ode.get_data_by_component(0);

  // Reverse the arrays to get correct order in elements 
  Vector x_array_tau(npts_rec_arrays);
  Vector tau_arr(npts_rec_arrays);
  for (int i=0; i < npts_rec_arrays; i++){
    x_array_tau[i] = -x_arr_reversed[npts_rec_arrays-1-i];
    tau_arr[i]     = tau_array_reversed[npts_rec_arrays-1-i];
  }
  
  tau_of_x_spline.create(x_array_tau, tau_arr, "τ");

  //=============================================================================
  // TODO: Compute visibility functions and spline everything
  //=============================================================================

  Vector g_tilde_array(npts_rec_arrays, 0.0);

  for (int i = 0; i < npts_rec_arrays; i++) g_tilde_array[i]=-dtaudx_of_x(x_array[i])*exp(-tau_of_x(x_array[i]));

  g_tilde_of_x_spline.create(x_array, g_tilde_array, "g̃");

  Utils::EndTiming("Optical depth (τ)");
}

void RecombinationHistory::solve_for_sound_horizon(){
  Utils::StartTiming("Sound horizon (s)");

  Vector x_array = Utils::linspace(x_start, x_end, npts_rec_arrays);

  // Set up ODE for sound horizon
  ODESolver sound_horizon_ode;
  ODEFunction dsdx = [&](double x, const double *s, double *dsdx){
    dsdx[0] = cs_of_x(x)/cosmo->Hp_of_x(x);
    return GSL_SUCCESS;
  };

  Vector s_ini{cs_of_x(x_start)/cosmo->Hp_of_x(x_start)};
  sound_horizon_ode.solve(dsdx, x_array, s_ini);
  auto s_array = sound_horizon_ode.get_data_by_component(0);

  // Spline result 
  sound_horizon_of_x_spline.create(x_array, s_array, "s");

  Utils::EndTiming("Sound horizon (s)");
}

//====================================================
// Get methods
//====================================================

double RecombinationHistory::get_number_density_H(double x) const{
  
  return (1.-Yp)*3.*pow(cosmo->get_H0(), 2.)*cosmo->get_OmegaB(0.0)/(8.*M_PI*Constants.G*Constants.m_H*exp(3.*x));
}

double RecombinationHistory::tau_of_x(double x) const{

  return tau_of_x_spline(x);
}

double RecombinationHistory::dtaudx_of_x(double x) const{

  return tau_of_x_spline.deriv_x(x);
}

double RecombinationHistory::ddtauddx_of_x(double x) const{

  return tau_of_x_spline.deriv_xx(x);
}

double RecombinationHistory::g_tilde_of_x(double x) const{

  return g_tilde_of_x_spline(x);
}

double RecombinationHistory::dgdx_tilde_of_x(double x) const{

  return g_tilde_of_x_spline.deriv_x(x);
}

double RecombinationHistory::ddgddx_tilde_of_x(double x) const{

  return g_tilde_of_x_spline.deriv_xx(x);
}

double RecombinationHistory::Xe_of_x(double x) const{

  return exp(log_Xe_of_x_spline(x));
}

double RecombinationHistory::Xe_of_x_Saha_approx(double x) const{
  return Xe_of_x_spline_only_Saha(x);
}

double RecombinationHistory::ne_of_x(double x) const{

  return exp(log_ne_of_x_spline(x));
}

double RecombinationHistory::get_Yp() const{

  return Yp;
}

double RecombinationHistory::R_of_x(double x) const{
  return 4.*cosmo->get_OmegaR(0.0)/(3.*cosmo->get_OmegaB(0.0));
}

double RecombinationHistory::cs_of_x(double x) const{
  double R=R_of_x(x);
  return Constants.c*sqrt(R/(3.*(1.+R)));
}
double RecombinationHistory::sound_horizon_of_x(double x) const{
  return sound_horizon_of_x_spline(x);
}

//====================================================
// Print some useful info about the class
//====================================================
void RecombinationHistory::info() const{
  std::cout<<"\n";
  std::cout<<"Info about recombination/reionization history class:\n";
  std::cout<<"Yp: "<<Yp<<"\n";
  std::cout<<std::endl;
}
void RecombinationHistory::sound_horizon() const{
  std::cout<<"Sound horizon at decoupling:\n";
  std::cout<<"r_s≡s(x_decoupling)="<<sound_horizon_of_x(x_Saha_to_Peebles)/Constants.Mpc<<" Mpc"<<"\n";
  std::cout<<std::endl;
} 

//====================================================
// Output the data computed to file
//====================================================
void RecombinationHistory::output(const std::string filename) const{
  std::ofstream fp(filename.c_str());

  Vector x_array = Utils::linspace(x_start, x_end, npts_rec_arrays);
  auto print_data = [&] (const double x) {
    fp << x                      << " ";
    fp << Xe_of_x(x)             << " ";
    fp << ne_of_x(x)             << " ";
    fp << tau_of_x(x)            << " ";
    fp << dtaudx_of_x(x)         << " ";
    fp << ddtauddx_of_x(x)       << " ";
    fp << g_tilde_of_x(x)        << " ";
    fp << dgdx_tilde_of_x(x)     << " ";
    fp << ddgddx_tilde_of_x(x)   << " ";
    fp << Xe_of_x_Saha_approx(x) << " ";
    fp << "\n";
  };
  std::for_each(x_array.begin(), x_array.end(), print_data);
}

