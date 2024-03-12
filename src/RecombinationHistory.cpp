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
   
  // Compute and spline tau, dtaudx, ddtauddx, g, dgdx, ddgddx, ...
  solve_for_optical_depth_tau();
}

//====================================================
// Solve for X_e and n_e using Saha and Peebles and spline the result
//====================================================

void RecombinationHistory::solve_number_density_electrons(){
  Utils::StartTiming("Xe");
  
  //=============================================================================
  // TODO: Set up x-array and make arrays to store X_e(x) and n_e(x) on
  //=============================================================================
  Vector x_array = Utils::linspace(x_start, x_end, npts_rec_arrays);
  Vector Xe_arr(npts_rec_arrays);
  
  Vector Xe_arr_only_Saha = Xe_arr;
  Vector ne_arr           = Xe_arr;

  // Save index of when leaving Saha regime and start using Peebles
  int idx_Peebles_transition   = 0;
  double Xe_Peebles_transition = 0;

  // Calculate recombination history
  bool saha_regime = true;

  // Used later if we want to not compute the solution using Saha all the way in later milestones
  bool break_Saha_solution = false;

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
    if(Xe_current < Xe_saha_limit)
      saha_regime = false;

    if(saha_regime){
      
      //=============================================================================
      // TODO: Store the result we got from the Saha equation
      //=============================================================================

      Xe_arr[i]           = Xe_current;
      ne_arr[i]           = ne_current;
      Xe_arr_only_Saha[i] = Xe_current; // Keep the solutions obtained by Saha eq only

    } else {

      //==============================================================
      // TODO: Compute X_e from current time until today by solving 
      // the Peebles equation (NB: if you solve all in one go remember to
      // exit the for-loop!)
      // Implement rhs_peebles_ode
      //==============================================================
      //...
      //...

      // The Peebles ODE equation
      ODESolver peebles_Xe_ode;
      ODEFunction dXedx = [&](double x, const double *Xe, double *dXedx){
        return rhs_peebles_ode(x, Xe, dXedx);
      };
      
      //=============================================================================
      // TODO: Set up IC, solve the ODE and fetch the result 
      //=============================================================================
      //...
      //...
    
    }
  }

  //=============================================================================
  // TODO: Spline the result. Implement and make sure the Xe_of_x, ne_of_x 
  // functions are working
  //=============================================================================
  
  // Spline the result in logarithmic form. Used in get Xe_of_x and ne_of_x methods
  Vector log_Xe_arr           = log(Xe_arr);
  Vector log_Xe_arr_only_Saha = log(Xe_arr_only_Saha);
  Vector log_ne_arr           = log(ne_arr);

  log_Xe_of_x_spline.create(x_array,log_Xe_arr,"log Xe");
  log_ne_of_x_spline.create(x_array,log_ne_arr,"log ne");

  Utils::EndTiming("Xe");
}

//====================================================
// Solve the Saha equation to get ne and Xe
//====================================================
std::pair<double,double> RecombinationHistory::electron_fraction_from_saha_equation(double x) const{
  const double a = exp(x);
 
  // Physical constants
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

  // Electron fraction and number density
  double Xe = 0.0;
  double ne = 0.0;
  
  //=============================================================================
  // TODO: Compute Xe and ne from the Saha equation
  //=============================================================================

  // Right hand side of Saha equation
  const double rhs_Saha = pow(m_e*T_B/(2*M_PI), 3./2.)*exp(-epsilon_0/(k_b*T_B))/nH;

  // Calculate Xe

  // If near endpoint, take care of instability and set solution to basically zero
  if (rhs_Saha < 1e-20) Xe = 1e-20;

  // Determine if we have to use the Taylor approximation in the second order equation
  else if (rhs_Saha > 1e+9) Xe = 1.0;

  else Xe = rhs_Saha*(-1 + sqrt(1.0+4.0/rhs_Saha))/2.0;

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
  const double m_H         = Constants.m_H;
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
  const double hbar_square = pow(hbar, 2);
  const double E_B = k_b*T_B;

  // Expression needed in dXedx
  const double alpha        = 1./137.;
  const double phi_2        = 0.448*log(epsilon_0/E_B);
  const double alpha_2      = 64*M_PI*pow(alpha/m_e, 2)*phi_2*sqrt(epsilon_0/(27*M_PI*E_B));
  const double beta         = alpha_2*E_B*pow(m_e*E_B/(2*M_PI), 3./2.)*exp(-epsilon_0/E_B);
  const double beta_2       = beta*exp(3*epsilon_0/(4*E_B));
  const double n_1s         = (1.0-X_e)*nH;
  const double lambda_alpha = H*pow(3*epsilon_0, 2)/(pow(8*M_PI, 2)*n_1s*hbar_square*hbar*pow(c, 3));
  const double C_r          = (lambda_2s1s+lambda_alpha)/(lambda_2s1s+lambda_alpha+beta_2);
  
  dXedx[0] = C_r/H*(beta*(1-X_e)-nH*alpha_2*pow(X_e, 2));

  return GSL_SUCCESS;
}

//====================================================
// Solve for the optical depth tau, compute the 
// visibility function and spline the result
//====================================================

void RecombinationHistory::solve_for_optical_depth_tau(){
  Utils::StartTiming("opticaldepth");

  // Set up x-arrays to integrate over. We split into three regions as we need extra points in reionisation
  const int npts = 1000;
  Vector x_array = Utils::linspace(x_start, x_end, npts);

  // The ODE system dtau/dx, dtau_noreion/dx and dtau_baryon/dx
  ODEFunction dtaudx = [&](double x, const double *tau, double *dtaudx){

    // Set the derivative for photon optical depth
    dtaudx[0] = -Constants.c*ne_of_x(x)*Constants.sigma_T/cosmo->H_of_x(x);

    return GSL_SUCCESS;
  };

  //=============================================================================
  // TODO: Set up and solve the ODE and make tau splines
  //=============================================================================
  //...
  //...

  //=============================================================================
  // TODO: Compute visibility functions and spline everything
  //=============================================================================
  //...
  //...

  Utils::EndTiming("opticaldepth");
}

//====================================================
// Get methods
//====================================================

double RecombinationHistory::get_number_density_H(double x) const{
  
  return (1-Yp)*3*pow(cosmo->get_H0(), 2)*cosmo->get_OmegaB(0.0)/(8*M_PI*Constants.G*Constants.m_H*exp(-3*x));
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

double RecombinationHistory::ne_of_x(double x) const{

  return exp(log_ne_of_x_spline(x));
}

double RecombinationHistory::get_Yp() const{
  return Yp;
}

//====================================================
// Print some useful info about the class
//====================================================
void RecombinationHistory::info() const{
  std::cout << "\n";
  std::cout << "Info about recombination/reionization history class:\n";
  std::cout << "Yp:          " << Yp          << "\n";
  std::cout << std::endl;
} 

//====================================================
// Output the data computed to file
//====================================================
void RecombinationHistory::output(const std::string filename) const{
  std::ofstream fp(filename.c_str());
  const int npts       = 5000;
  const double x_min   = x_start;
  const double x_max   = x_end;

  Vector x_array = Utils::linspace(x_min, x_max, npts);
  auto print_data = [&] (const double x) {
    fp << x                    << " ";
    fp << Xe_of_x(x)           << " ";
    fp << ne_of_x(x)           << " ";
    fp << tau_of_x(x)          << " ";
    fp << dtaudx_of_x(x)       << " ";
    fp << ddtauddx_of_x(x)     << " ";
    fp << g_tilde_of_x(x)      << " ";
    fp << dgdx_tilde_of_x(x)   << " ";
    fp << ddgddx_tilde_of_x(x) << " ";
    fp << "\n";
  };
  std::for_each(x_array.begin(), x_array.end(), print_data);
}

