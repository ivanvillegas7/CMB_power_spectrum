#include"PowerSpectrum.h"

//====================================================
// Constructors:
//====================================================

PowerSpectrum::PowerSpectrum(
    BackgroundCosmology *cosmo, 
    RecombinationHistory *rec, 
    Perturbations *pert,
    double A_s,
    double n_s,
    double kpivot_mpc) : 
  cosmo(cosmo), 
  rec(rec), 
  pert(pert),
  A_s(A_s),
  n_s(n_s),
  kpivot_mpc(kpivot_mpc)
{}

//====================================================
// Do all the solving:
//====================================================
void PowerSpectrum::solve(){

  //=========================================================================
  // TODO: Choose the range of k's and the resolution to compute Theta_ell(k).
  //=========================================================================
  double dk          = 2.*M_PI/(eta0*n_k);
  Vector k_array     = Utils::linspace(k_min, k_max, (k_max-k_min)/dk);
  Vector log_k_array = Utils::linspace(log(k_min), log(k_max), 2*(k_max-k_min)/dk);

  //=========================================================================
  // TODO: Make splines for j_ell. 
  // Implement generate_bessel_function_splines.
  //=========================================================================
  generate_bessel_function_splines();
  
  //=========================================================================
  // TODO: Line of sight integration to get Theta_ell(k).
  // Implement line_of_sight_integration.
  //=========================================================================
  line_of_sight_integration(k_array);
  
  //=========================================================================
  // TODO: Integration to get Cell by solving dCell^f/dlogk = Delta(k)*f_ell(k)^2
  // Implement solve_for_cell.
  //=========================================================================
  auto cell_TT = solve_for_cell(log_k_array, thetaT_ell_of_k_spline, thetaT_ell_of_k_spline);
  cell_TT_spline.create(ells, cell_TT, "Cell_TT_of_ell");
  
  //=========================================================================
  // TODO: Do the same for polarization...
  //=========================================================================
  if (Constants.polarization){
    auto cell_TE = solve_for_cell(log_k_array, thetaT_ell_of_k_spline, thetaE_ell_of_k_spline);
    cell_TE_spline.create(ells, cell_TE, "Cell_TE_of_ell");

    auto cell_EE = solve_for_cell(log_k_array, thetaE_ell_of_k_spline, thetaE_ell_of_k_spline);
    cell_EE_spline.create(ells, cell_EE, "Cell_EE_of_ell");
  }
}

//====================================================
// Generate splines of j_ell(z) needed for LOS integration.
//====================================================

void PowerSpectrum::generate_bessel_function_splines(){
  Utils::StartTiming("Bessel spline");

  int const N_ell = ells.size();
  
  // Make storage for the splines
  j_ell_splines = std::vector<Spline>(N_ell);

  // Determine argument interval
  double z_min = 0.0;
  double z_max = k_max*eta0;
  double dz = 2.*M_PI/n_bessel;
  Vector z_array = Utils::linspace(z_min, z_max, (z_max-z_min)/dz);
    
  //=============================================================================
  // TODO: Compute splines for bessel functions j_ell(z).
  // Choose a suitable range for each ell.
  // NB: you don't want to go larger than z ~ 40000, then the bessel routines
  // might break down. Use j_ell(z) = Utils::j_ell(ell, z).
  //=============================================================================

  for(size_t i_l = 0; i_l < N_ell; i_l++){
    const int ell = ells[i_l];

    // Make vector to store j_ell
    Vector j_ell_array(z_array.size());

    // Loop over z-array
    for(int i_z=0; i_z < z_array.size(); i_z++){
      j_ell_array[i_z] = Utils::j_ell(ell, z_array[i_z]);
    }

    // Make the j_ell_splines[i] spline.
    j_ell_splines[i_l].create(z_array, j_ell_array);
  }

  Utils::EndTiming("Bessel spline");
}

//====================================================
// Do the line of sight integration for a single
// source function.
//====================================================

Vector2D PowerSpectrum::line_of_sight_integration_single(
    Vector & k_array, 
    std::function<double(double,double)> &source_function){
  Utils::StartTiming("Line of sight");
  
  // Make storage for the results
  Vector2D result = Vector2D(ells.size(), Vector(k_array.size()));

  // Create arrays
  double dx = 2.*M_PI/n_x_los;
  Vector x_array = Utils::linspace(x_start_los, x_end_los, (x_end_los-x_start_los)/dx);
  
  std::cout << std::endl;
  for(size_t i_k = 0; i_k < k_array.size(); i_k++){

    //=============================================================================
    // TODO: Implement to solve for the general line of sight integral 
    // F_ell(k) = Int dx jell(k(eta-eta0))*S(x,k) for all the ell values 
    // for the given value of k.
    //=============================================================================
    
    // Progress bar...
    if(i_k*10/k_array.size() != ((i_k+1)*10)/k_array.size()){
      std::cout << "Progress: " << (i_k+1)*100/k_array.size() << "%\n" << std::flush;
      if(i_k == k_array.size()-1) std::cout << std::endl;
    }
    double k_value = k_array[i_k]; // k-value for each iteration
    for(int i_l=0; i_l < ells.size(); i_l++){
      //double ell = ells[i_l]; // ell-value for each iteration
      //std::cout<<"Problem is below this line when i_l="<<i_l<<std::endl;
      Vector integrand(x_array.size());
      for(int i=0; i < x_array.size(); i++){
        integrand[i] = source_function(x_array[i], k_value)*j_ell_splines[i_l](k_value*(eta0-cosmo->eta_of_x(x_array[i])));
        //std::cout<<"Problem is below this line when i="<<i<<std::endl;
      }

      // Store the result for Source_ell(k) in results[ell][ik].
      result[i_l][i_k] = integrate(dx, integrand);
    }
  }
  
  Utils::EndTiming("Line of sight");
  return result;
}

//====================================================
// Do the line of sight integration.
//====================================================
void PowerSpectrum::line_of_sight_integration(Vector & k_array){
  const int n_k        = k_array.size();
  const int n          = 100;
  const int nells      = ells.size();
  
  // Make storage for the splines we are to create.
  thetaT_ell_of_k_spline = std::vector<Spline>(nells);

  //============================================================================
  // TODO: Solve for Theta_ell(k) and spline the result.
  //============================================================================

  // Make a function returning the source function.
  std::function<double(double,double)> source_function_T = [&](double x, double k){
    return pert->get_Source_T(x,k);
  };
  
  // Do the line of sight integration.
  Vector2D thetaT_ell_of_k = line_of_sight_integration_single(k_array, source_function_T);
  
  // Spline the result and store it in thetaT_ell_of_k_spline.
  for(int i_l=0; i_l < nells; i_l++){
    thetaT_ell_of_k_spline[i_l].create(k_array, thetaT_ell_of_k[i_l]);
  }
  
  //============================================================================
  // TODO: Solve for ThetaE_ell(k) and spline.
  //============================================================================
  if(Constants.polarization){

    // Make storage for the splines we are to create.
    thetaE_ell_of_k_spline = std::vector<Spline>(nells);

    // Make a function returning the source function.
    std::function<double(double,double)> source_function_E = [&](double x, double k){
    return pert->get_Source_E(x,k);
    };

    // Do the line of sight integration.
    Vector2D thetaE_ell_of_k = line_of_sight_integration_single(k_array, source_function_E);

    // Spline the result and store it in thetaT_ell_of_k_spline.
    for(int i_l=0; i_l < nells; i_l++){
      thetaE_ell_of_k_spline[i_l].create(k_array, thetaE_ell_of_k[i_l]);
    }
  }
}

// Function to integrate.
double PowerSpectrum::integrate(double dx, Vector y_array){
  // Declare and define needed variables.
  double integral_value = 0;
  for(int i=0; i < y_array.size()-1; i++){
    integral_value += (y_array[i+1]+y_array[i])/2;
  }
  return integral_value*dx;
}

//====================================================
// Compute Cell (could be TT or TE or EE):
// Cell = Int_0^inf*4*π*P(k)*f_ell*g_ell*dk/k.
//====================================================
Vector PowerSpectrum::solve_for_cell(
    Vector & log_k_array,
    std::vector<Spline> & f_ell_spline,
    std::vector<Spline> & g_ell_spline){
  const int nells      = ells.size();

  //============================================================================
  // TODO: Integrate Cell = Int*4*π*P(k)*f_ell*g_ell*dk/k, or equivalently
  // solve the ODE system dCell/dlogk = 4*π*P(k)*f_ell*g_ell.
  //============================================================================
  
  Vector result(nells);
  int N        = log_k_array.size();
  double dlogk = (log_k_array[N-1]-log_k_array[0])/N;

  // Loop over and integrate for all ells.
  for(int i_l=0; i_l < nells; i_l++){
    Vector integrand(log_k_array.size());
    for(int i=0; i < log_k_array.size(); i++){
      double k_value = exp(log_k_array[i]);
      integrand[i]   = get_matter_power_spectrum(0.0, k_value)*pow(Constants.Mpc, -2)*f_ell_spline[i_l](k_value)*g_ell_spline[i_l](k_value);
      }

    result[i_l] = 4.*M_PI*integrate(dlogk, integrand);
  }

  return result;
}

//====================================================
// The primordial power-spectrum:
//====================================================

double PowerSpectrum::primordial_power_spectrum(const double k) const{
  return A_s*pow(k*Constants.Mpc/kpivot_mpc, n_s-1.)*2.*pow(M_PI, 2)/pow(k*Constants.Mpc, 3);
}

//====================================================
// P(k) in units of (Mpc)^3:
//====================================================

double PowerSpectrum::get_matter_power_spectrum(const double x, const double k_mpc) const{

  //=============================================================================
  // TODO: Compute the matter power spectrum.
  //=============================================================================

  // Variables:
  double OmegaM  = cosmo->get_OmegaB()+cosmo->get_OmegaCDM();
  double H0      = cosmo->get_H0();
  double Phi     = pert->get_Phi(x, k_mpc);

  // Calculate Delta_M:
  double Delta_M = 2.*pow(Constants.c*k_mpc/H0, 2.)*Phi*exp(x)/(3.*OmegaM);

  // Calculate P(k,x).
  double pofk   = pow(Delta_M, 2)*primordial_power_spectrum(k_mpc);

  return pofk;
}

//====================================================
// Get methods:
//====================================================
double PowerSpectrum::get_cell_TT(const double ell) const{
  return cell_TT_spline(ell);
}
double PowerSpectrum::get_cell_TE(const double ell) const{
  return cell_TE_spline(ell);
}
double PowerSpectrum::get_cell_EE(const double ell) const{
  return cell_EE_spline(ell);
}
double PowerSpectrum::get_ThetaT(const double ell, const double k) const {
  return thetaT_ell_of_k_spline[ell](k);
}
double PowerSpectrum::get_ThetaE(const double ell, const double k) const {
  return thetaE_ell_of_k_spline[ell](k);
}

//====================================================
// Output the cells to file.
//====================================================

void PowerSpectrum::output(std::string filename) const{
  // Output in standard units of μK^2
  std::ofstream fp(filename.c_str());
  const int ellmax = int(ells[ells.size()-1]);
  auto ellvalues   = Utils::linspace(2, ellmax, ellmax-1);
  auto print_data  = [&] (const double ell) {
    double normfactor  = (ell*(ell+1.))/(2.0*M_PI)*pow(1e6*cosmo->get_TCMB(), 2);
    double normfactorN = (ell*(ell+1.))/(2.0*M_PI)*pow(1e6*cosmo->get_TCMB()*pow(4./11., 1./3.), 2);
    double normfactorL = (ell*(ell+1.))*(ell*(ell+1.))/(2.*M_PI);
    fp << ell                               << " ";
    fp << cell_TT_spline(ell)*normfactor    << " ";
    if(Constants.polarization){
      fp << cell_EE_spline(ell)*normfactor  << " ";
      fp << cell_TE_spline(ell)*normfactor  << " ";
    }
    fp << "\n";
  };
  std::for_each(ellvalues.begin(), ellvalues.end(), print_data);
}

void PowerSpectrum::output_MPS(const std::string filename) const{
  
  Vector k_array = Utils::linspace(k_min, k_max, n_k);

  std::ofstream fp(filename.c_str());
  auto print_data = [&] (const double k) {
    fp << k*Constants.Mpc/cosmo->get_h()                        << " ";
    fp << get_matter_power_spectrum(0.0, k)*pow(cosmo->get_h(), 3) << " ";
    fp << k*eta0                                                          << " ";
    fp << get_ThetaT(0, k)                                                << " ";
    fp << get_ThetaT(10, k)                                               << " ";
    fp << get_ThetaT(24, k)                                               << " ";
    fp << get_ThetaT(62, k)                                               << " ";
    fp <<"\n";
  };
  std::for_each(k_array.begin(), k_array.end(), print_data);
}
