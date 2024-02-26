#include "BackgroundCosmology.h"
#include <cmath>

//====================================================
// Constructors
//====================================================
    
BackgroundCosmology::BackgroundCosmology(
    double h, 
    double OmegaB, 
    double OmegaCDM, 
    double OmegaK,
    double Neff, 
    double TCMB) :
  h(h),
  OmegaB(OmegaB),
  OmegaCDM(OmegaCDM),
  OmegaK(OmegaK),
  Neff(Neff), 
  TCMB(TCMB)
{

  //=============================================================================
  // TODO: Compute OmegaR, OmegaNu, OmegaLambda, H0, ...
  //=============================================================================

  const double pi     = M_PI;
  const double k_B    = Constants.k_b;                                                                   // J/K
  const double h_bar  = Constants.hbar;                                                                  // J*s
  const double c      = Constants.c;                                                                     // m/s
  const double G      = Constants.G;                                                                     // m³/(kg*s²)

  H0     = Constants.H0_over_h*h;                                                                        // The Hubble parameter today in 1/s

  OmegaR       = 2*(pow(pi, 2)/30)*(pow(k_B*TCMB, 4)/(pow(h_bar, 3)*pow(c, 5)))*(8*pi*G/(3*pow(H0, 2))); // Photon density today
  OmegaNu      = Neff*(7./8.)*pow(4./11., 4./3.)*OmegaR;                                                                                      // Neutrino density today
  OmegaLambda  = 1-(OmegaB+OmegaR+OmegaNu+OmegaCDM+OmegaK);                                              // Dark energy (Λ) density

}

//====================================================
// Do all the solving. Compute eta(x)
//====================================================

// Solve the background
void BackgroundCosmology::solve(){

  Utils::StartTiming("Eta");
    
  //=============================================================================
  // TODO: Set the range of x and the number of points for the splines
  // For this Utils::linspace(x_start, x_end, npts) is useful
  //=============================================================================

  double npts = 1000;

  Vector x_array = Utils::linspace(x_start, x_end, npts);

  // The ODE for deta/dx
  ODEFunction detadx = [&](double x, const double *eta, double *detadx){

    //=============================================================================
    // TODO: Set the rhs of the detadx ODE
    //=============================================================================

    const double c = Constants.c; // m/s

    detadx[0] = c/Hp_of_x(x);

    return GSL_SUCCESS;
  };

  //=============================================================================
  // TODO: Set the initial condition, set up the ODE system, solve and make
  // the spline eta_of_x_spline 
  //=============================================================================

  Vector eta_ini{0.0};

  ODESolver ode;
  ode.solve(detadx, x_array, eta_ini);

  auto eta_array = ode.get_data_by_component(0);

  eta_of_x_spline.create(x_array, eta_array, "eta");

  Utils::EndTiming("Eta");

  Utils::StartTiming("t");

  // The ODE for deta/dx
  ODEFunction dtdx = [&](double x, const double *t, double *dtdx){

    //=============================================================================
    // TODO: Set the rhs of the dtdx ODE
    //=============================================================================

    dtdx[0] = 1/H_of_x(x);

    return GSL_SUCCESS;
  };

  //=============================================================================
  // TODO: Set the initial condition, set up the ODE system, solve and make
  // the spline t_of_x_spline 
  //=============================================================================

  Vector t_ini{1./(2.*H_of_x(0.0))};
  
  ode.solve(dtdx, x_array, t_ini);

  auto t_array = ode.get_data_by_component(0);

  t_of_x_spline.create(x_array, t_array, "t");

  Utils::EndTiming("t");
  
}

//====================================================
// Get methods
//====================================================

double BackgroundCosmology::H_of_x(double x) const{

  double a = exp(x);

  return H0*pow((OmegaB+OmegaCDM)*pow(a, -3)+(OmegaR+OmegaNu)*pow(a, -4)+OmegaK*pow(a, -2)+OmegaLambda, 1./2.);
}

double BackgroundCosmology::Hp_of_x(double x) const{

  return exp(x)*H_of_x(x);
}

double BackgroundCosmology::dHpdx_of_x(double x) const{

  //=============================================================================
  // TODO: Implement...
  //=============================================================================
  
  double a = exp(x);

  double root = (OmegaB+OmegaCDM)*pow(a, -3)+(OmegaR+OmegaNu)*pow(a, -4)+OmegaK*pow(a, -2)+OmegaLambda;

  double drootdx = -3*(OmegaB+OmegaCDM)*pow(a, -3)-4*(OmegaR+OmegaNu)*pow(a, -4)-2*OmegaK*pow(a, -2);

  double ddrootddx = 9*(OmegaB+OmegaCDM)*pow(a, -3)+16*(OmegaR+OmegaNu)*pow(a, -4)+4*OmegaK*pow(a, -2);

  return Hp_of_x(x)+a*H0*drootdx/(2*pow(root, 1./2.));
}

double BackgroundCosmology::ddHpddx_of_x(double x) const{

  //=============================================================================
  // TODO: Implement...
  //=============================================================================

  double a = exp(x);

  double in_root = (OmegaB+OmegaCDM)*pow(a, -3)+(OmegaR+OmegaNu)*pow(a, -4)+OmegaK*pow(a, -2)+OmegaLambda;

  double din_rootdx = -3*(OmegaB+OmegaCDM)*pow(a, -3)-4*(OmegaR+OmegaNu)*pow(a, -4)-2*OmegaK*pow(a, -2);

  double ddin_rootddx = 9*(OmegaB+OmegaCDM)*pow(a, -3)+16*(OmegaR+OmegaNu)*pow(a, -4)+4*OmegaK*pow(a, -2);

  return 2*dHpdx_of_x(x)-Hp_of_x(x)+(a*H0/2)*(ddin_rootddx/pow(in_root, 1./2.)-din_rootdx/(2*pow(in_root, 3./2.)));
}

double BackgroundCosmology::get_OmegaB(double x) const{ 

  if(x == 0.0) return OmegaB;

  return OmegaB*pow(H0, 2)/(exp(x)*pow(Hp_of_x(x), 2));
}

double BackgroundCosmology::get_OmegaR(double x) const{ 

  if(x == 0.0) return OmegaR;
  
  return OmegaR*pow(H0, 2)/(exp(2*x)*pow(Hp_of_x(x), 2));
}

double BackgroundCosmology::get_OmegaNu(double x) const{ 

  if(x == 0.0) return OmegaNu;
  
  return OmegaNu*pow(H0, 2)/(exp(2*x)*pow(Hp_of_x(x), 2));
}

double BackgroundCosmology::get_OmegaCDM(double x) const{ 

  if(x == 0.0) return OmegaCDM;
  
  return OmegaCDM*pow(H0, 2)/(exp(x)*pow(Hp_of_x(x), 2));
}

double BackgroundCosmology::get_OmegaLambda(double x) const{ 

  if(x == 0.0) return OmegaLambda;
  
  return OmegaLambda*pow(H0, 2)/pow(H_of_x(x), 2);
}

double BackgroundCosmology::get_OmegaK(double x) const{ 

  if(x == 0.0) return OmegaK;
  
  return OmegaK*pow(H0, 2)/pow(Hp_of_x(x), 2);
}
    
double BackgroundCosmology::get_luminosity_distance_of_x(double x) const{
  
  return get_comoving_distance_of_x(x)*exp(-x);
}

double BackgroundCosmology::get_angular_distance_of_x(double x) const{
  
  return get_comoving_distance_of_x(x)*exp(-x);
}

double BackgroundCosmology::get_comoving_distance_of_x(double x) const{

  //=============================================================================
  // TODO: Implement...
  //=============================================================================  

  return eta_of_x(0)-eta_of_x(x);
}

double BackgroundCosmology::eta_of_x(double x) const{

  return eta_of_x_spline(x);
}

double BackgroundCosmology::t_of_x(double x) const{

  return t_of_x_spline(x);
}

double BackgroundCosmology::get_H0() const{ 

  return H0; 
}

double BackgroundCosmology::get_h() const{ 

  return h; 
}

double BackgroundCosmology::get_Neff() const{ 

  return Neff; 
}

double BackgroundCosmology::get_TCMB(double x) const{ 

  if(x == 0.0) return TCMB;

  return TCMB * exp(-x); 
}

//====================================================
// Print out info about the class
//====================================================
void BackgroundCosmology::info() const{ 
  std::cout << "\n";
  std::cout << "Info about cosmology class:\n";
  std::cout << "OmegaB:      " << OmegaB      << "\n";
  std::cout << "OmegaCDM:    " << OmegaCDM    << "\n";
  std::cout << "OmegaLambda: " << OmegaLambda << "\n";
  std::cout << "OmegaK:      " << OmegaK      << "\n";
  std::cout << "OmegaNu:     " << OmegaNu     << "\n";
  std::cout << "OmegaR:      " << OmegaR      << "\n";
  std::cout << "Neff:        " << Neff        << "\n";
  std::cout << "h:           " << h           << "\n";
  std::cout << "TCMB:        " << TCMB        << "\n";
  std::cout << std::endl;
} 

//====================================================
// Output some data to file
//====================================================
void BackgroundCosmology::output(const std::string filename) const{
  const double x_min = -10.0;
  const double x_max =  0.0;
  const int    n_pts =  1000;
  
  Vector x_array = Utils::linspace(x_min, x_max, n_pts);

  std::ofstream fp(filename.c_str());
  auto print_data = [&] (const double x) {
    fp << x                               << " ";
    fp << eta_of_x(x)                     << " ";
    fp << Hp_of_x(x)                      << " ";
    fp << dHpdx_of_x(x)                   << " ";
    fp << get_OmegaB(x)                   << " ";
    fp << get_OmegaCDM(x)                 << " ";
    fp << get_OmegaLambda(x)              << " ";
    fp << get_OmegaR(x)                   << " ";
    fp << get_OmegaNu(x)                   << " ";
    fp << get_luminosity_distance_of_x(x) << " ";
    fp << ddHpddx_of_x(x)                 << " ";
    fp << t_of_x(x)                       << " ";
    fp <<"\n";
  };
  std::for_each(x_array.begin(), x_array.end(), print_data);
}

