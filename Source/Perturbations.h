#ifndef _PERTURBATIONS_HEADER
#define _PERTURBATIONS_HEADER
#ifdef _USEOPENMP
#include <omp.h>
#endif
#include <vector>
#include <fstream>
#include <algorithm>
#include "Utils.h"
#include "BackgroundCosmology.h"
#include "RecombinationHistory.h"

using Vector   = std::vector<double>;
using Vector2D = std::vector<Vector>;

class Perturbations{
  private:

    BackgroundCosmology *cosmo = nullptr;
    RecombinationHistory *rec  = nullptr;
   
    // The scales we integrate over
    const int n_k        = 100;
    const double k_min   = Constants.k_min;
    const double k_max   = Constants.k_max;
    
    // Start and end of the time-integration
    const int n_x        = 1000;
    const double x_start = Constants.x_start;
    const double x_end   = Constants.x_end;

    // Below is a full list of splines you probably need, 
    // but you only need to make the splines you will need

    // Splines of scalar perturbations quantities
    Spline2D delta_CDM_spline{"δ_CDM spline"};
    Spline2D delta_B_spline{"δ_B spline"};
    Spline2D v_CDM_spline{"v_CDM spline"};
    Spline2D v_B_spline{"v_B spline"};
    Spline2D Phi_spline{"Φ spline"};
    Spline2D Pi_spline{"Π spline"};
    Spline2D Psi_spline{"Ψ spline"};
   
    // Splines of source functions (ST for temperature; SE for polarization)
    Spline2D ST_spline{"ST spline"};
    Spline2D SE_spline{"SE spline"};
    
    // Splines of mulipole quantities
    // NB: If you use there you have to allocate the container first
    // e.g. Theta_spline = std::vector<Spline2D>(n_ell_theta); before using it
    std::vector<Spline2D> Theta_spline;
    std::vector<Spline2D> Theta_p_spline;
    std::vector<Spline2D> Nu_spline;
    
    //==========================================================
    // [1] Tight coupling ODE system
    //==========================================================

    // Set the initial conditions at the start (which is in tight coupling)
    Vector set_ic(
        const double x, 
        const double k) const;
    
    // Right hand side of the ODE in the tight coupling regime
    int rhs_tight_coupling_ode(double x, double k, const double *y, double *dydx);
    
    // Compute the time when tight coupling ends
    double get_tight_coupling_time(const double k) const;
    
    //==========================================================
    // [2] The full ODE system 
    //==========================================================
    
    // Set initial condition after tight coupling
    Vector set_ic_after_tight_coupling(
        const Vector &y_tight_coupling, 
        const double x, 
        const double k) const;

    // Right hand side of the ODE in the full regime
    int rhs_full_ode(double x, double k, const double *y, double *dydx);
    
    //==========================================================
    // [3] Integrate the full system
    //==========================================================
    
    // Integrate perturbations and spline the result
    void integrate_perturbations();
    
    //==========================================================
    // [4] Compute source functions from the result
    //==========================================================
    
    // Compute source functions and spline the result
    void compute_source_functions();

  public:

    // Constructors
    Perturbations() = default;
    Perturbations(
        BackgroundCosmology *cosmo, 
        RecombinationHistory *rec); 

    // Do all the solving
    void solve();
    
    // Print some useful info about the class
    void info() const;

    // Output info to file
    void output(const double k, const std::string filename) const;

    // Get the quantities we have integrated
    double get_delta_CDM(const double x, const double k) const;
    double get_delta_B(const double x, const double k) const;
    double get_v_CDM(const double x, const double k) const;
    double get_v_B(const double x, const double k) const;
    double get_dv_Bdx(const double x, const double k) const;
    double get_Phi(const double x, const double k) const;
    double get_dPhidx(const double x, const double k) const;
    double get_Psi(const double x, const double k) const;
    double get_dPsidx(const double x, const double k) const;
    double get_Pi(const double x, const double k) const;
    double get_dPidx(const double x, const double k) const;
    double get_Theta(const double x, const double k, const int ell) const;
    double get_dThetadx(const double x, const double k, const int ell) const;
    double get_Theta_p(const double x, const double k, const int ell) const;
    double get_dTheta_pdx(const double x, const double k, const int ell) const;
    double get_Nu(const double x, const double k, const int ell) const;
    double get_Source_T(const double x, const double k) const;
    double get_Source_E(const double x, const double k) const;
};

#endif
