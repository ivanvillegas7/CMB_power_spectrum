#ifndef _RECOMBINATION_HISTORY_HEADER
#define _RECOMBINATION_HISTORY_HEADER
#include <vector>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include "Utils.h"
#include "BackgroundCosmology.h"

using Vector = std::vector<double>;

class RecombinationHistory{
  private:

    // The cosmology we use
    BackgroundCosmology *cosmo = nullptr;
    
    // Helium fraction
    double Yp;
 
    // The start and end points for recombination arrays (can be modified)
    const double x_start  = Constants.x_start;
    const double x_end    = Constants.x_end;
    
    // Adhoc method for finding the time of 
    double x_rec_end;

    // Numbers of points of Xe,ne array (modify as you see fit)
    const int npts_rec_arrays = 1e5;//4000;
  
    // Xe for when to switch between Saha and Peebles
    const double Xe_saha_limit = 0.99;

    // Store the x-value in transition between regimes
    double x_Saha_to_Peebles = 0;

    // Store the index in transition between regimes
    int idx_Peebles_transition = 0;

    //===============================================================
    // [1] Computation of Xe (Saha and Peebles equation)
    //===============================================================
 
    // Compute Xe from the Saha equation
    std::pair<double,double> electron_fraction_from_saha_equation(double x) const;
    
    // Right hand side of the dXedx Peebles equation
    int rhs_peebles_ode(double x, const double *y, double *dydx);
    
    // Solve for Xe 
    void solve_number_density_electrons();
    
    //===============================================================
    // [2] Compute τ and visibility functions
    //===============================================================

    // The two things we need to solve: Xe/ne and τ
    void solve_for_optical_depth_tau();

    // Solve for s
    void solve_for_sound_horizon();

    // Splines contained in this class
    Spline log_Xe_of_x_spline{"X_e"};
    Spline Xe_of_x_spline_only_Saha{"X_e Saha"};
    Spline log_ne_of_x_spline{"n_e"};
    Spline tau_of_x_spline{"τ"}; 
    Spline dtau_of_x_spline{"dτ"}; 
    Spline g_tilde_of_x_spline{"g̃"};
    Spline dg_tilde_of_x_spline{"dg̃"};
    Spline sound_horizon_of_x_spline{"s"};

  public:

    // Construtors
    RecombinationHistory() = delete;
    RecombinationHistory(
        BackgroundCosmology *cosmo, 
        double Yp);

    // Do all the solving
    void solve();
    
    // Print some useful info about the class
    void info() const;

    // Print sound at decoupling
    void sound_horizon() const;

    // Output some data to file
    void output(const std::string filename) const;

    // Get functions that we must implement
    double tau_of_x(double x) const;
    double dtaudx_of_x(double x) const;
    double ddtauddx_of_x(double x) const;
    double g_tilde_of_x(double x) const;
    double dgdx_tilde_of_x(double x) const;
    double ddgddx_tilde_of_x(double x) const;
    double Xe_of_x(double x) const;
    double Xe_of_x_Saha_approx(double x) const;
    double ne_of_x(double x) const;
    double get_Yp() const;
    double R_of_x(double x) const;
    double cs_of_x(double x) const;
    double sound_horizon_of_x(double x) const;
    double get_number_density_B(double x) const;
    double get_number_density_H(double x) const;
};

#endif
