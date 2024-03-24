#include "Utils.h"
#include "BackgroundCosmology.h"
#include "RecombinationHistory.h"
#include "Perturbations.h"
#include "PowerSpectrum.h"
#include "SupernovaFitting.h"

int main(int argc, char **argv){
  Utils::StartTiming("CMB power spectrum");
  Utils::StartTiming("Milestone I");

  //=========================================================================
  // Parameters
  //=========================================================================

  // Background parameters
  double h           = 0.67;
  double OmegaB      = 0.05;
  double OmegaCDM    = 0.267;
  double OmegaK      = 0.0;
  double Neff        = 3.046;
  double TCMB        = 2.7255;

  // Recombination parameters
  double Yp          = 0.0; //0.245;

  // Power-spectrum parameters
  double A_s         = 2.1e-9;
  double n_s         = 0.965;
  double kpivot_mpc  = 0.05;

  //=========================================================================
  // Module I
  //=========================================================================

  std::cout<<"\n";

  // Do the supernova fits.
  // Make sure you read the comments on the top of src/SupernovaFitting.h
  //mcmc_fit_to_supernova_data("data/supernovadata.txt", "Results/results_supernovafitting.txt");

  // Set up and solve the background
  BackgroundCosmology cosmo(h, OmegaB, OmegaCDM, OmegaK, Neff, TCMB);
  cosmo.solve();
  cosmo.info();
  
  // Output background evolution quantities
  cosmo.output("Results/cosmology.txt");

  Utils::EndTiming("Milestone I");
  std::cout<<"\n";

  //=========================================================================
  // Module II
  //=========================================================================
  
  Utils::StartTiming("Milestone II");

  // Solve the recombination history
  RecombinationHistory rec(&cosmo, Yp);
  rec.solve();
  rec.info();

  // Output recombination quantities and times
  rec.output("Results/recombination.txt");
  rec.save_time_results("Results/recombination_times.txt");
  // Print recombination times
  rec.print_time_results();
  
  Utils::EndTiming("Milestone II");
  std::cout<<"\n";

  // Remove when Module III is completed
  return 0;

  //=========================================================================
  // Module III
  //=========================================================================
  
  Utils::StartTiming("Milestone III");

  // Solve the perturbations
  Perturbations pert(&cosmo, &rec);
  pert.solve();
  pert.info();
  
  // Output perturbation quantities
  double kvalue = 0.01 / Constants.Mpc;
  pert.output(kvalue, "Results/perturbations_k0.01.txt");

  Utils::EndTiming("Milestone III");
  std::cout<<"\n";
  
  // Remove when Module IV is completed
  return 0;
  
  //=========================================================================
  // Module IV
  //=========================================================================

  Utils::StartTiming("Milestone IV");

  PowerSpectrum power(&cosmo, &rec, &pert, A_s, n_s, kpivot_mpc);
  power.solve();
  power.output("Results/cells.txt");

  Utils::EndTiming("Milestone IV");
  std::cout<<"\n";
  Utils::EndTiming("CMB power spectrum");
  
  // Remove when module is completed
  return 0;
}
