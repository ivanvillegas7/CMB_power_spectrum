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
  double Neff        = 0.0;
  if (Constants.neutrinos){
    Neff             = 3.046;
  }
  double TCMB        = 2.7255;

  // Recombination parameters
  double Yp          = 0.0;
  if (Constants.Helium)
  {
    Yp               = 0.245;
  }

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
  //mcmc_fit_to_supernova_data("Data/supernovadata.txt", "Results/results_supernovafitting.txt");

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
  
  //Print sound horizon at decoupling
  rec.sound_horizon();
  
  Utils::EndTiming("Milestone II");
  std::cout<<"\n";

  //=========================================================================
  // Module III
  //=========================================================================
  
  Utils::StartTiming("Milestone III");

  // Solve the perturbations
  Perturbations pert(&cosmo, &rec);
  pert.solve();
  pert.info();
  
  // Output perturbation quantities
  double kvalue1 = 0.1 / Constants.Mpc;
  double kvalue01 = 0.01 / Constants.Mpc;
  double kvalue001 = 0.001 / Constants.Mpc;
  double kvalue0001 = 0.0001 / Constants.Mpc;
  pert.output(kvalue1, "Results/perturbations_k1.txt");
  pert.output(kvalue01, "Results/perturbations_k01.txt");
  pert.output(kvalue001, "Results/perturbations_k001.txt");
  pert.output(kvalue0001, "Results/perturbations_k0001.txt");

  Utils::EndTiming("Milestone III");
  std::cout<<"\n";
  
  //=========================================================================
  // Module IV
  //=========================================================================

  Utils::StartTiming("Milestone IV");

  PowerSpectrum power(&cosmo, &rec, &pert, A_s, n_s, kpivot_mpc);
  power.solve();
  power.output("Results/cells.txt");
  power.output_MPS("Results/Matter_PS.txt");

  Utils::EndTiming("Milestone IV");
  std::cout<<"\n";
  Utils::EndTiming("CMB power spectrum");
  std::cout << std::endl;
  
  return 0;
}
