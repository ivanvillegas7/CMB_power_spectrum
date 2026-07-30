#include <cstdlib>
#include <string>
#include "Utils.h"
#include "BackgroundCosmology.h"
#include "RecombinationHistory.h"
#include "Perturbations.h"
#include "PowerSpectrum.h"
#include "SupernovaFitting.h"

// Run a Python plotting script right after the C++ module that generated its
// input has finished. The scripts use relative paths like "../Results/...",
// so they must be run with their working directory set to Python/.
void run_python(const std::string &script){
  // Make sure everything C++ has printed so far is actually on screen
  // before Python starts printing its own output, so the terminal always
  // shows things in the right order.
  std::cout.flush();
  std::string command = "cd Python && python3 " + script;
  int status = std::system(command.c_str());
  if(status != 0){
    std::cerr << "Warning: " << script << " exited with an error (code "
              << status << ")\n";
  }
}

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
  if(Constants.compute_supernova_MCMC){
    mcmc_fit_to_supernova_data("Data/supernovadata.txt", "Results/results_supernovafitting.txt");
    std::cout<<"\n";
  }

  // Set up and solve the background
  BackgroundCosmology cosmo(h, OmegaB, OmegaCDM, OmegaK, Neff, TCMB);
  cosmo.solve();
  cosmo.info();
  
  // Output background evolution quantities
  cosmo.output("Results/cosmology.txt");

  // Generate the Milestone I plots/tables straight away
  run_python("Milestone1.py");

  std::cout<<"\n";
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

  // Write flags + decoupling index to file (read automatically by Python,
  // no more copy-pasting "Index at decoupling" from the terminal)
  rec.output_info("Results/run_info.txt");

  // Generate the Milestone II plots/tables straight away
  run_python("Milestone2.py");

  std::cout<<"\n";
  Utils::EndTiming("Milestone II");
  std::cout<<"\n";

  //=========================================================================
  // Module III
  //=========================================================================
  
  Utils::StartTiming("Milestone III");

  std::cout << "Solving the perturbations for every k (this is the slowest step; progress below):\n\n";

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

  // Generate the Milestone III plots straight away
  run_python("Milestone3.py");

  std::cout<<"\n";
  Utils::EndTiming("Milestone III");
  std::cout<<"\n";
  
  //=========================================================================
  // Module IV
  //=========================================================================

  Utils::StartTiming("Milestone IV");

  PowerSpectrum power(&cosmo, &rec, &pert, A_s, n_s, kpivot_mpc);
  power.solve();
  power.output("Results/cells.txt");
  if(Constants.compute_PS_components){
    power.output_components("Results/cells_SW.txt", "Results/cells_ISW.txt",
                             "Results/cells_DOPPLER.txt", "Results/cells_POL.txt");
  }
  power.output_MPS("Results/Matter_PS.txt");

  // Generate the Milestone IV plots (CMB PS, matter PS, CMB map) straight away
  run_python("Milestone4.py");

  std::cout<<"\n";
  Utils::EndTiming("Milestone IV");
  std::cout<<"\n";

  std::cout<<"\n";
  Utils::EndTiming("CMB power spectrum");
  std::cout << std::endl;
  
  return 0;
}
