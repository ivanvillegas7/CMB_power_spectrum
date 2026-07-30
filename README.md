# Computing the CMB power spectrum

## The Cosmic Microwave Background and the Large Scale Structure of our Universe

This repository will store the code used to compute the CMB power spectrum, a project for the Master's degree course AST5220 - Cosmology II of the University of Oslo (UiO). The structure of the repository and the code is based on Hans A. Whinter's (@HAWinther), following his [webpage](https://cmb.wintherscoming.no/index.php). This repository contains all the needed files for runnig the project, solving the cosmological background, making a MCMC fit to the given supernova data, solving the recombination history and solving the evolution of the Universe, as well as computing the CMB power-spectrum, the matter power-spectrum and generating a CMB map.

The code now runs the full ("PhD") version: it includes Helium in the recombination history, reionization, massless neutrinos and photon polarization, in addition to everything already needed for the Master's degree deliverable. All four can still be switched off independently by setting the corresponding flag (`neutrinos`, `Helium`, `reionization`, `polarization`) to `false` in `Source/Utils.h`, in case you want to reproduce the simpler Master's-only results.

## Compiling and running

Edit the Makefile adding the right paths and compile the code running `make`. This builds the fast (release) `cmb` executable, with OpenMP enabled to parallelize over `k`. If you want a build with STL bounds-checking enabled instead (useful while debugging, at the cost of speed), run `make debug`, which produces `cmb_debug`.

The code runs from Main.cpp and then proceeds to go through the different milestones one by one untill we have the CMB power spectra in the end. Once compiled, run it as `./cmb`. This single command now does everything: it solves the background, the recombination history, the perturbations and the power-spectra, and automatically runs the corresponding Python plotting script after each milestone finishes (no need to run any Python script by hand, and no values printed to the terminal need to be copied anywhere — they are all written to and read from `Results/run_info.txt`). If you want to compile this on your computer you need to install the [GSL](ftp://ftp.gnu.org/gnu/gsl/) library first.

By default, plots are only saved to `Plots/` (no windows pop up), so that `./cmb` runs fully unattended. If you also want to see each plot in a floating window, run `CMB_SHOW_PLOTS=1 ./cmb` instead.

## How to install GSL

See [this](https://solarianprogrammer.com/) for how to install it on a Windows machine. On Linux or a Mac you can either use a package manager or install it directly as follows:

- Go the the home directory: `cd $HOME`

- Make a local folder to hold libraries: `mkdir local`

- Enter this directory: `cd local`

- Download the code (if you don't have wget you need to get the file to this dir by other means): `wget ftp://ftp.gnu.org/gnu/gsl/gsl-2.6.tar.gz` 

- Untar the code: `tar -xvf gsl-2.6.tar.gz`

- You should now have the gsl-2.6 folder. Enter it: `cd gsl-2.6`

- Run the configure script: `./configure --prefix=$HOME/local`

- Compile and install it: `make ; make install`

- In the CMB code Makefile change the include and lib paths to point to the library:

   `-INC  = -I$(HOME)/local/include`
  
   `-LIBS = -L$(HOME)/local/lib -lgsl -lgslcblas`

- If this fails with `libgsl.so not found` then run the command:
  `export LD\_LIBRARY\_PATH="$LD\_LIBRARY\_PATH:$HOME/local/lib"`

and try to run `./cmb` again and it should work. To avoid having to run this command every time you open a new terminal open the `$HOME/.bashrc` file and add this line to the end of the file and it will load everytime you open a new window.

## Plot the results

Running `./cmb` already takes care of generating all the plots: right after each milestone finishes in C++, the corresponding Python script (`Milestone1.py`, `Milestone2.py`, `Milestone3.py`, `Milestone4.py`) is run automatically, so you normally don't need to touch Python at all. If you want to regenerate every plot from Python alone (e.g. after `Results/` has already been populated by a previous `./cmb` run), you can still run `Main.py`. If you want to plot the results, you will need to install '[healpy](https://healpy.readthedocs.io/en/latest/install.html)'. You can do it simply running `conda install healpy` in your favourite conda distribution.

You can also only make the plots for each milestone independently by running the corresponding Python script.

## Telegram notification when the simulation finishes

By default, `./cmb` sends a Telegram message (with the CMB power-spectrum plot attached) as soon as the whole run finishes, using `Python/notification.py`. This is controlled by the `send_notification` flag in `Source/Utils.h` (set it to `false` to disable it).

To receive these notifications on your own Telegram account, you need a bot token and your chat ID:

1. **Create a bot and get its token.** In Telegram, open a chat with [@BotFather](https://t.me/BotFather) and send `/newbot`. Follow the prompts (choose a name and a username for your bot). BotFather will reply with a token that looks like `123456789:AAH...` — this is your `bot_token`.
2. **Start a chat with your new bot.** Search for the bot by the username you just gave it and send it any message (e.g. `/start`), so Telegram has a conversation to deliver messages to.
3. **Get your chat ID.** With the browser (or `curl`), open:
   ```
   https://api.telegram.org/bot<your_bot_token>/getUpdates
   ```
   replacing `<your_bot_token>` with the token from step 1. After sending your bot a message, this page will show a JSON response containing `"chat":{"id": ...}` — that number is your `bot_chatID`.
4. **Edit `Python/notification.py`** and replace the placeholder `bot_token` and `bot_chatID` values with your own.

## Bool flags controlling the pipeline

`Source/Utils.h` has several flags, all following the same pattern (edit, then `make clean && make`):

- `neutrinos`, `Helium`, `reionization`, `polarization`: which physics to include (see above).
- `compute_PS_components`: whether to also compute and output the individual Sachs-Wolfe / ISW / Doppler / polarization contributions to `C_ell^TT` (`cells_SW/ISW/DOPPLER/POL.txt`). This costs 4 extra line-of-sight integrations, the single most expensive step in the pipeline.
- `compute_supernova_MCMC`: whether to (re-)run the MCMC fit to the supernova data. Off by default; a sample results file already ships in `Results/`.
- `send_notification`: whether to send the Telegram notification described above.

## Further documentation

This is a project for a Master's degree course, and each deliverable was submitted as a separate report. The folder named `Reports` stores the report of each milestone: Milestone I (`I: Solving the cosmological background.pdf`), Milestone II (`II: Recombination history of the Universe.pdf`), Milestone III (`III: Evolution of structure in the Universe.pdf`) and Milestone IV (`IV: The CMB and matter power-spectra.pdf`).

There is also a final report, `Computing the CMB power spectrum.pdf`, joining the four milestones together into a single article, using the definitive version of the code (the full, "PhD"-level version described above).
