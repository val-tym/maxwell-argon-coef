#include <iostream>
#include <map>
#include <filesystem>
#include <fstream>
#include <cmath>

using namespace std;
namespace fs = filesystem;

const int n = 10;

const map<int, pair<int, int> > sims = {
    {1, pair<int, int>{85, 2000}},
    //{5, pair<int, int>{90, 2000}},
    {10, pair<int, int>{95, 2000}},
    //{15, pair<int, int>{105, 1540}},
    {20, pair<int, int>{105, 1150}},
    //{25, pair<int, int>{110, 920}},
    {30, pair<int, int>{115, 765}},
    //{35, pair<int, int>{115, 655}},
    {40, pair<int, int>{120, 570}},
    //{45, pair<int, int>{120, 505}},
    {50, pair<int, int>{120, 460}},
    //{55, pair<int, int>{125, 420}},
    {60, pair<int, int>{125, 380}},
    //{65, pair<int, int>{125, 355}},
    {70, pair<int, int>{130, 330}},
    //{75, pair<int, int>{130, 310}},
    {80, pair<int, int>{130, 290}},
    //{85, pair<int, int>{130, 275}},
    {90, pair<int, int>{130, 265}},
    //{95, pair<int, int>{135, 255}},
    {100, pair<int, int>{135, 240}},
    //{105, pair<int, int>{135, 230}},
    {110, pair<int, int>{135, 225}},
    //{115, pair<int, int>{135, 220}},
    {120, pair<int, int>{135, 210}},
    //{125, pair<int, int>{135, 205}},
    {130, pair<int, int>{140, 200}},
    //{135, pair<int, int>{145, 195}},
    {140, pair<int, int>{140, 190}},
    //{145, pair<int, int>{140, 190}},
    {150, pair<int, int>{140, 185}},
    //{155, pair<int, int>{140, 180}},
    {160, pair<int, int>{140, 180}},
    //{165, pair<int, int>{145, 175}},
    {170, pair<int, int>{145, 170}},
    //{175, pair<int, int>{145, 170}},
    {180, pair<int, int>{145, 170}},
    //{185, pair<int, int>{145, 165}},
    {190, pair<int, int>{145, 160}},
    //{195, pair<int, int>{150, 160}},
    {200, pair<int, int>{150, 160}}
};

int main() {
    for (const auto &sim: sims) {
        int rho = sim.first;
        int t_min = sim.second.first;
        int t_max = sim.second.second;

        if (!fs::exists(to_string(rho))) {
            if (!fs::create_directory(to_string(rho))) {
                cerr << "Failed to create folder." << endl;
                return 1;
            }
        }
        float l = pow(399478 / 6.02 / rho, 1.0f / 3.0f); //(M/Na/rho)^(1/3)
        for (int t = t_min; t <= t_max; t += max((t_max - t_min)/10,5)) {
            if (!fs::exists(to_string(rho) + "/" + to_string(t))) {
                if (!fs::create_directory(to_string(rho) + "/" + to_string(t))) {
                    cerr << "Failed to create folder." << endl;
                    return 1;
                }
            }
            ofstream argon_pdb("./" + to_string(rho) + "/" + to_string(t) + "/argon.pdb");
            if (!argon_pdb) {
                cerr << "Failed to create file." << endl;
                return 1;
            }
            argon_pdb << "CRYST1   " << to_string(l) << "   " << to_string(l) << "   " << to_string(l) <<
                    "  90.00  90.00  90.00 P 1           1\n"
                    "ATOM      1  Ar   Ar     1       0.000   0.000   0.000  1.00  0.00\n"
                    "TER\n"
                    "END";
            argon_pdb.close();

            ofstream argon_top("./" + to_string(rho) + "/" + to_string(t) + "/argon.top");
            if (!argon_top) {
                cerr << "Failed to create file." << endl;
                return 1;
            }
            argon_top << "[ defaults ]\n"
                    "; nbfunc        comb-rule       gen-pairs       fudgeLJ fudgeQQ\n"
                    "  1             1               no              1.0     1.0\n"
                    "\n"
                    "[ atomtypes ]\n"
                    "AR  39.948    0.0   A     0.00622127     9.69576e-06\n"
                    "\n"
                    "[ moleculetype ]\n"
                    "; molname       nrexcl\n"
                    "Ar              1\n"
                    "\n"
                    "[ atoms ]\n"
                    "; id    at type res nr  residu name     at name  cg nr   charge\n"
                    "1       AR       1       Ar              Ar        1       0\n"
                    "\n"
                    "[ system ]\n"
                    "; Name\n"
                    "Argon\n"
                    "\n"
                    "[ molecules ]\n"
                    "; Compound        #mols\n"
                    "AR                " << to_string(pow(n, 3));
            argon_top.close();

            ofstream md_mdp("./" + to_string(rho) + "/" + to_string(t) + "/md.mdp");
            if (!md_mdp) {
                cerr << "Failed to create file." << endl;
                return 1;
            }
            md_mdp << ";\n"
                    ";	File 'mdout.mdp' was generated\n"
                    ";	By user: bert (1001)\n"
                    ";	On host: bertp3\n"
                    ";	At date: Sat Dec  4 13:41:42 2004\n"
                    ";\n"
                    "\n"
                    "; VARIOUS PREPROCESSING OPTIONS\n"
                    "include                  = \n"
                    "define                   = \n"
                    "\n"
                    "; RUN CONTROL PARAMETERS\n"
                    "integrator               = md\n"
                    "; Start time and timestep in ps\n"
                    "tinit                    = 0\n"
                    "dt                       = 0.001\n"
                    "nsteps                   = 1000000\n"
                    "; For exact run continuation or redoing part of a run\n"
                    "init_step                = 0\n"
                    "; mode for center of mass motion removal\n"
                    "comm-mode                = Linear\n"
                    "; number of steps for center of mass motion removal\n"
                    "nstcomm                  = 1\n"
                    "; group(s) for center of mass motion removal\n"
                    "comm-grps                = \n"
                    "\n"
                    "; LANGEVIN DYNAMICS OPTIONS\n"
                    "; Temperature, friction coefficient (amu/ps) and random seed\n"
                    "bd-fric                  = 0\n"
                    "ld-seed                  = 1993\n"
                    "\n"
                    "; ENERGY MINIMIZATION OPTIONS\n"
                    "; Force tolerance and initial step-size\n"
                    "emtol                    = 100\n"
                    "emstep                   = 0.1\n"
                    "; Max number of iterations in relax_shells\n"
                    "niter                    = 20\n"
                    "; Step size (1/ps^2) for minimization of flexible constraints\n"
                    "fcstep                   = 0\n"
                    "; Frequency of steepest descents steps when doing CG\n"
                    "nstcgsteep               = 1000\n"
                    "nbfgscorr                = 10\n"
                    "\n"
                    "; OUTPUT CONTROL OPTIONS\n"
                    "; Output frequency for coords (x), velocities (v) and forces (f)\n"
                    "nstxout                  = 1000\n"
                    "nstvout                  = 1000\n"
                    "nstfout                  = 0\n"
                    "; Output frequency for energies to log file and energy file\n"
                    "nstlog                   = 10000\n"
                    "nstcalcenergy            = 1\n"
                    "nstenergy                = 10000\n"
                    "; Output frequency and precision for xtc file\n"
                    "nstxout-compressed       = 500\n"
                    "compressed-x-precision   = 1000\n"
                    "; This selects the subset of atoms for the xtc file. You can\n"
                    "; select multiple groups. By default all atoms will be written.\n"
                    "xtc-grps                 = \n"
                    "; Selection of energy groups\n"
                    "energygrps               = System\n"
                    "\n"
                    "; long-range cut-off for switched potentials\n"
                    "cutoff-scheme            = Verlet\n"
                    "\n"
                    "; NEIGHBORSEARCHING PARAMETERS\n"
                    "; nblist update frequency\n"
                    "nstlist                  = 50\n"
                    "; Periodic boundary conditions: xyz (default), no (vacuum)\n"
                    "; or full (infinite systems only)\n"
                    "pbc                      = xyz\n"
                    "; nblist cut-off        \n"
                    "rlist                    = 0.85\n"
                    "\n"
                    "\n"
                    "; OPTIONS FOR ELECTROSTATICS AND VDW\n"
                    "; Method for doing electrostatics\n"
                    "coulombtype              = Cut-off\n"
                    "rcoulomb-switch          = 0\n"
                    "rcoulomb                 = 0.85\n"
                    "; Dielectric constant (DC) for cut-off or DC of reaction field\n"
                    "epsilon-r                = 1\n"
                    "; Method for doing Van der Waals\n"
                    "vdw-type                 = Cut-off\n"
                    "; cut-off lengths       \n"
                    "rvdw-switch              = 0\n"
                    "rvdw                     = 0.85\n"
                    "; Apply long range dispersion corrections for Energy and Pressure\n"
                    "DispCorr                 = Enerpres\n"
                    "; Extension of the potential lookup tables beyond the cut-off\n"
                    "table-extension          = 1\n"
                    "; Spacing for the PME/PPPM FFT grid\n"
                    "fourierspacing           = 0.12\n"
                    "; FFT grid size, when a value is 0 fourierspacing will be used\n"
                    "fourier_nx               = 0\n"
                    "fourier_ny               = 0\n"
                    "fourier_nz               = 0\n"
                    "; EWALD/PME/PPPM parameters\n"
                    "pme_order                = 4\n"
                    "ewald_rtol               = 1e-05\n"
                    "ewald_geometry           = 3d\n"
                    "epsilon_surface          = 0\n"
                    "\n"
                    "\n"
                    "; IMPLICIT SOLVENT (for use with Generalized Born electrostatics)\n"
                    "implicit_solvent         = No\n"
                    "\n"
                    "; OPTIONS FOR WEAK COUPLING ALGORITHMS\n"
                    "; Temperature coupling  \n"
                    "Tcoupl                   = no\n"
                    "; Groups to couple separately\n"
                    "tc-grps                  = System\n"
                    "; Time constant (ps) and reference temperature (K)\n"
                    "tau_t                    = 0.1\n"
                    "ref_t                    = " << to_string(t) << "\n"
                    "\n"
                    "; Simulated Annealing (Linear)\n"
                    "annealing                = no\n"
                    "annealing_npoints        = 2\n"
                    "annealing_time           = 0 1000       ; Total simulation time (ps)\n"
                    "annealing_temp           = " << t << " " << t << "\n"
                    "\n"
                    "; Pressure coupling     \n"
                    "Pcoupl                   = no\n"
                    "Pcoupltype               = isotropic\n"
                    "; Time constant (ps), compressibility (1/bar) and reference P (bar)\n"
                    "tau_p                    = 0.5\n"
                    "compressibility          = 4.5e-5\n"
                    "ref_p                    = 1.0\n"
                    "\n"
                    "; SIMULATED ANNEALING  \n"
                    "; Type of annealing for each temperature group (no/single/periodic)\n"
                    ";annealing                =\n"
                    "; Number of time points to use for specifying annealing in each group\n"
                    ";annealing_npoints        =\n"
                    "; List of times at the annealing points for each group\n"
                    ";annealing_time           =\n"
                    "; Temp. at each annealing point, for each group.\n"
                    ";annealing_temp           =\n"
                    "\n"
                    "; GENERATE VELOCITIES FOR STARTUP RUN\n"
                    "gen_vel                  = yes\n"
                    "gen_temp                 = " << to_string(t) << "\n"
                    "gen_seed                 = 173529\n"
                    "\n"
                    "; OPTIONS FOR BONDS    \n"
                    "constraints              = all-bonds\n"
                    "; Type of constraint algorithm\n"
                    "constraint-algorithm     = Lincs\n"
                    "continuation      = no\n"
                    "; Use successive overrelaxation to reduce the number of shake iterations\n"
                    "Shake-SOR                = no\n"
                    "; Relative tolerance of shake\n"
                    "shake-tol                = 1e-04\n"
                    "; Highest order in the expansion of the constraint coupling matrix\n"
                    "lincs-order              = 4\n"
                    "; Number of iterations in the final step of LINCS. 1 is fine for\n"
                    "; normal simulations, but use 2 to conserve energy in NVE runs.\n"
                    "; For energy minimization with constraints it should be 4 to 8.\n"
                    "lincs-iter               = 1\n"
                    "; Lincs will write a warning to the stderr if in one step a bond\n"
                    "; rotates over more degrees than\n"
                    "lincs-warnangle          = 30\n"
                    "; Convert harmonic bonds to morse potentials\n"
                    "morse                    = no\n"
                    "\n"
                    "; ENERGY GROUP EXCLUSIONS\n"
                    "; Pairs of energy groups for which all non-bonded interactions are excluded\n"
                    "energygrp_excl           = \n"
                    "\n"
                    "; NMR refinement stuff \n"
                    "; Distance restraints type: No, Simple or Ensemble\n"
                    "disre                    = No\n"
                    "; Force weighting of pairs in one distance restraint: Conservative or Equal\n"
                    "disre-weighting          = Conservative\n"
                    "; Use sqrt of the time averaged times the instantaneous violation\n"
                    "disre-mixed              = no\n"
                    "disre-fc                 = 1000\n"
                    "disre-tau                = 0\n"
                    "; Output frequency for pair distances to energy file\n"
                    "nstdisreout              = 100\n"
                    "; Orientation restraints: No or Yes\n"
                    "orire                    = no\n"
                    "; Orientation restraints force constant and tau for time averaging\n"
                    "orire-fc                 = 0\n"
                    "orire-tau                = 0\n"
                    "orire-fitgrp             = \n"
                    "; Output frequency for trace(SD) to energy file\n"
                    "nstorireout              = 100\n"
                    "\n"
                    "; Free energy control stuff\n"
                    "free-energy              = no\n"
                    "init-lambda              = 0\n"
                    "delta-lambda             = 0\n"
                    "sc-alpha                 = 0\n"
                    "sc-sigma                 = 0.3\n"
                    "\n"
                    "; Non-equilibrium MD stuff\n"
                    "acc-grps                 = \n"
                    "accelerate               = \n"
                    "freezegrps               = \n"
                    "freezedim                = \n"
                    "cos-acceleration         = 0\n"
                    "\n"
                    "; Electric fields      \n"
                    "; Format is number of terms (int) and for all terms an amplitude (real)\n"
                    "; and a phase angle (real)\n"
                    "E-x                      = \n"
                    "E-xt                     = \n"
                    "E-y                      = \n"
                    "E-yt                     = \n"
                    "E-z                      = \n"
                    "E-zt                     = \n"
                    "\n"
                    "; User defined thingies\n"
                    "user1-grps               = \n"
                    "user2-grps               = \n"
                    "userint1                 = 0\n"
                    "userint2                 = 0\n"
                    "userint3                 = 0\n"
                    "userint4                 = 0\n"
                    "userreal1                = 0\n"
                    "userreal2                = 0\n"
                    "userreal3                = 0\n"
                    "userreal4                = 0";
            md_mdp.close();

            system(string("cd " + to_string(rho) + "/" + to_string(t) +
                          " && gmx genconf -f argon.pdb -o argon_start.pdb -nbox " + to_string(n) + " " + to_string(n) +
                          " "
                          + to_string(n) +
                          " && gmx grompp -f md.mdp -c argon_start.pdb -p argon.top" +
                          " && gmx mdrun -s topol.tpr -v -c argon_1ns.gro -nice 0" +
                          " && gmx energy -o temperature.xvg <<< 7" +
                          " && gmx energy -o kinetic.xvg <<< 5").c_str());
        }
    }

    return 0;
}
