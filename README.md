# Example Integral REactor for Nuclear Education (EIRENE) 

EIRENE is an non-proprietary reactor model somewhat resembling the IMSR-400 design of Terrestrial Energy Inc. https://aris.iaea.org/PDF/IMSR400.pdf

*Eirene* is a Greek goddess of Peace and the season of spring, a keeper of the gates of heaven along with her sisters Eunomia (Good Order) and Dike (Justice). (Hesiod, The Works and Days, verse 225)

# Description of subdirectories:
* EIRENE: Contains source model files for the EIRENE SCALE 6.3.1 model, including KENO-VI inputs for running neutronics analyses (critical uranium search and fuel and isothermal temperature feedback coefficients); TRITON inputs for depletion analyses (both with and without Sourdough refueling); SHIFT inputs for generating .h5 files from neutronic analyses; and nuclide vector data collected at each depletion step with the Sourdough refueling simulations using different refuel salt enrichment levels.
* Gd-SCRods: Contains source model files for the original EIRENE and SIRENE models in SCALE 6.3.1 and Serpent-2, respectively, where the SS316 shutdown rods were gadolinium burnable absorbers instead. Files include FTC and ITC model inputs for EIRENE and SIRENE, in addition to early EIRENE depletion inputs in SCALE 6.3.1 TRITON, SHIFT files, and the Python scripts used to create and run them.
* OpenEIRENE: Contains source model files for the OpenEIRENE OpenMC model, including the original OpenEIRENE as-is model and all model inputs utilized in calculating fuel and isothermal temperature feedback coefficients.
* SIRENE: Contains source model files for the SIRENE Serpent-2 model, including model input files for neutronic analyses (critical uranium search and fuel and isothermal temperature feedback coefficient); depletion analyses with and without insoluble fission product removal; and the Python scripts used to create and run these model files.
* ThEIRENE: Contains source model files for the Th-EIRENE SCALE 6.3.1 model, a structurally and gemoetrically identical variation of the EIRENE model which uses a mixture of fertile thorium and fissile HALEU fuel. The directory includes model input files for performing neutronics analyses (critical uranium search and fuel and isothermal temperature feedback coefficients); TRITON inputs for depletion analyses (both with and without Sourdough refueling); and nuclide vector data collected at each depletion step with the Sourdough refueling simulations using 19.75% enriched refuel salt.
* util: Contains Python scripts used in the creation of the EIRENE, SIRENE, and OpenEIRENE models and shell scripts for executing the model input files on a computational cluster.

# Published papers:



