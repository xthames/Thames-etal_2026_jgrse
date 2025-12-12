[![DOI](https://zenodo.org/badge/265254045.svg)](https://zenodo.org/doi/10.5281/zenodo.10442485)

<!-- Get rid of the metarepo instructions (the two sections below this) once you're done. -->

# metarepo
## [Check out the website for instructions](https://immm-sfa.github.io/metarepo)
`metarepo` is short for meta-repository, a GitHub repository that contains instructions to reproduce results in a published work. This repo is a template for creating your own metarepo.

## Purpose
A meta-repository creates a single point of access for someone to find all of the components that were used to create a published work for the purpose of reproducibility. This repository should contain references to all minted data and software as well as any ancillary code used to transform the source data, create figures for your publication, conduct the experiment, and / or execute the contributing software.

<!-- Get rid of the metarepo instructions (the two sections above this) once you're done. -->

# Thames-etal_2025_jgrse

**Reconciling Coupled Thermal-Water Evolution Models of Earth with Observations through Variable Regassing Efficiency**

Alexander B. Thames<sup>1</sup> and Bradford J. Foley<sup>1\*</sup>

<sup>1 </sup>Department of Geosciences, The Pennsylvania State University, University Park, Pennsylvania 16802, USA.


\* corresponding author:  bjf5382@psu.edu

## Abstract
Reconciling Earth's thermal evolution with geochemical observations remains a fundamental challenge in geophysics. Parameterized convection models typically predict high present-day Urey ratios that are inconsistent with geochemical estimates. Previous studies have considered mantle water concentration as an important factor in modulating Earth's thermal evolution as it could weaken the feedbacks between mantle temperature, viscosity, and heat loss. Dehydration of the mantle when temperatures are high followed by rehydration as the mantle cools can keep mantle viscosity roughly constant in time, allowing primordial heat to be retained such that the present-day Urey ratio is within geochemical constraints. However, previous coupled thermal-water evolution models either still produced Urey ratios that are too high or calculated mantle water concentration *a posteriori* rather than through a model of surface-interior water exchange. We propose that time variability in the efficiency of the return of water to the mantle (regassing) provides a physical mechanism to help resolve these discrepancies. Using coupled thermal-water evolution models with Latin hypercube sampling across uncertain parameters, we incorporate this time-dependent regassing efficiency into two possible cases of water evolution, one where the upper mantle water mass is fixed relative to the total mantle water mass and one where it is variable. Both models successfully produce realizations that match present-day constraints on mantle temperature, mantle water concentration, heat loss, plate speed, and the Urey ratio while aligning with petrological estimates of Archean-Proterozoic mantle temperatures and Phanerozoic freeboard constancy. Variable regassing efficiency allows mantle water concentration to vary inversely with temperature by 1-2 orders of magnitude, keeping viscosity approximately constant such that primordial heat can be retained because the feedback between mantle temperature and heat loss is significantly weakened. The mantle water history we model means that the surface ocean mass is larger than the present day mass during the Archean and Proterozoic, with net regassing over the last ~1 Gyr lowering surface water mass to the present-day value. We find that as a result Earth's total water budget is constrained to 1.25-1.7 ocean masses across more than 95\% of successful realizations, with the maximum surface water mass seen in Earth's past being ~1.4 ocean masses. Low total water budgets are critical to the success of these variable regassing models as they allow the upper mantle to sufficiently dehydrate during the Archean as the mantle warms, while still being able to regas the excess surface water to reach the present-day ocean mass.

## Journal reference
Thames, A. B. & Foley, B. J. Reconciling Coupled Thermal-Water Evolution Models of Earth with Observations through Variable Regassing Efficiency. *Journal of Geophysical Research: Solid Earth*. (In Preparation).

## Data reference
### Output data
Thames, A. B. & Foley, B. J. (2025). Output Data for Thames & Foley -- Reconciling Coupled Thermal-Water Evolution Models (v1.0.0) [Dataset]. Zenodo. https://doi.org/10.5281/zenodo.17903920

## Reproduce my experiment
1. Download all scripts from `workflow` to a common directory 
2. Review scripts `SubmitExploreAndSolve.sh` and `ExploreAndSolveGeneralTXum_2025XI.m`. Add filepaths to this common directory such that each filepath in these scripts points to a valid location. Change usernames, path separators, etc where/if appropriate. The minimum viable directory hierarchy looks like:
     * `ExploreAndSolveGeneralTXum_2015XI.m`
     * `ExploreAndSolveTargetedTXum_2015XI.m`
     * `SubmitExploreAndSolve.sh`
     * `VolatileProcessing_2025XI.m`
     * `models/`
         * `fwd/`
             * `base/`
             * `var/`
         * `rev/`
     * `plots/`
     * `supplemental/`

   The experiment can be run in two formats:
    * through an HPC cluster with `MATLAB2023a` or higher and a workload manager/job scheduler installed (e.g., using SLURM in a Linux environment)
    * on a local machine with `MATLAB2023a` or higher installed

   Both formats produce identical results, but using the HPC cluster performs the experiment approximately an order of magnitude faster

3. **If using the HPC format**, submit `SubmitExploreAndSolve.sh` as a job using the appropriate syntax (*note this example considers a Linux/SLURM environment*):
    | Script Name | Description | How to Run |
    | --- | --- | --- |
    | `SubmitExploreAndSolve.sh` | Shell script to submit experiment as job | `sbatch SubmitExploreAndSolve.sh scriptType timeDirection viscDep viscStrength modelType varyfR varyphiRum` |
    * `scriptType`: can be `general` or `targeted`. Determines the script to run, either `ExploreAndSolveGeneralTXum_2025XI.m` or `ExploreAndSolveTargetedTXum_2025XI.m`
    * `timeDirection`: can be `fwd` or `rev` (*note that while accepted `rev` has legacy functionality and is not used in this experiment*). Determines if the model starts at t=0 Gyr and runs forwards or at t=4.54 Gyr and runs backwards
    * `viscDep`: can be `conc` or `fug`. Determines if viscosity is dependent on mantle water concentration or mantle water fugacity
    * `viscStrength`: can be `strong` or `weak`. Determines the value of the power law for water in viscosity
    * `modelType`: can be `base` or `var`. Determines if the model can accept variable regassing (`var`) or not (`base`)
    * `varyfR`: can be `y` or `n`. Only necessary if `modelType=var`. Activates variable regassing
    * `varyphiRum`: can be `y` or `n`. Only necessary if `modelType=var`. Switches water mass conservation equation from EQ2 in the paper (default) to EQ3    
4. **If using a local machine**, the same environment variables listed in 3. can be directly set in the section\
`%% CONTROLLING WHAT TYPE OF MODEL TO RUN`
    | Script Name | Description | How to Run |
    | --- | --- | --- |
    | `ExploreAndSolveGeneralTXum_2025XI.m` | Script to run experiment with unrestricted parameter space | Execute script in IDE |
    | `ExploreAndSolveTargetedTXum_2025XI.m` | Script to run experiment with targeted parameter space | Execute script in IDE |
5. Both a file of chosen inputs/calculated present-day outputs (with suffix `MC`) and a file of time-series data (with suffix `TS`) are created by running the above scripts and saved to the corresponding folder in the hierarchy. Naming outputs happens automatically based on the environment variables:
    | Environment Variable | Option | Output Dataset Identifier |
    | --- | --- | --- |
    | `scriptType` | `general`, `targeted` | -, `Targeted` |
    | `timeDirection` | `fwd`, `rev` | `FWD`, `REV` |
    | `viscDep` | `conc`, `fug` | `CONC`, `FUG` |
    | `viscStrength` | `strong`, `weak` | `STRONG`, `WEAK` |
    | `modelType` | `base`, `var` | `BASE`, - |
    | `varyfR` | `n`, `y` | -, `VarfR` |
    | `varyphiRum` | `n`, `y` | -, `VarphiRum` |
    
    However, additional distinguishing language is sometimes needed for the outputs. When considering models with variable regassing efficiency, multiple runs may be necessary to find a sufficient number of realizations that align with present-day estimates. If so, you can manually distinguish between each "set" by including `_set#_` after the `_viscStrength_` identifier and before `_VarfR_` in the output file name (see `VolatileProcesing_2025XI.m` for specific examples)

6. The targeted parameter space used in the paper has been left in the corresponding script; determining a new targeted parameter space can be performed in the following way:
    | Script Name | Description | Section | How to Run |
    | --- | --- | --- | --- |
    | `VolatileProcessing_2025XI.m` | Identify parameter ranges discussed in Supporting Information | `%% FIND THE RANGES OF PARAMETERS TO TARGET` | Comment out other sections, execute script in IDE |
    | `VolatileProcessing_2025XI.m` | Identify parameter ranges discussed in main paper | `%% TARGETED OUTPUT ASSESSMENT USING JUST BEST-KNOWN CONSTRAINTS` | Comment out other sections, execute script in IDE |

    Additional distinguishing language is sometimes needed for these outputs, too. When running `ExploreAndSolveTargetedTXum_2025XI.m`, using the reduced set of present-day observations or the full set (see the paper's Supporting Information for more) can be tagged by including either `Essential` and `Full` in the filename (see `VolatileProcesing_2025XI.m` for specific examples) 

## Reproduce my figures
After creating the output data from the steps outlined above -- or downloading the output data used in the paper from the Zenodo repository and placing them in the hierarchy listed above -- to recreate the figures simply proceed as follows:

 1. Download supplementary data from [Herzberg et al. (2010)](https://doi.org/10.1016/j.epsl.2010.01.022). Rename to `Herzberg_etal_2010_Supplemental.xls` and store in `supplemental/`
 2. Download supplementary data from [Condie et al. (2016)](https://doi.org/10.1016/j.gsf.2016.01.006). Rename to `Condie_etal_2016_Supplemental.xlsx` and store in `supplemental/`

| Figure Number(s) | Script Name | Description | Section | How to Run |
| --- | --- | --- | --- | --- |
| 2-3 | `VolatileProcessing_2025XI.m` | Time-series of coupled thermal-water evolution of Earth's mantle using constant regassing efficiency | `%% LOAD IN THE FWD BASE TIMESERIES DATA` | Comment out other sections, execute script in IDE |
| 4 | `VolatileProcessing_2025XI.m` | Investigating shape of variable regassing efficiency from temperature and plate speed | `%% LOAD IN THE FWD BASE TIMESERIES DATA` | Comment out other sections, execute script in IDE |
| 5 | `VolatileProcessing_2025XI.m` | Investigating experimental Magni et al. (2014) parameters | `%% LOAD IN THE FWD BASE TIMESERIES DATA` | Comment out other sections, execute script in IDE |
| 6-9, S1-2, S4-5 | `VolatileProcessing_2025XI.m` | Time-series of coupled thermal-water evolution of Earth's mantle using variable regassing efficiency | `%% LOAD IN ALL-CONSTRAINTS TARGETED VARIANT CASES FOR TIME-SERIES ANALYSIS` | Comment out other sections, execute script in IDE |
| 10 | `VolatileProcessing_2025XI.m` | Time-series of coupled thermal-water evolution of Earth's mantle using variable regassing efficiency | `%% SUCCESS INPUT FACTOR MAPPING` | Comment out other sections, execute script in IDE |
| S3, S6 | `VolatileProcessing_2025XI.m` | Multipanel plots of water-related variables compared to present-day observations | `%% TARGETED OUTPUT ASSESSMENT USING JUST BEST-KNOWN CONSTRAINTS` | Comment out other sections, execute script in IDE |

These recreated figures should mirror the files found in `figures/`.

