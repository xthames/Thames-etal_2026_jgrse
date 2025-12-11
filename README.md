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
Reference for each minted data source for your output data.  For example:

Human, I.M. (2021). My output dataset name [Data set]. DataHub. https://doi.org/some-doi-number

_your output data references here_

## Reproduce my experiment
1. Download all scripts from `workflow` to a common directory 
2. Review scripts `SubmitExploreAndSolve.sh` and `ExploreAndSolveGeneralTXum_2025XI.m`. Create/restructure this directory such that each existing filepath points to a valid location. The essential directory hierarchy looks like:
 * `models/`
     * `fwd/`
         * `base/`
         * `var/`
     * `rev/`
 * `supplemental/`
The experiment can be run in two formats:
 * through an HPC cluster with MATLAB and a workload manager/job scheduler installed (e.g., using SLURM in a Linux environment)
 * on a local machine with MATLAB installed
Because this experiment can be run in two formats, the appropriate directory hierarchy will depend on the format the experiment is run in. Both formats produce identical results, but using an HPC cluster is appreciably faster

### Exploring models with unrestricted parameter space
1. Install the software components required to conduct the experiment from [contributing modeling software](#contributing-modeling-software)
2. Download and install the supporting [input data](#input-data) required to conduct the experiment
3. Run the following scripts in the `workflow` directory to recreate this experiment:

| Script Name | Description | How to Run |
| --- | --- | --- |
| `step_one.py` | Script to run the first part of my experiment | `python3 step_one.py -f /path/to/inputdata/file_one.csv` |
| `step_two.py` | Script to run the second part of my experiment | `python3 step_two.py -o /path/to/my/outputdir` |

4. Download and unzip the [output data](#output-data) from my experiment 
5. Run the following scripts in the `workflow` directory to compare my outputs to those from the publication

| Script Name | Description | How to Run |
| --- | --- | --- |
| `compare.py` | Script to compare my outputs to the original | `python3 compare.py --orig /path/to/original/data.csv --new /path/to/new/data.csv` |

## Reproduce my figures
Use the scripts found in the `figures` directory to reproduce the figures used in this publication.

| Figure Number(s) | Script Name | Description | How to Run |
| --- | --- | --- | --- |
| 1, 2 | `generate_plot.py` | Description of figure, ie. "Plots the difference between our two scenarios" | `python3 generate_plot.py -input /path/to/inputs -output /path/to/outuptdir` |
| 3 | `generate_figure.py` | Description of figure, ie. "Shows how the mean and peak differences are calculated" | `python3 generate_figure.py -input /path/to/inputs -output /path/to/outuptdir` |

