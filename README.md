# SLiMsims_StabilizingSelection
Code to generate simulations of stabilizing selection on complex traits in SLiM and process output data. 

This repository contains code associated with [Yair and Coop (2022)](https://doi.org/10.1101/2021.09.10.459833 ). We provide an overview of the workflow and contents of each directory below. See those directories and comments in scripts for more details. 

# Simulation framework
Each simulation has a burn-in of 60,000 generations (which 6N generations for N=1e4, the population size we simulated). Since the burn-in is a time-consuming process, we split the burn-in into 2 parts: (1) a burn-in "shared" between groups of 10 simulation replicates for the first 40,000 generations, then (2) an independent 20,000 generation burn-in for each simulation replicate ("unshared"). 

After the burn-in, we track the divergence in the genetic basis of trait variation between two populations that independently diverge from the common ancestor (burn-in state). Since we do not consider migration between populations, we repeat the divergence from the common ancestor twice to represent each population. 

We need to transfer information about the simulation state between each of the three above steps (shared burn-in, unshared burn-in, population divergence). SLiM provides a function to record the state of the simulation, which includes information about each segregating mutation and the individuals that carry it. To reduce the memory required, we let the simulation ignore mutations that become fixed in the simulation, then referred to as substitutions. However, we need information about the effects of these substitutions because they still contribute to an individual's phenotype, which we need when evaluating an individual's fitness. Therefore, we separately record the additive contribution of all substitutions to an individual's phenotype to a text file, and read that information when beginning the next simulation step.

During the divergence simulation, we record data on each mutation to a text file. We merge this data from each population and process it in R. This R code may not be particularly efficient or streamlined, but it works. 

We first describe the baseline scenario that we simulate, and then 5 variations of this scenario. Each simulation type has its own associated Snakefile, which should help you assess the workflow and run these simulations at scale. Sivan Yair has scripts that break down each simulation stage into smaller pieces, which may help with running these simulations on a cluster; contact Sivan if you would like access to those. 

### Baseline
Phenotypic stabilizing selection on a single trait. Mutations are assigned effect sizes from a normal distribution; effect sizes do not vary between populations. Effect sizes are stored as the selection coefficients of mutations because selection coefficients get saved in the simulation state recorded by SLiM. 

Relevant Parameters:
- strength of stabilizing selection (w)
- standard deviation of the mutant effect size distribution with mean zero (sd_a)

This scenario is the only one that we used to add directional selection to (see more details on that variation below). Because we could initiate a change in the optimum at any point in the divergence between populations, we recorded simulation data at multiple points leading up to the end of the generation. At all points where we recorded simulation output, we also wrote the standard deviation of the phenotype distribution to a text file, as this information is used to determine the extent of the optimum shift. 

### Neutral
We did not require the strength of stabilizing selection and did not calculate fitness every generation.

### Heavy Tail
This is the same as the baseline scenario, except we use a mutant effect size distribution with a heavy tail, produced as a mixture of three normal distributions. We choose to simulate a single mutant effect size distribution, and so do not use the standard deviation parameter in the baseline scenario.

### Directional Selection
This uses the burn-in from the baseline scenario. It continues from the population divergence of the baseline scenario. The additional parameter tells us the number of standard deviations of the phenotype that the optimum will move (num_sd_shift)

### GxE
Throughout the simulation, we record the effect of a mutation in three populations (the common ancestor and each descendant population). During the burn-in we use the effects of the common ancestor, and then these effects change during divergence depending on which population is being simulated. We record the effect of a mutation in each population in the "values" attribute of the mutation. This information is not saved by SLiM, and therefore we record this and read this from a separate text file. The additional parameter is the correlation of effect sizes between the three populations (corr). 

### Pleiotropy
When a mutation arises, we assign its effect on each trait independently from the same mutant effect size distribution. An individual's fitness is evaluated based on their Euclidean distance from the optimum. We record the effect of a mutation on each trait in the "values" attribute of the mutation. This information is not saved by SLiM, and therefore we record this and read this from a separate text file. The additional parameter is the number of traits under selection (n). 

## Incorporating Physical Linkage 
All of the above scenarios involve free recombination among QTLs that consist of 1bp loci. We also simulated QTLs of length 1kb, where there was low recombination within a QTL. Scripts for this scenario are provided. They are a replicate of the Baseline scenario, the only difference being the lines in the initialize() block that set up chromosome regions and recombination rates. Since these scripts are so similar to the Baseline, we do not provide a Snakefile, since just the input directory of scripts would need to change. 

# Processing output
We show how to process simulation output to get mean polygenic scores, Qst/Fst, and the reduction in prediction accuracy when using polygenic scores instead of additive genetic values, for different ascertainment schemes and generations. R scripts work with the output of all simulations of a particular variation to produce a single data frame that is saved as an RDS file and can be used elsewhere for further investigation.  This step is included in the Snakefiles. R scripts can be found in the Analysis Prep folder. The Snakefiles assume there's an existing directory for processed output that is called "Prepped_Data" with subdirectories, "Qst_Fst", "prediction_accuracy", and "mean_scores"



