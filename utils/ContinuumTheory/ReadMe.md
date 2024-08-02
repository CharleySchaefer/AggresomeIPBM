# Continuum Theory

## Thermodynamics

### Flory-Huggins Theory

The library in Thermodynamics/FloryHuggins/C includes functions written in C to calculate Flory-Huggins free energies and generate binary phase diagrams.

<img src="Thermodynamics/FloryHuggins/C/Demo/Fig_PhaseDiagramB.png" alt="PhaseDiagram" width="200"/>

The phase diagram was calculated using Demo\_BinaryFH.sh, and the raw data was stored in PhaseDiagram_50.0.out. The phase diagram above was generated in eps format using gnuplot with the script mkfig\_PhaseDiagramB.sh. The png version was generated using eps2png.sh. 

Using the same data in PhaseDiagram_50.0.out and the eps2png.sh script, one can generate information about droplet compaction using mkfig\_DropletCompaction.sh:


<img src="Thermodynamics/FloryHuggins/C/Demo/Fig_DropletCompaction.png" alt="DropletCompaction" width="200"/>

and about RNA partitioning from the dilute phase to the droplets using  mkfig\_RNA\_partitioning.sh:

<img src="Thermodynamics/FloryHuggins/C/Demo/Fig_RNA_partitioning.png" alt="RNA_partitioning" width="200"/>

### ScreenTape Analysis

In the ScreenTape\_analysis directory, the Flory-Huggins model for RNA partitioning is compared to measured ScreenTape data. See ReadMe.txt in that directory for details.


<img src="Thermodynamics/ScreenTape_analysis/Fig_tot_conc.png" alt="fit" width="200"/>
<img src="Thermodynamics/ScreenTape_analysis/Fig_RNA_length_cutoff30.png" alt="distribution" width="200"/>




## Dynamics

In the Dynamics directory the time-dependent RNA length distribution inside the droplets is calculated. See ReadMe.txt in that directory for details.

<img src="Dynamics/Fig_Distr.png" alt="timedistribution" width="200"/>
<img src="Dynamics/Fig_SizeVsTime.png" alt="timelength" width="200"/>


