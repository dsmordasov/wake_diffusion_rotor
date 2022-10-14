A suite of Python scientific computing tool scripts running and postprocessing PyWakeEllipSys (PWES) and HAWC2S simulations, investigating the 'wake-diffusion' rotor concept for my MSc thesis. 

You can download the thesis in its entirety here: [Numerical investigation of the wake-diffusion rotor concept](http://resolver.tudelft.nl/uuid:f4a34434-a186-4127-aedd-d9d8f5e62c2a)


Tl;dr changing the geometry of the blade can, at the cost of individual wind turbine power generation loss, produce higher power generation over multiple wind turbines.

### Abstract
With the growing global utilisation of wind energy, lowering of its levelised cost of energy is pursued. This effort is hindered by wind turbine wakes and their detrimental effects on wind farm profitability. Most currently researched wake mitigation methods are based on wind farm system-level interventions, and there is a research gap on individual wind turbine design wake alleviation measures.

This thesis investigates the wake-diffusion rotor concept, a wind turbine design with an outboard shifted thrust distribution along the blade, and its effects on wake mitigation and power generation in wind farms.
Three wind turbine rotors with equivalent thrust coefficient are modelled as actuator disks in a RANS CFD software PyWakeEllipSys, corresponding to an inboard shifted thrust distribution, a conventional thrust distribution, and an outboard shifted thrust distribution representing the wake-diffusion rotor concept. Investigated scenarios include a single wind turbine and rows of three and ten wind turbines aligned to the inflow.
Compared to the baseline case, outboard thrust distribution produces lower wake velocity deficits, increased momentum transfer and increased turbulence generation. The wake-diffusion rotor concept increases the wind farm power generation by 3.8 % and 2.5 % in the three and ten wind turbine row scenarios respectively, with the largest relative gain in wind turbine power of 13.9 % for the second wind turbine in the row. Consecutive wind turbines in the row experience diminishing gains. The inboard shifted thrust distribution had opposite effects for all of these aspects.
The benefits in wake mitigation and power generation in wind farms from using the wake-diffusion rotor concept are concluded to come from higher turbulence mixing caused by the outboard shifted thrust distribution. Recommendations on its use for maximum effectiveness are given. Future research could focus on the design implementation of this concept, investigate combinations with other wake mitigation methods and perform parametric wind farm AEP studies.

### Engineering software used
PyWakeEllipSys is a CFD package software, which combines a RANS based wind farm flow model with PyWake, an open source AEP calculator, with the RANS model being based on the general purpose flow solver EllipSys3D. This is a Fortran/MPI based closed source licensed software, the development of which was initiated at DTU, and now continued at DTU Wind and Energy Systems. 

HAWCStab2 is a frequency based aeroservoelastic code for steady states computation, most often used for stability analysis of wind turbines, developed at DTU. HAWC2S is a command line version of this code

### Script descriptions
`config.py` contains quality of life imports, visualising scripts and settings for running PyWakeEllipSys simulations on the DTU HPC cluster and local post-processing scripts.

`h2s.py` contains the same for local HAWC2S simulations and posprocessing.

`blade_design_tool.py` contains code that allows for local wind turbine blade design, and subsequent evaluation using HAWC2S.


Files beginning with `pwes` are PyWakeEllipSys simulation scripts to be ran on the DTU HPC cluster.

Files beginning with `pp` are local post-processing scripts, and (most often) require the data from the DTU HPC cluster simulations. This allows for analysis and visualisation of the produced CFD and BEM data. The graphs producible with these codes include: 
- wind turbine wake velocity deficits
- turbulence intensity
- lateral and vertical momentum transport 
- power production along inflow-aligned wind turbine rows
- load distributions along the wind turbine blade for varying configurations
- ... and more, as can be viewed in the thesis. 

Big thanks to my thesis co-supervisor Mads Christian Baungaard for providing me with his highly relevant code regarding PyWakeEllipSys CFD data post-processing.


