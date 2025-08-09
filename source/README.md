# VISKY-MD (V1.1).

## Features:

- Harmonic bonds
- Harmonic angles
- Lennar-Jones pair interactions
- DPD
- Hybrid particle-field molecular dynamics 
-> with classical interactions and Gaussian convolution of atomic denisities with FFT
- Shear flow (Non-eq.)
- Thermostats: Berendsen, Andersen, Lowe-Andersen
- Multi particle collision dynamics

### In progress:

- Different particle types

## Tools

- RDF
- MSD
- E2E autocorrelation & relaxation time (polymers)

## Running the code:

### Julia packages to install before running:
import Pkg; Pkg.add("Formatting"); Pkg.add("Distributions"); Pkg.add("ProgressMeter"); Pkg.add("Debugger"); Pkg.add("Statistics"); Pkg.add("LinearAlgebra"); Pkg.add("Random"); Pkg.add("FFTW")

### Runing the code.
You will need an input file, typically called input.data and a control file, typically named control.txt:

- input.data
- control.txt

Run with the following command:

julia /your/path/VISKY/source/start.jl
