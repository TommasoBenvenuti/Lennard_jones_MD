# Lennard_jones_MD 


## Purpose
The code is intended to run classical molecular dynamics (MD) simulations using the Lennard-Jones potential to study multi-component mixture (i.e. methan-hydrogen), both in NVE and NVT ensemble
It produces quantities such as radial distribution function, velocity auto and cross correlation function, Spectra of VACF.
Bussi-Parrinello-Donadio Thermostat has been implemented and tested. Equations of motion are solved by position Verlet algorithm.


## Main files and routine descriptions
- `input.f90`
  - Handles reading simulation parameters and general configuration.
  - Contains logic for reading parameters from `input.txt` (e.g., number of particles, density, target temperature, time step, number of steps, etc.).

- `input.txt`
  - Example input file with the numeric parameters required by the code.

- `Inizializza.f90`
  - Builds the initial  particle random configuration (positions) inside a periodic simulation box
  - Sets velocities accordingly to Boltzmann distribution.

- `Forze.f90`
  - Computes  forces using the Lennard-Jones potential.
  - Implements periodic boundary conditions and a potential cutoff.
  - Returns potential energy and forces on each particle.

- `Integrazione.f90`
  - Contains the time integration algorithm (e.g., Verlet) to update positions and velocities.
  - Also handles kinetic and total energy calculation at each step, and gathers data for analysis and output.

- `BPDThermos.f90 
  - Contains routines to control the system temperature and apply velocity scaling during the simulation (Bussi Parrinello Donadio)

- `Scaling.f90`
  - Utility routines for scaling positions for a high density system

- `Spectrum.f90`
  - Routines for spectral analysis: computing a power spectrum of VACF from time series produced by the simulation.

- `Vacfblocking.f90`
  - Computes the velocity autocorrelation function (VACF) and Cross correlation using optimization techniques.

## Input format
See `input.txt` for a concrete example. Typical parameters (as read by `input.txt`) may include:
- N: number of particles
- density (rho) in reduced units 
- initial temperature 
- dt: time step
- nsteps: number of integration steps
- ensemble (NVT or NVE)
- chemical species with the number of particles for each species 

Check `input.f90` for the exact parameter names and ordering used by this implementation.

## Expected output
Depending on the configuration, the code typically produces:
- Time series files of total energy, potential energy, kinetic energy, and temperature.
- Files for VACF and spectral data (e.g., time vs VACF, or frequency vs spectral amplitude).
- Optional snapshot files with positions/velocities for visualization.
Check the output/print statements in `Integrazione.f90`, `Vacfblocking.f90`, and `Spectrum.f90` for the exact filenames and formats.

## How to use 

Compile using fopenmp 

