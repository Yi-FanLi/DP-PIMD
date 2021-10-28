# A Document on the fix style dp_pimd in LAMMPS

## 1. Input script keywords
- method: the method to use
  - pimd: pimd in the Cartesian coordinate (temporarily disabled, you should not use this)
  - nmpimd (default): pimd in the normal mode coordinate (only this is supported)
  - cmd: centroid molecular dynamics (disabled, not supported)
- integrator: the splitting scheme of the integrator
  - obabo (default): OBABO scheme (please use this currently)
  - baoab: BAOAB scheme
    - O stands for thermostatting step; 
    - B stands for exerting the forcefield; 
    - A stands for the free ring polymer evolution.
- ensemble: the MD simulation ensemble
  - nve: the microcanonical ensemble
  - nvt (default): the canonical ensemble
  - nph: the isobaric-isoenthalpic ensemble
  - npt: the isothermal-isobaric ensemble
- temp (default is 298.15): the target temperature to be sampled
- thermostat: the thermostat to control the temperature
  - PILE_L (default): the local Langevin thermostat
  - SVR: the stochastic velocity rescaling thermostat
  - PILE_G: the global Langevin thermostat
- tau (default is 1.0): the damping time of the Lagevin thermostat coupled to the centroid mode
- scale (default is 1.0): the scaling factor of the damping time of the Langevin thermostats coupled to the non-centroid modes (damping time = characteristic time / scale)
- press: the external pressure exterted via the barostat
- taup: the damping time of the barostat
- fmmode: the way of choosing fictitious mass
  - physical (default): the physical mass of the atoms
  - normal: the pnysical mass times eigenvalues of each normal mode

## 2. Output properties
- 1: normal mode kinetic energy, prints the kinetic energy of one indivisual normal mode in each log file
- 2: normal mode spring elastic energy, prints the ring polymer elastic energy of one indivisual normal mode in each log file
- 3: average potential energy over all beads
- 4: total ring polymer Hamiltonian
- 5: primitive quantum kinetic energy estimator (smae as kinetic_td in i-PI)
- 6: virial quantum kinetic energy estimator
- 7: centroid-virial quantum kinetic energy estimator (same as kinetic_cv in i-PI)
- 8: primitive pressure estimator
- 9: virial pressure estimator
- 10: centroid-virial presure estimator (same as pressure_cv in i-PI)
- 11: velocity of the volume degree of freedom; only useful in the NPH and NPT ensembles
- 12: total enthalpy, Hamiltonian + PV