cola_halo
=========

Parallel COLA cosmological simulation + 2LPT initial condition
generator + FoF halo finder

This is a parallel cosmological N-body simulation code running the
following on the fly, designed to generate hundreds of realizations:

1. Generates random Gaussian initial condition with 2LPT --
second-order Lagrangian Perturbation Theory (Scoccimarro 1998; Crocce,
Pueblas, & Scoccimarro 2006,[[1]]).
2. Time evolve N-body particles with the COmoving Lagrangian
Acceleration (COLA) method (Tassev, Zaldarriaga, & Eisenstein
2013,[[2]])
3. Finds dark-matter haloes with the Friends-of-Friends algorithm
(Davis et al. 1985, [Washington University N-body shop][3]).


This code is based on publicly available codes ([1], [2], [3]),
modified and redistributed under GPLv3.

# Compile

This code requires an MPI C compiler `mpicc`, GNU Scientific Library
[GSL](https://www.gnu.org/software/gsl/) and
[FFTW](http://www.fftw.org/). If they are installed in nonstandard
locations, set the variables `GSL_DIR` and `FFTW3_DIR` in the
`Makefile`.

```bash
% make
```

should compile the code and create the executable `cola_halo`.

# Run


- The simulation requires a linear **power spectrum** extrapolated to
  redshift z=0, specified by the `powerspectrum` variable in the
  parameter file `param.lua`. The file should give a list of *k*
  [h/Mpc] and *P(k)* [(Mpc/h)<sup>3</sup>], which can be the
  `matterpower.dat` file generate by [CAMB](http://camb.info/). The
  power spectrum should have the &sigma;<sub>8</sub>; specified in the
  parameter file (the code does not rescale the amplitude to the given
  &sigma;<sub>8</sub>). A power spectrum file consistent with the
  default parameter file is
  [here](http://www13273u.sakura.ne.jp/wizcola/camb0_matterpower.dat).

- Edit the parameter file `param.lua` and set the number of particles,
  box size, number of time steps, etc.

- Run the code with the name of the parameter file

```
   mpirun -n 4 cola_halo param.lua
```
Where 4 is the number of MPI nodes you use. It can be any number, but
we recommend it is a divisor of the number of particle per dimension
`nc` for effective domain decomposition.

# Outputs

File names are formatted as fof00100a.txt, where 00100 is the
`random_seed` for the initial condition, and a,b,c,... corresponds to
first, second, third, ..., `output_redshifts`.

- The N-body particles `snp` are in the
  [GADGET](http://www.mpa-garching.mpg.de/gadget/) binary format.
- The halo catalogue `fof` is a ascii file:

```
nfof x y z vx vy vz
```

where `nfof` is the number of N-body particles in the FoF group, (x y
z) is the centre-of-mass position in comoving 1/h Mpc, (vx vy vz) is
the average peculiar (physical) velocity in km/s.

- The matter density field `grid` and subsampled particles `sub` are
  simple binary files; see `tools/density_grid.c` and
  `tools/subsampled_particles.c`.

# Options

- You can do simple calculation in the parameter file, e.g,

```
coarse_grid_nc = nc/4
```

with the [LUA programming launguage](http://www.lua.org/) by turning
on the `USE_LUA` option in the parameter file.

- You can turn on Hybrid MPI+OpenMP parallelisation by setting the
  `OPENMP` variable. The FoF halo finder, however, is not parallelised
  by OpenMP.

# Related codes

- l-picola: Another parallel COLA implementation by Cullan Howlett and
  Marc Manera [4], including on-the-fly lightcone generation.

- QrPM: A fork of this code. The parallelisation improved by using 2D
  domain decomposition. [5]


# Authors

This code is assembled and parallelized (COLA & FOF) by Jun Koda. See
reference 5 below about the application to the mock generation of
WiggleZ Dark Energy Survey and its accuracy. Feel free to [open an
issue](https://github.com/junkoda/cola_halo/issues) for bug reports or
questions.

# References 
1. Crocce, M., Pueblas, S., & Scoccimarro, R. 2006, MNRAS, 373, 369
2. Davis, M., Efstathiou, G., Frenk, C.~S., & White, S.D.M. 1985, ApJ, 292, 371 
3. Scoccimarro, R. 1998, MNRAS, 299, 1097
4. Tassev, S., Zaldarriaga, M., & Eisenstein, D.~J. 2013, JCAP, 6, 36 
5. Koda, J., et al 2015, [arXiv:1507.05329](http://arxiv.org/abs/1507.05329), MNRAS submitted

[1]: http://cosmo.nyu.edu/roman/2LPT/
[2]: https://bitbucket.org/tassev/colacode/
[3]: http://www-hpcc.astro.washington.edu/tools/fof.html
[4]: https://github.com/CullanHowlett/l-picola
[5]: https://github.com/rainwoodman/QrPM
