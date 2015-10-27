#ifndef WRITE_H
#define WRITE_H 1

#include "particle.h"

//void write_snapshot(const char filebase[], Snapshot const * const snapshot,
//		    double h, int use_longid);
void write_snapshot(const char filebase[], Snapshot const * const snapshot,
		    int use_longid);
void write_snapshot1(const char filename[], Snapshot const * const snapshot);
//void write_snapshot(const char filebase[], Snapshot const * const snapshot,
//		    double boxsize, double omega_m, double h);

void write_particles_binary(const char filename[], Snapshot const * const snapshot);
void write_particles_binary_mpi(const char filename[], ParticleSubsample const * const p, const int np, float const * const header);

void write_force(const char filebase[], Particles const * const particles);
#endif
