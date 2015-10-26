#ifndef PARTICLE_H
#define PARTICLE_H

typedef float float3[3];

typedef struct {
  float x[3];
  float dx1[3]; // ZA displacement
  float dx2[3]; // 2LPT displacement
  float v[3];   // velocity
  long long id;
} Particle;

typedef struct {
  Particle* p;
  float3* force;
  float a_x, a_v;

  int np_local, np_allocated;
  long long np_total;
  float np_average;
} Particles;

typedef struct {
  float x[3];
  float v[3];
  long long id;
} ParticleMinimum;

typedef struct {
  ParticleMinimum* p;
  int np_local;
  int np_allocated;
  long long np_total;
  float np_average;
  float a;
  float boxsize;
  int nc;
  float omega_m, h;
  int seed;
  char* filename;
} Snapshot;

typedef struct {
  float x[3];
  float v[3];
  //float f[3];
} ParticleSubsample;

#endif
