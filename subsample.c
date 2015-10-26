//
// Subsample the N-body particle and write to file 
//
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <mpi.h>
#include <gsl/gsl_rng.h>

#include "particle.h"
#include "msg.h"
#include "comm.h"
#include "write.h"

static gsl_rng* Random_generator= 0;
static double SubsampleFactor= 0.0;


void subsample_init(const double subsample_factor, const unsigned int seed)
{
  Random_generator = gsl_rng_alloc(gsl_rng_ranlxd1);
  SubsampleFactor= subsample_factor;

  const int this_node= comm_this_node();
  gsl_rng_set(Random_generator, 2*seed + 100*this_node);

  msg_printf(verbose, "Subsampling initialized. Factor= %le, seed= %u\n", 
	     subsample_factor, 2*seed);  
}

void subsample_finalize(void)
{
   gsl_rng_free(Random_generator);
}


// Oldversion: collect subsample to node 0 and then write to file
// newer version below (write_random_subsample()) uses MPI_File_write
void write_random_sabsample1(const char filename[], Snapshot const * const snapshot, void* const mem, size_t mem_size)
{
  ParticleMinimum const * const p= snapshot->p;
  const int np= snapshot->np_local;
  const int nbuf = mem_size/sizeof(ParticleMinimum);
  
  if(SubsampleFactor*np + 5*sqrt(SubsampleFactor*np) > nbuf) { 
    msg_abort(9300, 
	      "Not enough memory for local subsampling ~ %.2lf particles\n",
	      SubsampleFactor*np);
  }

  ParticleMinimum* out= (ParticleMinimum*) mem;
  int nsub= 0; // number of subsampled particles in this node
  for(int i=0; i<np; i++) {
    if(gsl_rng_uniform(Random_generator) < SubsampleFactor)
      out[nsub++]= p[i];
  }

  //
  // Prepare to gather particles. Get numbers.
  // 
  const int this_node= comm_this_node();
  const int nnode= comm_nnode();

  int* const nsub_recv= malloc(sizeof(int)*nnode*2); assert(nsub_recv);
  int* disp= nsub_recv + nnode;
  int ret= 
    MPI_Gather(&nsub, 1, MPI_INT, nsub_recv, 1, MPI_INT, 0, MPI_COMM_WORLD);
  assert(ret == MPI_SUCCESS);

  ParticleMinimum* buf1= 0;
  int ns= 0; // Total number of subsample particles
  if(this_node == 0) {
    for(int i=0; i<nnode; i++) {
      disp[i] = ns*sizeof(ParticleMinimum);
      ns  += nsub_recv[i];
      nsub_recv[i] *= sizeof(ParticleMinimum);
    }
    if(nbuf < nsub+ns)
      msg_abort(9100, "Not enough space to gather subsample particles %d (local) + %d (global)\n", nsub, ns);

    buf1= out+nsub; // Buffer for all particles

    msg_printf(verbose, "Subsampled particles %d (average %.2lf)\n", 
	       ns, SubsampleFactor*snapshot->np_total);
  }

    
  // Gather subsample particles to node 0
  ret= MPI_Gatherv(out, nsub*sizeof(ParticleMinimum), MPI_BYTE,
		   buf1, nsub_recv, disp, MPI_BYTE,
		   0, MPI_COMM_WORLD);
  assert(ret == MPI_SUCCESS);
  
  if(this_node == 0) {
    Snapshot subsample;
    subsample.p= buf1;
    subsample.np_local= ns;
    subsample.np_allocated= ns;
    subsample.np_total= ns;
    subsample.np_average= ns;
    subsample.a= snapshot->a;
    subsample.nc= 0;
    subsample.boxsize= snapshot->boxsize;
    subsample.omega_m= snapshot->omega_m;
    subsample.h= snapshot->h;
    subsample.seed= snapshot->seed;

    write_particles_binary(filename, &subsample);
    //write_snapshot1(filename, &subsample); // Gadget Format
  }
  free(nsub_recv);
}

void write_random_sabsample(const char filename[], Snapshot const * const snapshot, void* const mem, size_t mem_size)
{
  // Assumption: Number of subsampled particles are < max(int) (2^31 for 4-byte)
  ParticleMinimum const * const p= snapshot->p;
  const int np= snapshot->np_local;
  const int nbuf = mem_size/sizeof(ParticleSubsample);
  
  if(SubsampleFactor*np + 5*sqrt(SubsampleFactor*np) > nbuf) { 
    msg_abort(9300, 
	      "Not enough memory for local subsampling ~ %.2lf particles\n",
	      SubsampleFactor*np);
  }

  ParticleSubsample* out= (ParticleSubsample*) mem;
  int nsub= 0; // number of subsampled particles in this node
  for(int i=0; i<np; i++) {
    if(gsl_rng_uniform(Random_generator) < SubsampleFactor) {
      out[nsub].x[0]= p[i].x[0];
      out[nsub].x[1]= p[i].x[1];
      out[nsub].x[2]= p[i].x[2];
      out[nsub].v[0]= p[i].v[0];
      out[nsub].v[1]= p[i].v[1];
      out[nsub].v[2]= p[i].v[2];
      //out[nsub].f[0]= p[i].f[0];
      //out[nsub].f[1]= p[i].f[1];
      //out[nsub].f[2]= p[i].f[2];
      nsub++;
    }
  }

  // Subsample file header
  const double rho_crit = 27.7455;
  const float redshift= 1.0/snapshot->a - 1.0;

  float header[6];
  header[0]= snapshot->boxsize;
  header[1]= snapshot->omega_m*rho_crit*pow(snapshot->boxsize, 3.0)/np;
             // particle mass
  header[2]= snapshot->omega_m;
  header[3]= snapshot->h;
  header[4]= snapshot->a;
  header[5]= redshift;

  write_particles_binary_mpi(filename, out, nsub, header);
}


// Regular subsampling (not used, but has a benefit of not having shot noise)

/*
void write_subsample(const char filename[], const int fac, Snapshot const * const snapshot, void* const mem, size_t mem_size)
{
  // Subsample particles regularly and write as Gadget file
  // fac: resampling factor per dimension
  ParticleMinimum const * const p= snapshot->p;
  const int np= snapshot->np_local;
  const int nc= snapshot->nc;

  ParticleMinimum* out= (ParticleMinimum*) mem;
  int nbuf = mem_size/sizeof(ParticleMinimum);
  int nsub= 0; // number of subsampled particles in this node

  for(int i=0; i<np; i++) {
    // subsampling every "fac" particles per dimension
    long long n= p[i].id - 1;   // id is starting with 1
    int nz= n % nc;
    int ny= (n/nc) % nc;
    int nx= n/(nc*nc);
    if(nx % fac == 0 && ny % fac == 0 && nz % fac == 0)
      out[nsub++]= p[i];
  }

  //
  // Prepare to gather particles. Get numbers.
  // 
  const int this_node= comm_this_node();
  const int nnode= comm_nnode();

  int* const nsub_recv= malloc(sizeof(int)*nnode*2); assert(nsub_recv);
  int* disp= nsub_recv + nnode;
  int ret= 
    MPI_Gather(&nsub, 1, MPI_INT, nsub_recv, 1, MPI_INT, 0, MPI_COMM_WORLD);
  assert(ret == MPI_SUCCESS);

  ParticleMinimum* buf1= 0;
  int ns= 0; // Total number of subsample particles
  if(this_node == 0) {
    for(int i=0; i<nnode; i++) {
      disp[i] = ns*sizeof(ParticleMinimum);
      ns  += nsub_recv[i];
      nsub_recv[i] *= sizeof(ParticleMinimum);
    }
    if(nbuf < nsub+ns)
      msg_abort(9100, "Not enough space for subsample particles %d (local) + %d (global)\n", nsub, ns);


    buf1= out+nsub; // Buffer for all particles

  msg_printf(verbose, "Subsampled particles %d (%d)\n", ns, (nc/fac)*(nc/fac)*(nc/fac));

    // This can be dropped if nc % fac != 0 is intentional
    assert(ns == (nc/fac)*(nc/fac)*(nc/fac)); 
  }


  // Gather subsample particles to node 0
  ret= MPI_Gatherv(out, nsub*sizeof(ParticleMinimum), MPI_BYTE,
		   buf1, nsub_recv, disp, MPI_BYTE,
		   0, MPI_COMM_WORLD);
  assert(ret == MPI_SUCCESS);

  if(this_node == 0) {
    Snapshot subsample;
    subsample.p= buf1;
    subsample.np_local= ns;
    subsample.np_allocated= ns;
    subsample.np_total= ns;
    subsample.np_average= ns;
    subsample.a= snapshot->a;
    subsample.nc= nc/fac;
    subsample.boxsize= snapshot->boxsize;
    subsample.omega_m= snapshot->omega_m;
    subsample.h= snapshot->h;
    subsample.seed= snapshot->seed;

    write_snapshot1(filename, &subsample);
  }

  free(nsub_recv);
}
*/
