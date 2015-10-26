//
// Writes Gadget binary file
//
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "msg.h"
#include "comm.h"
#include "write.h"
#include "gadget_file.h"

void write_snapshot(const char filebase[], Snapshot const * const snapshot,
		    int use_long_id)
{
  const double h= snapshot->h;
  char filename[256];
  sprintf(filename, "%s.%d", filebase, comm_this_node());


  FILE* fp= fopen(filename, "w");
  if(fp == 0)
    msg_abort(9000, "Error: Unable to write to file: %s\n", filename);

  ParticleMinimum* const p= snapshot->p;
  const int np= snapshot->np_local;
  const double boxsize= snapshot->boxsize;
  const double omega_m= snapshot->omega_m;

  if(use_long_id)
    msg_printf(normal, "Longid is used for GADGET snapshot. %d-byte.\n", 
	       sizeof(unsigned long long));
  else
    msg_printf(normal, "ID is %d-byte unsigned int\n", sizeof(unsigned int));


  long long np_send= np, np_total;
  MPI_Reduce(&np_send, &np_total, 1, MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Bcast(&np_total, 1, MPI_LONG_LONG, 0, MPI_COMM_WORLD);

  GadgetHeader header; assert(sizeof(GadgetHeader) == 256);
  memset(&header, 0, sizeof(GadgetHeader));

  const double rho_crit = 27.7455;
  const double m= omega_m*rho_crit*pow(boxsize, 3.0)/np_total;
  
  header.np[1]= np;
  header.mass[1]= m;
  header.time= snapshot->a;
  header.redshift= 1.0/header.time - 1;
  header.np_total[1]= (unsigned int) np_total;
  header.np_total_highword[1]= (unsigned int) (np_total >> 32);
  header.num_files= comm_nnode();
  header.boxsize= boxsize;
  header.omega0= omega_m;
  header.omega_lambda= 1.0 - omega_m;
  header.hubble_param= h;


  int blklen= sizeof(GadgetHeader);
  fwrite(&blklen, sizeof(blklen), 1, fp);
  fwrite(&header, sizeof(GadgetHeader), 1, fp);
  fwrite(&blklen, sizeof(blklen), 1, fp);

  // position
  blklen= np*sizeof(float)*3;
  fwrite(&blklen, sizeof(blklen), 1, fp);
  for(int i=0; i<np; i++)
    fwrite(p[i].x, sizeof(float), 3, fp);
  fwrite(&blklen, sizeof(blklen), 1, fp);

  // velocity
  const float vfac= 1.0/sqrt(snapshot->a); // Gadget convention

  fwrite(&blklen, sizeof(blklen), 1, fp);
  for(int i=0; i<np; i++) {
    float vout[]= {vfac*p[i].v[0], vfac*p[i].v[1], vfac*p[i].v[2]};
    fwrite(vout, sizeof(float), 3, fp);
  }
  fwrite(&blklen, sizeof(blklen), 1, fp);

  // id
  if(use_long_id) {
    blklen= np*sizeof(unsigned long long);
    fwrite(&blklen, sizeof(blklen), 1, fp);
    for(int i=0; i<np; i++) {
      unsigned long long id_out= p[i].id;
      fwrite(&id_out, sizeof(unsigned long long), 1, fp); 
    }
  }
  else {
    blklen= np*sizeof(unsigned int);
    fwrite(&blklen, sizeof(blklen), 1, fp);
    for(int i=0; i<np; i++) {
      unsigned int id_out= p[i].id;
      fwrite(&id_out, sizeof(unsigned int), 1, fp); 
    }
  }
  fwrite(&blklen, sizeof(blklen), 1, fp);


  
  fclose(fp);  

  msg_printf(normal, "snapshot %s written\n", filebase);
}

// Writing Gadget file for subsampled particles
void write_snapshot1(const char filename[], Snapshot const * const snapshot)
{
  FILE* fp= fopen(filename, "w");
  if(fp == 0)
    msg_abort(9000, "Error: Unable to write to file: %s\n", filename);

  ParticleMinimum* const p= snapshot->p;
  const int np= snapshot->np_local;
  const double boxsize= snapshot->boxsize;
  const double omega_m= snapshot->omega_m;

  GadgetHeader header; assert(sizeof(GadgetHeader) == 256);
  memset(&header, 0, sizeof(GadgetHeader));

  const double rho_crit = 27.7455;
  const double m= omega_m*rho_crit*pow(boxsize, 3.0)/np;
  
  header.np[1]= np;
  header.mass[1]= m;
  header.time= snapshot->a;
  header.redshift= 1.0/header.time - 1;
  header.np_total[1]= (unsigned int) np;
  header.np_total_highword[1]= 0;
  header.num_files= comm_nnode();
  header.boxsize= boxsize;
  header.omega0= omega_m;
  header.omega_lambda= 1.0 - omega_m;
  header.hubble_param= snapshot->h;


  int blklen= sizeof(GadgetHeader);
  fwrite(&blklen, sizeof(blklen), 1, fp);
  fwrite(&header, sizeof(GadgetHeader), 1, fp);
  fwrite(&blklen, sizeof(blklen), 1, fp);

  // position
  blklen= np*sizeof(float)*3;
  fwrite(&blklen, sizeof(blklen), 1, fp);
  for(int i=0; i<np; i++)
    fwrite(p[i].x, sizeof(float), 3, fp);
  fwrite(&blklen, sizeof(blklen), 1, fp);

  // velocity
  const float vfac= 1.0/sqrt(snapshot->a); // Gadget convention

  fwrite(&blklen, sizeof(blklen), 1, fp);
  for(int i=0; i<np; i++) {
    float vout[]= {vfac*p[i].v[0], vfac*p[i].v[1], vfac*p[i].v[2]};
    fwrite(vout, sizeof(float), 3, fp);
  }
  fwrite(&blklen, sizeof(blklen), 1, fp);

  // id
  blklen= np*sizeof(int);
  fwrite(&blklen, sizeof(blklen), 1, fp);
  for(int i=0; i<np; i++) {
    int id_out= p[i].id;
    fwrite(&id_out, sizeof(int), 1, fp); 
  }
  fwrite(&blklen, sizeof(blklen), 1, fp);
  
  fclose(fp);  

  msg_printf(normal, "subsample %s written\n", filename);
}


//
// Write particle force
//

void write_force(const char filebase[], Particles const * const particles)
{
  char filename[256];
  int inode= comm_this_node();

  sprintf(filename, "%s.%d", filebase, inode);
  FILE* fp= fopen(filename, "w");
  if(fp == 0) {
    msg_abort(9010, "Unable to write force to %s\n", filename);
  }

  Particle const * const p= particles->p;
  float3 * const f= particles->force;
  const int np= particles->np_local;

  fwrite(&np, sizeof(int), 1, fp);
  for(int i=0; i<np; i++) {
    int id= (int) p[i].id;
    fwrite(&id, sizeof(int), 1, fp);
    fwrite(p[i].x, sizeof(float), 3, fp);
    fwrite(f[i], sizeof(float), 3, fp);
  }
  fwrite(&np, sizeof(int), 1, fp);

  fclose(fp);
}

// Writes binary file for subsampled particles. Old version, not in use.
void write_particles_binary(const char filename[], Snapshot const * const snapshot)
{
  FILE* fp= fopen(filename, "w");
  if(fp == 0)
    msg_abort(9000, "Error: Unable to write to file: %s\n", filename);

  ParticleMinimum* const p= snapshot->p;
  const int np= snapshot->np_local;
  const float boxsize= snapshot->boxsize; 
  const float omega_m= snapshot->omega_m;
  const double rho_crit = 27.7455;
  const float m= omega_m*rho_crit*pow(boxsize, 3.0)/np;  
  const float redshift= 1.0/snapshot->a - 1.0;

  // Header 6 floats
  fwrite(&snapshot->boxsize, sizeof(float), 1, fp);
  fwrite(&m, sizeof(float), 1, fp);
  fwrite(&snapshot->omega_m, sizeof(float), 1, fp);
  fwrite(&snapshot->h, sizeof(float), 1, fp);
  fwrite(&snapshot->a, sizeof(float), 1, fp);
  fwrite(&redshift, sizeof(float), 1, fp);

  fwrite(&np, sizeof(int), 1, fp);

  // positions, velocities
  for(int i=0; i<np; i++) {
    fwrite(p[i].x, sizeof(float), 3, fp);
    fwrite(p[i].v, sizeof(float), 3, fp);
  }

  fwrite(&np, sizeof(int), 1, fp);
  
  int ret= fclose(fp); assert(ret == 0);

  msg_printf(normal, "subsample binary %s written\n", filename);
}

// Writes binary file for subsampled particles
// Newer version using Parallel write MPI_File_write_at
void write_particles_binary_mpi(const char filename[], ParticleSubsample const * const p, const int np, float const * const header)
{
  //
  // Gather number of particles to compute the offset for writing
  // 
  const int this_node= comm_this_node();
  const int nnode= comm_nnode();

  int* const np_local= malloc(sizeof(int)*nnode); assert(np_local);

  int ret= 
    MPI_Gather(&np, 1, MPI_INT, np_local, 1, MPI_INT, 0, MPI_COMM_WORLD);
  assert(ret == MPI_SUCCESS);

  int np_total= 0; // Total number of subsample particles
  if(this_node == 0) {
    for(int i=0; i<nnode; i++) {
      int np_i= np_local[i];
      np_local[i]= np_total;
        // np_local <- partial_sum = Sum_{inode}^{this node - 1} np
      np_total  += np_i;
    }
  }

  int np_partial_sum= 0;
  ret= MPI_Scatter(np_local, 1, MPI_INT,
		   &np_partial_sum, 1, MPI_INT, 0, MPI_COMM_WORLD);
  assert(ret == MPI_SUCCESS);

  const size_t nfloat_header= 6;
  size_t offset= sizeof(float)*nfloat_header + sizeof(int) +
                 sizeof(ParticleSubsample)*np_partial_sum;

  // Write particles to a file
  MPI_File fh;
  MPI_Status stat;
  ret= MPI_File_open(MPI_COMM_WORLD, filename,
		     MPI_MODE_CREATE | MPI_MODE_WRONLY,
		     MPI_INFO_NULL, &fh); assert(ret == MPI_SUCCESS);

  if(this_node == 0) {
    // header: boxsize, m_particle, omega_m, h, a, redshift & np_total
    ret= MPI_File_write_at(fh, 0, header, sizeof(float)*nfloat_header,
		      MPI_BYTE, &stat); assert(ret == MPI_SUCCESS);
    ret= MPI_File_write_at(fh, sizeof(float)*nfloat_header, &np_total,
		     sizeof(int), MPI_BYTE, &stat); assert(ret == MPI_SUCCESS);
  }

  assert(sizeof(ParticleSubsample) % sizeof(float) == 0);
  MPI_File_write_at(fh, offset, p, sizeof(ParticleSubsample)*np,
		    MPI_BYTE, &stat);

  if(this_node == 0) {
    // footer: np_total
    size_t offset_footer= sizeof(float)*nfloat_header + sizeof(int) +
                          sizeof(ParticleSubsample)*np_total;
    ret= MPI_File_write_at(fh, offset_footer, &np_total, sizeof(int),
			   MPI_BYTE, &stat); assert(ret == MPI_SUCCESS);
  }

  ret= MPI_File_close(&fh); assert(ret == MPI_SUCCESS);
  
  msg_printf(verbose, "Subsampled particles %d written\n", np_total);
}
