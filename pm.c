//
// Particle Mesh gravitational force calculation
//
// This code is a modification to the original serial COLA code
// by Svetlin Tassev. See below.
//

/*
    Copyright (c) 2011-2013       Svetlin Tassev
                           Harvard University, Princeton University
 
    This file is part of COLAcode.

    COLAcode is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    COLAcode is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with COLAcode.  If not, see <http://www.gnu.org/licenses/>.
*/


// This file contains some standard functions for a PM code. 
// Nothing COLA-specific.

#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <fftw3-mpi.h>
#include "pm.h"
#include "msg.h"
#include "comm.h"
#include "timer.h"

#ifndef M_PI
#define M_PI 3.1415926535897932384626433832795
#endif

static int ThisNode, NNode, LeftNode=-1, RightNode=-1;

static int Ngrid;
static size_t NgridL; // for index to avoid 4byte integer overflow
static int PM_factor;
static int Local_nx, Local_x_start;
static int Local_ny_td, Local_y_start_td;  // transposed
static float BoxSize;

static fftwf_complex* fftdata;
static fftwf_plan p0, p11;
static fftwf_complex* density_k;


static void *BufMove;
static int   NMoveBuffer;

static long long NParticleTotal;

typedef struct {
  int nbuf;
  int nsend, nrecv;
  float *vec, *vrecv;
  int* index;
} BufferVec3;

static BufferVec3 BufPos;

static inline void WRtPlus(float * const d, 
		    const int i, const int j, const int k, const float f)
{
#ifdef _OPENMP
  #pragma omp atomic
#endif
  d[k + 2*(NgridL/2 + 1)*(j + NgridL*i)] += f;
}

static inline float REd(float const * const d, const int i, const int j, const int k)
{
  return d[k + 2*(NgridL/2 + 1)*(j + NgridL * i)];
}

void pm_init(const int nc_pm, const int nc_pm_factor, const float boxsize,
	     const float np_alloc_factor,
	     void* const mem1, const size_t size1,
	     void* const mem2, const size_t size2)
{
  // Assume FFTW3 is initialized in lpt.c

  Ngrid= nc_pm; NgridL= nc_pm;
  BoxSize= boxsize;
  PM_factor= nc_pm_factor;

  ptrdiff_t local_nx, local_x_start, local_ny, local_y_start;
  ptrdiff_t total_size= 
    fftwf_mpi_local_size_3d_transposed(Ngrid, Ngrid, Ngrid/2+1, MPI_COMM_WORLD,
	                 &local_nx, &local_x_start, &local_ny, &local_y_start);

  Local_nx= local_nx;
  Local_x_start= local_x_start;
  Local_ny_td= local_ny;
  Local_y_start_td= local_y_start;

  size_t bytes= sizeof(fftwf_complex)*total_size;

  if(mem1 == 0) {
    fftdata=  fftwf_alloc_complex(total_size);
    msg_printf(info, 
	       "%d Mbytes allocated for particle mesh (density)\n", 
	       (int)(bytes/(1024*1024)));
  }
  else {
    assert(size1 >= total_size*sizeof(fftwf_complex));
    fftdata= (fftwf_complex*) mem1;
  }

  if(mem2 == 0) {
    density_k= fftwf_alloc_complex((NgridL/2+1)*NgridL*Local_ny_td);
  }
  else {
    assert(size2 >= sizeof(fftwf_complex)*(NgridL/2+1)*NgridL*Local_ny_td);
    density_k= (fftwf_complex*) mem2;
  }

  if(fftdata == 0)
    msg_abort(4001, "Unable to allocate memory for particle mesh: %d Mbytes\n",
	      (int)(2*bytes/(1024*1024)));

  double pm_index_max= (double) Ngrid*Ngrid*(Ngrid/2+1);
  if(pm_index_max > pow(2.0, 31))
    msg_printf(info, 
	       "Local number of PM Mesh exceeds 4-byte integer limit %ld\n", 
	       2*NgridL*NgridL*(NgridL/2+1));

  msg_printf(verbose, "Setting up FFTW3 plans (FFTW_MEASURE)\n");

  p0=  fftwf_mpi_plan_dft_r2c_3d(Ngrid, Ngrid, Ngrid, 
				 (float*) fftdata, fftdata,
		    MPI_COMM_WORLD, FFTW_MEASURE | FFTW_MPI_TRANSPOSED_OUT);
  
  // inverse FFT
  p11= fftwf_mpi_plan_dft_c2r_3d(Ngrid, Ngrid, Ngrid, fftdata, (float*)fftdata,
		     MPI_COMM_WORLD, FFTW_MEASURE | FFTW_MPI_TRANSPOSED_IN);

  // Find Adjacent Nodes
  MPI_Comm_rank(MPI_COMM_WORLD, &ThisNode);
  MPI_Comm_size(MPI_COMM_WORLD, &NNode);

  int* local_x_table  = malloc(sizeof(int)*NNode*2);
  int* local_nx_table = local_x_table + NNode;

  MPI_Allgather(&Local_nx, 1, MPI_INT, local_nx_table, 1, MPI_INT, 
		MPI_COMM_WORLD);
  MPI_Allgather(&Local_x_start, 1, MPI_INT, local_x_table, 1, MPI_INT, 
		MPI_COMM_WORLD);

  for(int i=0; i<NNode; i++) {
    if((local_x_table[i] + local_nx_table[i]) % Ngrid == Local_x_start)
      LeftNode= i;
    if((Local_x_start + Local_nx) % Ngrid == local_x_table[i])
      RightNode= i;
  }
  assert(LeftNode >= 0 && RightNode >= 0);
  
  free(local_x_table);

  // workspace for move_particles
  BufMove= mem1;
  NMoveBuffer= size1/sizeof(Particle);
  msg_printf(debug, "NMove buffer %d\n", NMoveBuffer);

  // buffer for particle transfer
  const int ncp= nc_pm/nc_pm_factor;
  NParticleTotal= (long long) ncp*ncp*ncp;

  // If you encounter "sendrecv_buffer" shortage,
  // this number should be increased
  const int nbuf=
    (int) ceil(np_alloc_factor*ncp*ncp + 5*sqrt(np_alloc_factor*ncp*ncp));

  BufPos.nbuf= nbuf;
  BufPos.vec  = malloc(sizeof(float)*3*nbuf);
  BufPos.vrecv= malloc(sizeof(float)*3*nbuf);
  BufPos.index= malloc(sizeof(int)*nbuf);

  if(BufPos.vec == 0 || BufPos.vrecv == 0 || BufPos.index == 0)
    msg_abort(4002, "Error: Unable to allocate memory for PM buffer");

  msg_printf(info, "%d Mbytes allocated for particle buffer\n",
	     (sizeof(float)*6 + sizeof(int))*nbuf/(1024*1024));

  //NP_Average= (float) pow((double)nc_pm/nc_pm_factor, 3.0)/NNode;

}

void pm_finalize(void)
{
  fftwf_free(fftdata);

  fftwf_destroy_plan(p0);
  fftwf_destroy_plan(p11);
}


void PtoMesh(const Particle P[], const int np, float* const density)
{
  // particles are assumed to be periodiclly wraped up in y,z direction

  msg_printf(verbose, "Calculating PtoMesh\n");

  const float scaleBox=((float) Ngrid)/((float) BoxSize);
  const float WPAR= pow(PM_factor, 3);

#ifdef _OPENMP
  #pragma omp parallel for default(shared)
#endif
  for(int i = 0; i < Local_nx; i++)
    for(int j = 0; j < Ngrid; j++)
      for(int k = 0; k < Ngrid; k++)
	density[(i*NgridL + j)*2*(NgridL/2 + 1) + k] = -1.0f;

#ifdef _OPENMP
  #pragma omp parallel for default(shared)
#endif
  for(int i=0; i<np; i++) {
    float X=P[i].x[0]*scaleBox;
    float Y=P[i].x[1]*scaleBox;
    float Z=P[i].x[2]*scaleBox;
            
    int iI=(int) floorf(X); // without floor, -1 < X < 0 is mapped to iI=0
    int J=(int) Y;          // Assumes Y,Z are positive
    int K=(int) Z;
    float D1=X-((float) iI);
    float D2=Y-((float) J);
    float D3=Z-((float) K);
    float T1=1.-D1;
    float T2=1.-D2;
    float T3=1.-D3;

    float T2W =T2*WPAR;
    float D2W =D2*WPAR;

#ifdef CHECK
    assert(Y >= 0.0f && Z >= 0.0f);
#endif
            
    // No periodic wrapup in x direction. 
    // Buffer particles are copied from adjacent nodes, instead
    if(J >= Ngrid) J=0; 
    if(K >= Ngrid) K=0;
            
    int I1=iI+1;
    int J1=J+1; if(J1 >= Ngrid) J1=0; // assumes y,z < BoxSize
    int K1=K+1; if(K1 >= Ngrid) K1=0;

    iI -= Local_x_start;
    I1 -= Local_x_start;

    if(0 <= iI && iI < Local_nx) {
      WRtPlus(density, iI, J,  K,  T3*T1*T2W);
      WRtPlus(density, iI, J,  K1, D3*T1*T2W);
      WRtPlus(density, iI, J1, K,  T3*T1*D2W);
      WRtPlus(density, iI, J1, K1, D3*T1*D2W);
    }

    if(0 <= I1 && I1 < Local_nx) {
      WRtPlus(density, I1, J,  K,  T3*D1*T2W);
      WRtPlus(density, I1, J,  K1, D3*D1*T2W);
      WRtPlus(density, I1, J1, K,  T3*D1*D2W);
      WRtPlus(density, I1, J1, K1, D3*D1*D2W);
    }
  }

  msg_printf(verbose, "CIC density assignment finished.\n");
}

// FFT density mesh and copy it to density_k
void compute_density_k(void)
{
  // FFT density(x) mesh -> density(k)
  fftwf_mpi_execute_dft_r2c(p0, (float*) fftdata, fftdata);

  // copy density(k) in fftdata to density_k

#ifdef _OPENMP
  #pragma omp parallel for default(shared)
#endif
  for(int Jl=0; Jl<Local_ny_td; Jl++) {
    for(int iI=0; iI<Ngrid; iI++) {
      for(int K=0; K<Ngrid/2+1; K++){
	size_t index= K + (NgridL/2+1)*(iI + NgridL*Jl);
	density_k[index][0]= fftdata[index][0];
	density_k[index][1]= fftdata[index][1];
      }
    }
  }
}

// Calculate one component of force mesh from precalculated density(k)
void compute_force_mesh(const int axes)
{
                                                       timer_start(force_mesh);

  fftwf_complex* const FN11= fftdata;
  fftwf_complex* const P3D= density_k;

  //k=0 zero mode force is zero
  FN11[0][0]= 0.0f;
  FN11[0][1]= 0.0f;

  const float scale=2.0*M_PI/BoxSize;

  const float f1= -1.0/pow(Ngrid, 3.0)/scale;

#ifdef _OPENMP
#pragma omp parallel for default(shared)
#endif
  for(int Jl=0; Jl<Local_ny_td; Jl++) {
    int J= Jl + Local_y_start_td;

    int J0= J <= (Ngrid/2) ? J : J - Ngrid;

    float di[3];
    di[1]= (float) J0;

    for(int iI=0; iI<Ngrid; iI++) {
      int I0= iI <= (Ngrid/2) ? iI : iI - Ngrid;
      di[0]= (float) I0;

      int KMIN= (iI==0 && J==0); // skip (0,0,0) because FN=0 for k=(0,0,0)

      for(int K=KMIN; K<Ngrid/2+1; K++){
	di[2]= (float) K;

	float f2= f1/(K*K + I0*I0 + J0*J0)*di[axes];

	size_t index= K + (NgridL/2+1)*(iI + NgridL*Jl);
	FN11[index][0]= -f2*P3D[index][1];
	FN11[index][1]=  f2*P3D[index][0];
      }
    }
  }
                                                        timer_stop(force_mesh);


                                                              timer_start(fft);

  // Force(k) -> Force(x)
  fftwf_mpi_execute_dft_c2r(p11, fftdata, (float*) fftdata);
                                                              timer_stop(fft);
}

// Does 3-linear interpolation
// particles= Values of mesh at particle positions P.x
void force_at_particle_locations(const Particle P[], const int np, 
				 const int axes, 
		                 const float fmesh[], float3 f[])
{
  const float scaleBox=((float) Ngrid)/((float) BoxSize);
     
  msg_printf(verbose, "Calculating MP...\n");

#ifdef _OPENMP
  #pragma omp parallel for default(shared)     
#endif
  for(int i=0; i<np; i++) {
    float X=P[i].x[0]*scaleBox;
    float Y=P[i].x[1]*scaleBox;
    float Z=P[i].x[2]*scaleBox;
            
    int iI= (int) floorf(X);
    int J=  (int) Y;
    int K=  (int) Z;
    float D1= X - (float) iI;
    float D2= Y - (float) J;
    float D3= Z - (float) K;
    float T1= 1.0f - D1;
    float T2= 1.0f - D2;
    float T3= 1.0f - D3;

    if(J >= Ngrid) J=0;
    if(K >= Ngrid) K=0;
            
    int I1=iI+1;
    int J1=J+1; if(J1 >= Ngrid) J1=0;
    int K1=K+1; if(K1 >= Ngrid) K1=0;

    iI -= Local_x_start;
    I1 -= Local_x_start;

    f[i][axes]= 0.0f;

    if(0 <= iI && iI < Local_nx) {
      f[i][axes] += 
	REd(fmesh, iI, J,  K )*T3*T1*T2 +
	REd(fmesh, iI, J,  K1)*D3*T1*T2 +
	REd(fmesh, iI, J1, K )*T3*T1*D2 +
	REd(fmesh, iI, J1, K1)*D3*T1*D2;
    }
    if(0 <= I1 && I1 < Local_nx) {
      f[i][axes] += 
	REd(fmesh, I1, J,  K )*T3*D1*T2 +
	REd(fmesh, I1, J,  K1)*D3*D1*T2 +
	REd(fmesh, I1, J1, K )*T3*D1*D2 +
	REd(fmesh, I1, J1, K1)*D3*D1*D2;
    }
  }
}


// send particle positions near boundary for CIC density assignment
int send_buffer_positions(Particles* const particles, BufferVec3* buf)
{
  Particle* const p= particles->p;
  const int np_local= particles->np_local;
  const int np_alloc= particles->np_allocated;

  const float eps= BoxSize/Ngrid;
  const float x_right= BoxSize/Ngrid*(Local_x_start+Local_nx) - eps;

  const int nbuf= buf->nbuf;
  float* pos= buf->vec;
  int* index= buf->index;
  float* pos_recv= buf->vrecv;

  int nsend= 0;

  msg_printf(verbose, "Exchanging buffer positions.\n");

  float shift_r= Local_x_start + Local_nx == Ngrid ? - BoxSize : 0.0f;

  // Select particles to export
  for(int i=0; i< np_local; i++) {
    if(p[i].x[0] > x_right) {
      // Send this particle to the right node
      pos[3*nsend  ]= p[i].x[0] + shift_r;
      pos[3*nsend+1]= p[i].x[1];
      pos[3*nsend+2]= p[i].x[2];
      index[nsend]= i; // record index to receive force for this particle later
      nsend++;
    }

    if(nsend >= nbuf)
      msg_abort(6100, "Error: Not enough space for sendrecv_buffer: "
		      "%d buffer particles\n", nsend);
  }
  buf->nsend= nsend;


  int copy_to= np_local;
  MPI_Status status;

  int tag= 100;

  // Get number to receive
  int nrecv;
  MPI_Sendrecv(&nsend, 1, MPI_INT, RightNode, tag, 
	       &nrecv, 1, MPI_INT, LeftNode, tag, 
	       MPI_COMM_WORLD, &status); 
  if(copy_to + nrecv > np_alloc)
    msg_abort(6110, "Error: Not enough space for buffer particles: "
	            "%d particles in Node %d", copy_to+nrecv, ThisNode);

  // Get positions
  tag++;
  MPI_Sendrecv(buf->vec, 3*nsend, MPI_FLOAT, RightNode, tag,
		 buf->vrecv,  3*nrecv, MPI_FLOAT, LeftNode, tag,
		 MPI_COMM_WORLD, &status);

  buf->nrecv= nrecv;

  // Set positions to particle data
  for(int i=0; i<nrecv; ++i) {
    p[copy_to + i].id  = -1;         // buffer particles have positions only
    p[copy_to + i].x[0]= pos_recv[3*i  ];
    p[copy_to + i].x[1]= pos_recv[3*i+1];
    p[copy_to + i].x[2]= pos_recv[3*i+2];

    if(!(BoxSize/Ngrid*Local_x_start - eps - 1.0e-5*BoxSize <= pos_recv[3*i]&&
         pos_recv[3*i] < BoxSize/Ngrid*Local_x_start + 1.0e-5*BoxSize)) {
      printf("Warning: Buffer position error? %e not in rage %e %e\n",
	     pos_recv[3*i], 
	     BoxSize/Ngrid*Local_x_start - eps, 
	     BoxSize/Ngrid*Local_x_start
	     );
    }
  }

  int nsend_global;
  MPI_Reduce(&nsend, &nsend_global, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  msg_printf(info, "%d particle positions copied for PM assignment. "
             "%d expected on average.\n", 
	     nsend_global, 
	     (Ngrid/PM_factor)*(Ngrid/PM_factor)*NNode/(PM_factor));


  return np_local + nrecv;
}

// send force information back to the particle
void add_buffer_forces(Particles* const particles, BufferVec3* buf)
{
  const int np_local= particles->np_local;
  float3* const force= particles->force;

  MPI_Status status;

  int tag= 100;
  float* const frecv= buf->vrecv;

  float3* fsend= particles->force + np_local;

  int nsend= buf->nrecv;
  int nrecv= buf->nsend;

  int ret= MPI_Sendrecv(fsend, 3*nsend, MPI_FLOAT, LeftNode, tag,
			frecv,  3*nrecv, MPI_FLOAT, RightNode, tag,
			MPI_COMM_WORLD, &status);
  assert(ret == MPI_SUCCESS);

  int* index= buf->index;
  for(int i=0; i<nrecv; ++i) {
    assert(0 <= index[i] && index[i] < np_local);
    force[index[i]][0] += frecv[3*i];
    force[index[i]][1] += frecv[3*i+1];
    force[index[i]][2] += frecv[3*i+2];
  }
}


void check_total_density(float const * const density)
{
  double sum= 0.0;

  for(int i = 0; i < Local_nx; i++)
    for(int j = 0; j < Ngrid; j++)
      for(int k = 0; k < Ngrid; k++)
	sum += density[(i*NgridL + j)*2*(NgridL/2 + 1) + k];

  double sum_global;
  MPI_Reduce(&sum, &sum_global, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);


  if(ThisNode == 0) {
    double tol= 1.0e-7*pow(PM_factor, 3)*NParticleTotal;

    if(fabs(sum_global) > tol)
      msg_abort(6040, "Error: total CIC density error is large: %le > %le\n", 
		sum_global, tol);

    msg_printf(debug, 
	      "Total CIC density OK within machine precision: %lf (< %.2lf).\n",
	       sum_global, tol);


  }

}
  
// Calculates force on particles, particles->f, using particle mesh method
void pm_calculate_forces(Particles* particles)
{  
  int np_plus_buffer= send_buffer_positions(particles, &BufPos);

                                                            timer_start(assign);
  // x_i -> density(x) = fftdata
  PtoMesh(particles->p, np_plus_buffer, (float*) fftdata);
                                                            timer_stop(assign);
                                                            timer_start(check);
  check_total_density((float*) fftdata);
                                                            timer_stop(check);


                                                            timer_start(fft);
  // density(x) -> density(k)
  compute_density_k();
                                                            timer_stop(fft);

  for(int axes=0; axes<3; axes++) {
    // density(k) -> f(x_i) [fftdata]
    compute_force_mesh(axes);

							    timer_start(pforce);
    force_at_particle_locations(particles->p, np_plus_buffer, axes,
				(float*) fftdata, particles->force);
                                                            timer_stop(pforce);
  }

                                                            timer_start(comm);
  add_buffer_forces(particles, &BufPos);
                                                            timer_stop(comm);
}

