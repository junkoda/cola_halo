//
// Computes random Gaussian 2LPT displacement on regular particle grid
//
// Based on the code by Roman Scoccimaro, Sebastian Pueblas, Marc Manera et al
// http://cosmo.nyu.edu/roman/2LPT/
//
// link with -lfftw3f_mpi -lfftw3f -lm for single precision FFTW3
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include <fftw3-mpi.h>
#include <gsl/gsl_rng.h>

#include "msg.h"
#include "power.h"
#include "particle.h"

#ifndef M_PI
#define M_PI 3.1415926535897932384626433832795
#endif

static fftwf_plan Inverse_plan[6], Forward_plan;
static fftwf_plan Disp_plan[3], Disp2_plan[3];

static int Nmesh, Nsample;
static int Local_nx, Local_x_start;

static double Omega, OmegaLambda;

static unsigned int *seedtable;

fftwf_complex *(cdisp[3]), *(cdisp2[3]); // ZA and 2nd order displacements
float         *(disp[3]), *(disp2[3]);
fftwf_complex *(cdigrad[6]);
float         *(digrad[6]);


double F_Omega(const double a);
double F2_Omega(const double a);

// Setup variables for 2LPT initial condition
void lpt_init(const int nc, const void* mem, const size_t size)
{
  // nc: number of mesh per dimension

  ptrdiff_t local_nx, local_x_start;
  ptrdiff_t total_size=
    fftwf_mpi_local_size_3d(nc, nc, nc/2+1, MPI_COMM_WORLD,
			    &local_nx, &local_x_start);
  
  Local_nx= local_nx;
  Local_x_start= local_x_start;

  //
  // Allocate memory
  //

  if(mem == 0) {
    // allocate memory here
    size_t bytes= sizeof(fftwf_complex)*total_size;
    int allocation_failed= 0;

    // 1&2 displacement
    for(int axes=0; axes < 3; axes++) {
      cdisp[axes]= fftwf_alloc_complex(total_size);
      disp[axes] = (float*) cdisp[axes];
      
      cdisp2[axes]= fftwf_alloc_complex(total_size);
      disp2[axes] = (float*) cdisp2[axes];
      bytes += 2*sizeof(fftwf_complex)*total_size;
      
      allocation_failed = allocation_failed ||
	(cdisp[axes] == 0) || (cdisp2[axes] == 0);
    } 
    
    // 2LPT
    for(int i=0; i<6; i++) {
      cdigrad[i] = (fftwf_complex *) fftwf_alloc_complex(total_size);
      digrad[i] = (float*) cdigrad[i];
      
      bytes += sizeof(fftwf_complex)*total_size;
      allocation_failed = allocation_failed || (digrad[i] == 0);
    } 
    
   if(allocation_failed)
     msg_abort(2003, "Error: Failed to allocate memory for 2LPT."
	       "Tried to allocate %d Mbytes\n", (int)(bytes/(1024*1024)));
   
   msg_printf(info, "%d Mbytes allocated for LPT\n", (int)(bytes/(1024*1024)));
  }
  else {
    size_t bytes= 0;
    fftwf_complex* p= (fftwf_complex*) mem;
    
    for(int axes=0; axes<3; axes++) {
      cdisp[axes]= p;
      disp[axes]= (float*) p;
      bytes += sizeof(fftwf_complex)*total_size*2;

      p += total_size;
    }
    for(int i=0; i<6; i++) {
      cdigrad[i]= p;
      digrad[i]= (float*) p;
      bytes += sizeof(fftwf_complex)*total_size;
      p += total_size;
    }
    assert(bytes <= size);
  }

  //
  // FFTW3 plans
  //
  for(int i=0; i<6; ++i)
    Inverse_plan[i]=
      fftwf_mpi_plan_dft_c2r_3d(nc, nc, nc, cdigrad[i], digrad[i],
				MPI_COMM_WORLD, FFTW_ESTIMATE);

  Forward_plan=
    fftwf_mpi_plan_dft_r2c_3d(nc, nc, nc, digrad[3], cdigrad[3],
			      MPI_COMM_WORLD, FFTW_ESTIMATE);

  for(int i=0; i<3; ++i) {
    Disp_plan[i]=
      fftwf_mpi_plan_dft_c2r_3d(nc, nc, nc, cdisp[i], disp[i],
				MPI_COMM_WORLD, FFTW_ESTIMATE);
    Disp2_plan[i]=
      fftwf_mpi_plan_dft_c2r_3d(nc, nc, nc, cdisp2[i], disp2[i],
				MPI_COMM_WORLD, FFTW_ESTIMATE);
  }

  // FFTW_MPI_TRANSPOSED_IN/FFTW_MPI_TRANSPOSED_OUT would be faster
  // FFTW_MEASURE is probably better for multiple realization  

  // misc data
  Nmesh= nc;
  Nsample= nc;
  seedtable = malloc(Nmesh * Nmesh * sizeof(unsigned int)); assert(seedtable);
}

int lpt_set_displacement(const double InitTime, const double omega_m, const int Seed, const double Box, Particles* particles)
{
  msg_printf(verbose, "Computing LPT displacement fields...\n");
  msg_printf(info, "Random Seed = %d\n", Seed);

  //
  // Setting constant parameters
  //
  Omega= omega_m;
  OmegaLambda= 1.0 - omega_m;

  


  static const double Hubble= 3.085678e24/1.0e5*3.2407789e-18; 
  // Hubble = 100 [km/s(/Mpc/h)]
  const double hubble_a= Hubble * sqrt(Omega / pow(InitTime, 3) + (1 - Omega - OmegaLambda) / pow(InitTime, 2) + OmegaLambda);

  double vel_prefac = InitTime * hubble_a * F_Omega(InitTime);
  double vel_prefac2 = InitTime * hubble_a * F2_Omega(InitTime);

  const double fac = pow(2*M_PI/Box, 1.5);

  //
  // Setup random seeds 
  //  complicated way for backward compatibility?
  //
  gsl_rng* random_generator = gsl_rng_alloc(gsl_rng_ranlxd1);
  gsl_rng_set(random_generator, Seed);

  assert(seedtable);

  for(int i=0; i<Nmesh/2; i++) {
    for(int j=0; j<i; j++)
      seedtable[i * Nmesh + j] = 0x7fffffff * gsl_rng_uniform(random_generator);

    for(int j=0; j<i+1; j++)
      seedtable[j * Nmesh + i] = 0x7fffffff * gsl_rng_uniform(random_generator);

    for(int j=0; j<i; j++)
      seedtable[(Nmesh - 1 - i) * Nmesh + j] = 
	0x7fffffff * gsl_rng_uniform(random_generator);

    for(int j=0; j<i+1; j++)
      seedtable[(Nmesh - 1 - j) * Nmesh + i] = 
	0x7fffffff * gsl_rng_uniform(random_generator);

    for(int j=0; j<i; j++)
      seedtable[i * Nmesh + (Nmesh - 1 - j)] = 
	0x7fffffff * gsl_rng_uniform(random_generator);

    for(int j=0; j<i+1; j++)
      seedtable[j * Nmesh + (Nmesh - 1 - i)] = 
	0x7fffffff * gsl_rng_uniform(random_generator);

    for(int j=0; j<i; j++)
      seedtable[(Nmesh - 1 - i) * Nmesh + (Nmesh - 1 - j)] = 
	0x7fffffff * gsl_rng_uniform(random_generator);

    for(int j=0; j<i+1; j++)
      seedtable[(Nmesh - 1 - j) * Nmesh + (Nmesh - 1 - i)] = 
	0x7fffffff * gsl_rng_uniform(random_generator);
  }

  // clean the array
  for(int i = 0; i < Local_nx; i++)
    for(int j = 0; j < Nmesh; j++)
      for(int k = 0; k <= Nmesh / 2; k++)
	for(int axes = 0; axes < 3; axes++) {
	  cdisp[axes][(i * Nmesh + j) * (Nmesh / 2 + 1) + k][0] = 0.0f;
	  cdisp[axes][(i * Nmesh + j) * (Nmesh / 2 + 1) + k][1] = 0.0f;
	}

  double kvec[3];
  for(int i = 0; i < Nmesh; i++) {
    int ii = Nmesh - i;
    if(ii == Nmesh)
      ii = 0;
    if((i >= Local_x_start && i < (Local_x_start + Local_nx)) ||
       (ii >= Local_x_start && ii < (Local_x_start + Local_nx))) {
      for(int j = 0; j < Nmesh; j++) {
	gsl_rng_set(random_generator, seedtable[i * Nmesh + j]);
	
	for(int k = 0; k < Nmesh / 2; k++) {
	  double phase = gsl_rng_uniform(random_generator) * 2 * M_PI;
	  double ampl;
	  do
	    ampl = gsl_rng_uniform(random_generator);
	  while(ampl == 0.0);
	  
	  if(i == Nmesh / 2 || j == Nmesh / 2 || k == Nmesh / 2)
	    continue;
	  if(i == 0 && j == 0 && k == 0)
	    continue;
	  
	  if(i < Nmesh / 2)
	    kvec[0] = i * 2 * M_PI / Box;
	  else
	    kvec[0] = -(Nmesh - i) * 2 * M_PI / Box;
	  
	  if(j < Nmesh / 2)
	    kvec[1] = j * 2 * M_PI / Box;
	  else
	    kvec[1] = -(Nmesh - j) * 2 * M_PI / Box;
	  
	  if(k < Nmesh / 2)
	    kvec[2] = k * 2 * M_PI / Box;
	  else
	    kvec[2] = -(Nmesh - k) * 2 * M_PI / Box;
	  
	  double kmag2 = kvec[0]*kvec[0] + kvec[1]*kvec[1] + kvec[2]*kvec[2];
	  double kmag = sqrt(kmag2);
	  
#ifdef SPHEREMODE
	  // select a sphere in k-space
	  if(kmag * Box / (2 * M_PI) > Nsample / 2)
	    continue;
#else
	  if(fabs(kvec[0]) * Box / (2 * M_PI) > Nsample / 2)
	    continue;
	  if(fabs(kvec[1]) * Box / (2 * M_PI) > Nsample / 2)
	    continue;
	  if(fabs(kvec[2]) * Box / (2 * M_PI) > Nsample / 2)
	    continue;
#endif
		      
	  double p_of_k = PowerSpec(kmag);
		      
	  p_of_k *= -log(ampl);
		      
	  double delta = fac * sqrt(p_of_k); // / Dplus;	
	  // Displacement is extrapolated to a=1
		      
	  if(k > 0) {
	    if(i >= Local_x_start && i < (Local_x_start + Local_nx))
	      for(int axes = 0; axes < 3; axes++) {
	        cdisp[axes][((i - Local_x_start) * Nmesh + j) * (Nmesh / 2 + 1)
			    + k][0] =
		                      -kvec[axes] / kmag2 * delta * sin(phase);
		cdisp[axes][((i - Local_x_start) * Nmesh + j) * (Nmesh / 2 + 1)
			    + k][1] =
		                       kvec[axes] / kmag2 * delta * cos(phase);
	      }
	  }
	  else { // k=0 plane needs special treatment
	    if(i == 0) {
	      if(j >= Nmesh / 2)
		continue;
	      else {
		if(i >= Local_x_start && i < (Local_x_start + Local_nx)) {
		  int jj = Nmesh - j; // note: j!=0 surely holds at this point
		  
		  for(int axes = 0; axes < 3; axes++) {
		    cdisp[axes][((i - Local_x_start) * Nmesh + j) * 
				(Nmesh / 2 + 1) + k][0] =
		                      -kvec[axes] / kmag2 * delta * sin(phase);
		    cdisp[axes][((i - Local_x_start) * Nmesh + j) * 
				(Nmesh / 2 + 1) + k][1] =
		                       kvec[axes] / kmag2 * delta * cos(phase);
					  
		    cdisp[axes][((i - Local_x_start) * Nmesh + jj) * 
				(Nmesh / 2 + 1) + k][0] =
			      	    -kvec[axes] / kmag2 * delta * sin(phase);
		    cdisp[axes][((i - Local_x_start) * Nmesh + jj) * 
				(Nmesh / 2 + 1) + k][1] =
				    -kvec[axes] / kmag2 * delta * cos(phase);
		  }
		}
	      }
	    }
	    else { // here comes i!=0 : conjugate can be on other processor!
	      if(i >= Nmesh / 2)
		continue;
	      else {
		ii = Nmesh - i;
		if(ii == Nmesh)
		  ii = 0;
		int jj = Nmesh - j;
		if(jj == Nmesh)
		  jj = 0;
		
		if(i >= Local_x_start && i < (Local_x_start + Local_nx))
		  for(int axes = 0; axes < 3; axes++) {
		    cdisp[axes][((i - Local_x_start) * Nmesh + j) * 
				(Nmesh / 2 + 1) + k][0] =
		                      -kvec[axes] / kmag2 * delta * sin(phase);
		    cdisp[axes][((i - Local_x_start) * Nmesh + j) * 
				(Nmesh / 2 + 1) + k][1] =
				       kvec[axes] / kmag2 * delta * cos(phase);
		  }
		
		if(ii >= Local_x_start && ii < (Local_x_start + Local_nx))
		  for(int axes = 0; axes < 3; axes++) {
		    cdisp[axes][((ii - Local_x_start) * Nmesh + jj) * 
				(Nmesh / 2 + 1) + k][0] = 
		                      -kvec[axes] / kmag2 * delta * sin(phase);
		    cdisp[axes][((ii - Local_x_start) * Nmesh + jj) * 
				(Nmesh / 2 + 1) + k][1] = 
		                      -kvec[axes] / kmag2 * delta * cos(phase);
		  }
	      }
	    }
	  }
	}
      }
    }
  }

  
  //
  // 2nd order LPT
  //
      
  for(int i = 0; i < Local_nx; i++) {
    for(int j = 0; j < Nmesh; j++) {
      for(int k = 0; k <= Nmesh / 2; k++) {
	int coord = (i * Nmesh + j) * (Nmesh / 2 + 1) + k;
	if((i + Local_x_start) < Nmesh / 2)
	  kvec[0] = (i + Local_x_start) * 2 * M_PI / Box;
	else
	  kvec[0] = -(Nmesh - (i + Local_x_start)) * 2 * M_PI / Box;
	      
	if(j < Nmesh / 2)
	  kvec[1] = j * 2 * M_PI / Box;
	else
	  kvec[1] = -(Nmesh - j) * 2 * M_PI / Box;
	      
	if(k < Nmesh / 2)
	  kvec[2] = k * 2 * M_PI / Box;
	else
	  kvec[2] = -(Nmesh - k) * 2 * M_PI / Box;
	      
	// Derivatives of ZA displacement
	// d(dis_i)/d(q_j)  -> sqrt(-1) k_j dis_i
	cdigrad[0][coord][0] = -cdisp[0][coord][1] * kvec[0]; // disp0,0
	cdigrad[0][coord][1] =  cdisp[0][coord][0] * kvec[0];

	cdigrad[1][coord][0] = -cdisp[0][coord][1] * kvec[1]; // disp0,1
	cdigrad[1][coord][1] =  cdisp[0][coord][0] * kvec[1];

	cdigrad[2][coord][0] = -cdisp[0][coord][1] * kvec[2]; // disp0,2
	cdigrad[2][coord][1] =  cdisp[0][coord][0] * kvec[2];
	      
	cdigrad[3][coord][0] = -cdisp[1][coord][1] * kvec[1]; // disp1,1
	cdigrad[3][coord][1] =  cdisp[1][coord][0] * kvec[1];

	cdigrad[4][coord][0] = -cdisp[1][coord][1] * kvec[2]; // disp1,2
	cdigrad[4][coord][1] =  cdisp[1][coord][0] * kvec[2];

	cdigrad[5][coord][0] = -cdisp[2][coord][1] * kvec[2]; // disp2,2
	cdigrad[5][coord][1] =  cdisp[2][coord][0] * kvec[2];
      }
    }
  }

  msg_printf(verbose, "Fourier transforming displacement gradient...");

  for(int i = 0; i < 6; i++)
    fftwf_mpi_execute_dft_c2r(Inverse_plan[i], cdigrad[i], digrad[i]);

  msg_printf(verbose, "Done.\n");

  // Compute second order source and store it in digrad[3]

  for(int i = 0; i < Local_nx; i++) {
    for(int j = 0; j < Nmesh; j++) {
      for(int k = 0; k < Nmesh; k++) {
	int coord = (i * Nmesh + j) * (2 * (Nmesh / 2 + 1)) + k;

	digrad[3][coord] =
	  digrad[0][coord]*(digrad[3][coord] + digrad[5][coord]) +
	  digrad[3][coord]*digrad[5][coord] - 
          digrad[1][coord]*digrad[1][coord] -
          digrad[2][coord]*digrad[2][coord] -
	  digrad[4][coord]*digrad[4][coord];
      }
    }
  }

  msg_printf(verbose, "Fourier transforming second order source...\n");
  
  //fft_execute(Forward_plan);
  fftwf_mpi_execute_dft_r2c(Forward_plan, digrad[3], cdigrad[3]);

  // The memory allocated for cdigrad[0], [1], and [2] will be used for 
  // 2nd order displacements
  // cdigrad[3] has 2nd order displacement source.

  for(int axes = 0; axes < 3; axes++) {
    cdisp2[axes] = cdigrad[axes]; 
    disp2[axes] = (float*) cdisp2[axes];
  }

  // Solve Poisson eq. and calculate 2nd order displacements

  for(int i = 0; i < Local_nx; i++) {
    for(int j = 0; j < Nmesh; j++) {
      for(int k = 0; k <= Nmesh / 2; k++) {
	int coord = (i * Nmesh + j) * (Nmesh / 2 + 1) + k;
	if((i + Local_x_start) < Nmesh / 2)
	  kvec[0] = (i + Local_x_start) * 2 * M_PI / Box;
	else
	  kvec[0] = -(Nmesh - (i + Local_x_start)) * 2 * M_PI / Box;
	
	if(j < Nmesh / 2)
	  kvec[1] = j * 2 * M_PI / Box;
	else
	  kvec[1] = -(Nmesh - j) * 2 * M_PI / Box;
	
	if(k < Nmesh / 2)
	  kvec[2] = k * 2 * M_PI / Box;
	else
	  kvec[2] = -(Nmesh - k) * 2 * M_PI / Box;
	
	double kmag2= kvec[0] * kvec[0] + kvec[1] * kvec[1] + kvec[2] * kvec[2];
#ifdef CORRECT_CIC
	// calculate smooth factor for deconvolution of CIC interpolation
	fx = fy = fz = 1;
	if(kvec[0] != 0) {
	  fx = (kvec[0] * Box / 2) / Nmesh;
	  fx = sin(fx) / fx;
	}
	if(kvec[1] != 0) {
	  fy = (kvec[1] * Box / 2) / Nmesh;
	  fy = sin(fy) / fy;
	}
	if(kvec[2] != 0) {
	  fz = (kvec[2] * Box / 2) / Nmesh;
	  fz = sin(fz) / fz;
	}
	ff = 1 / (fx * fy * fz);
	smth = ff * ff;
#endif
	    
	// cdisp2 = source * k / (sqrt(-1) k^2)
	for(int axes = 0; axes < 3; axes++) {
	  if(kmag2 > 0.0) {
	    cdisp2[axes][coord][0] =  cdigrad[3][coord][1] * kvec[axes] / kmag2;
	    cdisp2[axes][coord][1] = -cdigrad[3][coord][0] * kvec[axes] / kmag2;
	  }
	  else 
	    cdisp2[axes][coord][0] = cdisp2[axes][coord][1] = 0.0;
#ifdef CORRECT_CIC
	  cdisp[axes][coord][0]  *= smth;  cdisp[axes][coord][1]  *= smth;
	  cdisp2[axes][coord][0] *= smth;  cdisp2[axes][coord][1] *= smth;
#endif
	}
      }
    }
  }
      
  // Now, both cdisp, and cdisp2 have the ZA and 2nd order displacements

  for(int axes = 0; axes < 3; axes++) {  
    msg_printf(verbose, "Fourier transforming displacements, axis %d.\n",axes);

    fftwf_mpi_execute_dft_c2r(Disp_plan[axes], cdisp[axes], disp[axes]);
    fftwf_mpi_execute_dft_c2r(Disp2_plan[axes], cdisp2[axes], disp2[axes]);
  }

  msg_printf(verbose, "Setting particle grid and displacements\n");

  float x[3];
  const float dx= Box/Nmesh;
  const float Dplus = 1.0/GrowthFactor(InitTime, 1.0);
  Particle* p= particles->p;
  float maxdisp= 0.0f;



  double nmesh3_inv= 1.0/pow((double)Nmesh, 3.0);
  long long id= (long long) Local_x_start*Nmesh*Nmesh + 1;

  // small change to make it consistent with cola
  // use consistent growth factor... **

  const double a= InitTime;
  const double omega=Omega/(Omega + (1.0 - Omega)*a*a*a);
  const double D2= Dplus*Dplus*pow(omega/Omega, -1.0/143.0);
  const double D20= pow(Omega, -1.0/143.0);

  msg_printf(verbose, "Inital growth factor (a=%g), D= %g D2= %g\n", 
	     InitTime, Dplus, D2);
  msg_printf(debug, "2LPT vel_prefac (a=%5.3f) %e %e\n", InitTime, Dplus*vel_prefac, Dplus*Dplus*vel_prefac2);



  for(int i=0; i<Local_nx; i++) {
    x[0]= (Local_x_start + i + 0.5f)*dx;
    for(int j=0; j<Nmesh; j++) {
      x[1]= (j + 0.5f)*dx;
      for(int k=0; k<Nmesh; k++) {
	x[2]= (k + 0.5f)*dx;

	for(int axes = 0; axes < 3; axes++) {
	  float dis= disp[axes][(i * Nmesh + j) * (2 * (Nmesh / 2 + 1)) + k];
	  float dis2= nmesh3_inv*
	            disp2[axes][(i * Nmesh + j) * (2 * (Nmesh / 2 + 1)) + k];

	  //p->x[axes]= x[axes] + dis*Dplus - 3.0/7.0*dis2*Dplus*Dplus;
	  p->x[axes]= x[axes] + dis*Dplus - 3.0/7.0*D20*dis2*D2;
	  p->dx1[axes]= dis;                 // 1LPT extrapolated to a=1
	  p->dx2[axes]= -3.0/7.0*D20*dis2;   // 2LPT extrapolated to a=1
	  p->v[axes]= 0.0f;                  // velocity in comoving 2LPT
	  
	  //2LPT velocity
	  //P[n].Vel[axes] = dis * vel_prefac - 3.0/7.0 * dis2 * vel_prefac2;
	  //dis & dis2 has growth factor
	  float dis_mag= (float) fabs(dis*Dplus - 3.0/7.0*D20*dis2*D2);
	  if(dis_mag > maxdisp)
	    maxdisp= dis_mag;
	}
	p->id= id++;
	p++;
      }
    }
  }

  gsl_rng_free(random_generator);

  float max_disp_glob;
  MPI_Reduce(&maxdisp, &max_disp_glob, 1, MPI_FLOAT, MPI_MAX, 0, 
	     MPI_COMM_WORLD);

  msg_printf(verbose, 
	     "Maximum displacement: %f, in units of the part-spacing= %f\n",
	     max_disp_glob, max_disp_glob / dx);

  particles->np_local= Local_nx*Nmesh*Nmesh;
  particles->a_x= InitTime;
  particles->a_v= 0.0;

  return 0;
}


void lpt_finalize(void)
{
  for(int i=0; i<3; i++) {
    fftwf_free(cdisp[i]);
    fftwf_free(cdisp2[i]);
  }
  for(int i=0; i<6; i++)
    fftwf_free(cdigrad[i]);

  fftwf_mpi_cleanup();
}


double F_Omega(const double a)
{
  double omega_a = Omega / 
    (Omega + a*(1.0 - Omega - OmegaLambda) + a*a*a*OmegaLambda);

  return pow(omega_a, 0.6);
}


double F2_Omega(const double a)
{
  double omega_a = Omega / 
    (Omega + a*(1 - Omega - OmegaLambda) + a*a*a*OmegaLambda);

  return 2.0 * pow(omega_a, 4.0/7.0);
}

int lpt_get_local_nx(void)
{
  return Local_nx;
}
