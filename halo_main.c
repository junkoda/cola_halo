//
// Main file for program "halo"
//   Read Gadget snapshot and does data analysis only, 
//   FoF halo finding/density field/particle subsampling
//
// 1. Snapshot must be in directory written in `dirs.txt`
// 2. Snapshot filename must be snp_000, snp_001, ... for multiple redshifts (not a b c)
//    Modify sprintf(filebase, "%s/%s_%03d", dirname, param.snapshot_filename, iout) for other format
//    Snapshot can be distributed to multiple files snp_000.0, snp_000.1, ...
//    `snp` is specified in param.lua snapshot= "snp"
//
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <fftw3-mpi.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#include "parameters.h"
#include "lpt.h"
#include "msg.h"
#include "comm.h"
#include "fof.h"
#include "write.h"
#include "timer.h"
#include "mem.h"
#include "move.h"
#include "subsample.h"
#include "coarse_grid.h"
#include "read.h"

// TODO: update this with parameters in param
const double OmegaM      = 0.273;
const double Hubble      = 0.705;

static const int nc_factor= 3;

int mpi_init(int* p_argc, char*** p_argv);
void fft_init(int threads_ok);
void snapshot_time(const float aout, const int iout, 
       Particles const * const particles, 
       Snapshot * const snapshot, 
       const char subsample_filename[], 
       const char cgrid_filename[], const int cgrid_nc,
       void* const mem1, const size_t size1,
       const int write_longid);

int main(int argc, char* argv[])
{
  const int multi_thread= mpi_init(&argc, &argv);
  msg_init();
                                                      timer_set_category(Init);

  //
  // Initialization / Memory allocation
  //						      
  Parameters param;
  read_parameters(argc, argv, &param);
  msg_set_loglevel(param.loglevel);

  fft_init(multi_thread);
  comm_init(nc_factor*param.nc, param.nc, param.boxsize);

  Memory mem; 
  allocate_shared_memory(param.nc, nc_factor, param.np_alloc_factor, &mem); 

  const int nc= param.nc;
  ptrdiff_t local_nx, local_x_start;
  fftwf_mpi_local_size_3d(nc, nc, nc/2+1, MPI_COMM_WORLD,
        &local_nx, &local_x_start);

  Particles* particles= 
     allocate_particles(param.nc, local_nx, param.np_alloc_factor);

  Snapshot* snapshot= allocate_snapshot(param.nc, local_nx, 
          particles->np_allocated, mem.mem2, mem.size2);
  snapshot->boxsize= param.boxsize;
  snapshot->omega_m= OmegaM;
  snapshot->h= Hubble;
  snapshot->seed= 0;//***

  
  //snapshot->filename= malloc(64); // never free (memory leak)
  //strncpy(snapshot->filename, param.snapshot_filename, 64);
  snapshot->filename= param.snapshot_filename;

  fof_init(particles->np_allocated, param.nc, mem.mem1, mem.size1);
  subsample_init(param.subsample_factor, param.random_seed);

  const int this_node= comm_this_node();

  // for each directories given
  char dirname[256];

  FILE* fp_dir= 0;
  if(this_node == 0) {
    fp_dir= fopen("dirs.txt", "r");
    if(fp_dir == 0)
      msg_abort(0, "Unable to open list of directories dirs.txt\n");
  }
                                                       timer_set_category(Snp);
  while(1) {
    if(this_node == 0) {
      char* sget= fgets(dirname, 255, fp_dir);
      if(sget == NULL)
        dirname[0]= '\0';
      else {
        const int last= strlen(dirname)-1;
        if(dirname[last] == '\n')
          dirname[last]= '\0';
      }
      if(dirname[0] != '\0')
        printf("Enter directory: %s\n", dirname);
    }

    MPI_Bcast(dirname, 255, MPI_BYTE, 0, MPI_COMM_WORLD);

    if(dirname[0] == '\0')
      break;

    msg_printf(verbose, "Entering directory %s\n", dirname);

    char filebase[256];
    for(int iout=0; iout<100000; iout++) { 
      if(strlen(snapshot->filename) == 0)
        msg_abort(1, "filename snapshot is not set in the parameter file\n");

      sprintf(filebase, "%s/%s_%03d", dirname, param.snapshot_filename, iout);

      if(this_node == 0)
        printf("Read %s\n", filebase);

      int file_found= read_snapshot(filebase, snapshot, mem.mem1, mem.size1);

      if(file_found == 0) {
        if(iout == 0)
          msg_abort(2, "Unable to read snapshot %s\n", filebase);
        break;
      }
                   
      const float boxsize= snapshot->boxsize; assert(boxsize > 0.0f);

                                                              timer_start(sub);
      // subsample
      if(param.subsample_filename) {
        sprintf(filebase, "%s/%s_%03d.b", dirname, param.subsample_filename, iout);
        //sprintf(filebase, "%s%05d%c.b", param.subsample_filename, snapshot->seed, suffix);
        //write_subsample(filebase, subsample_factor, snapshot, mem1, size1);
        write_random_sabsample(filebase, snapshot, mem.mem1, mem.size1);
      }

      // coarse mesh
      if(param.cgrid_filename) {
        sprintf(filebase, "%s/%s_%03d.b", dirname, param.cgrid_filename, iout);
        //sprintf(filebase, "%s%05d%c.b", param.cgrid_filename, snapshot->seed, suffix);
        coarse_grid2(filebase, snapshot, param.cgrid_nc, mem.mem1, mem.size1);
      }

                                                               timer_stop(sub);

      const float ll= 0.2*boxsize/nc; // FOF linking length
      fof_find_halos(snapshot, ll);

      // text file of fof halo
      sprintf(filebase, "%s/fof_%03d.txt", dirname, iout);
      //sprintf(filebase, "fof%05d%c.txt", snapshot->seed, suffix);
      fof_write_halos(filebase);

                                                      timer_set_category(COLA);

      timer_print();
    }
  }
  
  if(this_node == 0)
    fclose(fp_dir);
  
  msg_printf(normal, "Processed all directories. halo done.\n");

  MPI_Finalize();
  return 0;
}

int mpi_init(int* p_argc, char*** p_argv)
{
  // MPI+OpenMP paralellization: MPI_THREAD_FUNNELED
  // supported by mpich2 1.4.1, but now by openmpi 1.2.8

#ifdef _OPENMP
  int thread_level, hybrid_parallel;
  MPI_Init_thread(p_argc, p_argv, MPI_THREAD_FUNNELED, &thread_level);
  hybrid_parallel = (thread_level >= MPI_THREAD_FUNNELED);

  int myrank; MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

  if(myrank == 0) {
    if(hybrid_parallel)
      printf("MPI + multi thread supported (MPI_THREAD_FUNNELED).\n");
    else
      printf("Warning: MPI + multi thread not supported. 1 thread per node.\n");
  }
  
  return hybrid_parallel;
#else
  MPI_Init(p_argc, p_argv);
  int myrank; MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  if(myrank == 0)
    printf("MPI only without OpenMP\n");
  return 0;
#endif

}

void fft_init(int threads_ok)
{
  // Initialize FFTW3

#ifdef _OPENMP
  if(threads_ok)
    threads_ok= fftwf_init_threads();
  if(!threads_ok)
    msg_printf(warn, "Multi-thread FFTW not supported.\n");
#endif
    
  fftwf_mpi_init();

#ifdef _OPENMP
  if(threads_ok) {
    int nthreads= omp_get_max_threads();
    fftwf_plan_with_nthreads(nthreads);
    msg_printf(info, "Multi-threaded FFTW: %d threads\n", nthreads);
  }
#endif

}

    
