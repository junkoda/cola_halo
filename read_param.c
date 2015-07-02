//
// This code reads simple parameter files
// -- comment
// variable_number = 123
// variable_string = "string"
// variable_number_array = {1,2,3} 

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <assert.h>
#include <stdbool.h>
#include <mpi.h>
#include "parameters.h"

#include "msg.h"

static int myrank_;
//static const buf_len= 1024; // ToDo: malloc buf

int read_parameter_file(const char filename[], Parameters* const param);
static void bcast_string(char** string, int* len);
static void bcast_array_double(double** parray, int* len);

// Utilities
static void remove_comments(char* str)
{
  for(char* p=str; *p != '\0'; p++) {
    if(*p == '-' && *(p+1) == '-') {
      *p='\0';
      return;
    }
    else if(*p == '\n') {
      *p='\0';
      return;
    }
  }
}

static char* skip_spaces(char* p)
{
  while(*p == ' ' || *p == '\t') p++;
  return p;
}

static char* get_name(char* p, char* name)
{
  p= skip_spaces(p);

  while(*p != ' ' && *p != '=' && p != '\0')
    *name++ = *p++;

  *name='\0';
  
  return p+1;
}

static char* get_operator(char* p)
{
  p= skip_spaces(p);

  while(*p == '=')
    p++;

  return p;
}

static int get_int(char* p)
{
  p= skip_spaces(p);

  if(isdigit(*p)) 
    return atoi(p);

  return 0;
}

static double get_double(char* p)
{
  p= skip_spaces(p);

  if(isdigit(*p)) 
    return atof(p);

  return 0.0;
}

static char* get_string(char* p, int* len)
{
  // warning: no protection for string (filename) over 1024 characters
  char const * const original_p= p;
  char buf_str[1024];
  char* dest= buf_str;
  p= skip_spaces(p);

  if(*p == '"') {
    // Copy the string in "" to buf_str.
    p++;
    while(*p != '"' && *p != '\0' && *p != '\n')
      *dest++ = *p++;
  }
  else {
    msg_abort(1003, "Error: string value not starting from \": %s\n",
	      original_p);
  }

  *dest= '\0';

  *len= strlen(buf_str) + 1; // +1 for the terminating null '\0'
  char* val= malloc(*len); assert(val);
  strcpy(val, buf_str);
  
  return val;
}

static int get_bool(char* p)
{
  p= skip_spaces(p);

  if(*p == 't' || *p == 'T' || *p == '1')
    return 1;
  else if(*p == 'f' || *p == 'F' || *p == '0')
    return 0;

  msg_abort(1002, "Error: unknown bool value %s\n", p);
  abort();

  return 0;
}


static double* get_numarray(char* p, int* len)
{
  char const * const original_p= p;
  int nalloc= 10;
  double* array= (double*) malloc(sizeof(double)*nalloc);
  
  while(*p == ' ' || *p == '\t') p++;


  int n=0;
  if(*p == '{') {
    p++;

    do {
      if(n >= nalloc) {
	// increase the size of the array if necessary
	nalloc *= 2;
	array= (double*) realloc(array, sizeof(double)*nalloc);
      }

      char numbuf[64]; char* q= numbuf;

      p= skip_spaces(p);
	
      while(*p != ',' && *p != '}' && *p != ' ' && *p != '\0')
	*q++ = *p++;
      *q= '\0';

      if(!isdigit(numbuf[0]) && numbuf[0] != '-' && numbuf[0] != '+' &&
	 numbuf[0] != '.')
	msg_abort(1004, "Error: unable to parse array: %s\n", original_p);
		  
      array[n++]= atof(numbuf);

      if(*p == ',')
	p++;
      else if(*p == '}') {
	*len= n;
	return array;
      }
    } while(*p != '\0');
  }

  msg_abort(1003, "Error in parsing the array: %s\n", original_p);
  
  free(array);
  return 0;
}

int read_parameters(const int argc, char * argv[], Parameters* const param)
{
  if(argc < 2)
    msg_abort(1, "Error: Parameter file not specified. cola_code param.lua\n");

  char const * const filename= argv[argc-1];

  //int myrank;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank_);
  if(myrank_ == 0) {
    int ret= read_parameter_file(filename, param);
    if(ret != 0)
      msg_abort(1001, "Error: Unable to read parameter file: %s\n", filename);
  }

  // Share parameters with other nodes
  MPI_Bcast(param, sizeof(Parameters), MPI_BYTE, 0, MPI_COMM_WORLD);

  bcast_string(&param->power_spectrum_filename, 
	       &param->strlen_power_spectrum_filename);
  bcast_string(&param->fof_filename,        &param->strlen_fof_filename);
  bcast_string(&param->snapshot_filename,   &param->strlen_snapshot_filename);
  bcast_string(&param->subsample_filename,  &param->strlen_subsample_filename);
  bcast_string(&param->cgrid_filename,      &param->strlen_cgrid_filename);
  bcast_string(&param->init_filename,       &param->strlen_init_filename);

  bcast_array_double(&param->zout,          &param->n_zout);

  return 0;
}

int read_parameter_file(const char filename[], Parameters* const param)
{
  FILE* fp= fopen(filename, "r");
  if(fp == 0) return -1;
  
  char buf[1024];
  char name[64];

  param->power_spectrum_filename   = 0;
  param->snapshot_filename         = 0;
  param->fof_filename              = 0;
  param->cgrid_filename            = 0;
  param->subsample_filename        = 0;
  param->init_filename             = 0;

  param->strlen_fof_filename       = 0;
  param->strlen_snapshot_filename  = 0;
  param->strlen_cgrid_filename     = 0;
  param->strlen_subsample_filename = 0;
  param->strlen_init_filename      = 0;


  param->omega_m = -1.0;
  param->sigma8  = -1.0;
  param->h       = -1.0;
  param->boxsize = -1.0;
  
  param->a_final =  1.0;

  param->pm_nc_factor= 3;
  param->random_seed=  1;
  param->nrealization= 1;
  param->loglevel=     0;
  param->ntimestep=    0;

  param->zout=         0;
  param->n_zout=       0;

  param->np_alloc_factor= 1.25;
  param->fof_linking_factor= 0.2;
  param->subsample_factor= 0.0;
  param->cgrid_nc= 0;
  param->write_longid= 0; 

  while(fgets(buf, 1023, fp)) {
    remove_comments(buf);

    char* p= get_name(buf, name);

    if(name[0] == '\0')
      continue;

    p= get_operator(p);
    
    if(strcmp(name, "nc") == 0)
      param->nc= (int) get_int(p);
    else if(strcmp(name, "boxsize") == 0)
      param->boxsize= get_double(p);
    else if(strcmp(name, "a_final") == 0)
      param->a_final= get_double(p);
    else if(strcmp(name, "ntimestep") == 0)
      param->ntimestep= get_int(p);
    else if(strcmp(name, "output_redshifts") == 0)
      param->zout= get_numarray(p, &param->n_zout);
    else if(strcmp(name, "random_seed") == 0)
      param->random_seed= get_int(p);
    else if(strcmp(name, "omega_m") == 0)
      param->omega_m= get_double(p);
    else if(strcmp(name, "nrealization") == 0)
      param->nrealization= get_int(p);
    else if(strcmp(name, "h") == 0)
      param->h= get_double(p);
    else if(strcmp(name, "sigma8") == 0)
      param->sigma8= get_double(p);
    else if(strcmp(name, "pm_nc_factor") == 0)
      param->pm_nc_factor= get_int(p);
    else if(strcmp(name, "np_alloc_factor") == 0)
      param->np_alloc_factor= get_double(p);
    else if(strcmp(name, "loglevel") == 0)
      param->loglevel= get_int(p);
    else if(strcmp(name, "powerspectrum") == 0)
      param->power_spectrum_filename=
	get_string(p, &param->strlen_power_spectrum_filename);
    else if(strcmp(name, "fof") == 0)
      param->fof_filename= get_string(p, &param->strlen_fof_filename);
    else if(strcmp(name, "linking_factor") == 0)
      param->fof_linking_factor= get_double(p);
    else if(strcmp(name,  "snapshot") == 0)
      param->snapshot_filename=	get_string(p, &param->strlen_snapshot_filename);
    else if(strcmp(name, "subsample") == 0)
      param->subsample_filename=
	get_string(p, &param->strlen_subsample_filename);
    else if(strcmp(name, "subsample_factor") == 0)
      param->subsample_factor= get_double(p);
    else if(strcmp(name, "coarse_grid") == 0)
      param->cgrid_filename= get_string(p, &param->strlen_cgrid_filename);
    else if(strcmp(name, "coarse_grid_nc") == 0)
      param->cgrid_nc= get_int(p);
    else if(strcmp(name, "initial") == 0)
        param->init_filename= get_string(p, &param->strlen_init_filename);
    else if(strcmp(name, "write_longid") == 0)
      param->write_longid= get_bool(p);
    else
      msg_abort(1007,
		 "Error: unknown variable name in parameter file= %s\n", name);
  }

  return 0;
}

void bcast_string(char** pstring, int* len)
{
  const int ret1= MPI_Bcast(len, 1, MPI_INT, 0, MPI_COMM_WORLD);
    assert(ret1 == MPI_SUCCESS);

  const int n= *len;

  if(n == 0) {
    *pstring= 0;
    return;
  }

  if(myrank_ != 0) {
    *pstring= malloc(sizeof(char)*n);
  }
    assert(*pstring);

  const int ret2= MPI_Bcast(*pstring, n, MPI_CHAR, 0, MPI_COMM_WORLD);
    assert(ret2 == MPI_SUCCESS);
}

void bcast_array_double(double** parray, int* len)
{
  const int ret1= MPI_Bcast(len, 1, MPI_INT, 0, MPI_COMM_WORLD);
    assert(ret1 == MPI_SUCCESS);

  const int n= *len;

  if(n == 0) {
    *parray= 0;
    return;
  }

  if(myrank_ != 0) {
    *parray= malloc(sizeof(double)*n);
  }
    assert(*parray);

  const int ret2= MPI_Bcast(*parray, n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    assert(ret2 == MPI_SUCCESS);
}
