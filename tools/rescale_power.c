//
// Reads CAMB matter power spectrum k P(k) and
// print the rescaled power spectrum
//
// Based on N-GenIC power.c by Volker Springel
//   http://www.mpa-garching.mpg.de/gadget/right.html#ICcode
//
// My cola_halo code is not renormalising the power spectrum.
// If your power spectrum is not normalised, you can use this code
// to renormalise the power spectrum file, in the same way as N-GenIC
//
// usage: rescale_power <matterpower file name> <sigma8>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>

// Code copied from power.c to have exactly the same sigma8 as N-GenIC

#define WORKSIZE 100000

static double Norm;
static int NPowerTable;
static double r_tophat;


static struct pow_table
{
  double logk, logD;
} *PowerTable;

void read_power_table_camb(const char filename[]);
double TopHatSigma2(double R);



void read_power_table_camb(const char filename[])
{
  double k, p;
  double fac= 1.0/(2.0*M_PI*M_PI);
  Norm= 1.0;

  FILE* fp= fopen(filename, "r");
  if(!fp)
    fprintf(stderr, "Error: unable to read input power spectrum: %s",
	      filename);

  NPowerTable = 0;
  do {
    if(fscanf(fp, "%lg %lg ", &k, &p) == 2)
      NPowerTable++;
    else
      break;
  } while(1);

  fprintf(stderr,
	 "Found %d pairs of values in input spectrum table\n", NPowerTable);

  PowerTable = malloc(NPowerTable * sizeof(struct pow_table));

  rewind(fp);

  int n = 0;
  do {
    if(fscanf(fp, " %lg %lg ", &k, &p) == 2) {
      PowerTable[n].logk = log10(k);
      PowerTable[n].logD = log10(fac*k*k*k*p);
      n++;
    }
    else
      break;
  } while(1);
  assert(NPowerTable == n);

  fclose(fp);
}


double PowerSpec(const double k)
{
  const double logk = log10(k);

  if(logk < PowerTable[0].logk || logk > PowerTable[NPowerTable - 1].logk)
    return 0;

  int binlow = 0;
  int binhigh = NPowerTable - 1;

  while(binhigh - binlow > 1) {
    int binmid = (binhigh + binlow) / 2;
    if(logk < PowerTable[binmid].logk)
      binhigh = binmid;
    else
      binlow = binmid;
  }

  const double dlogk = PowerTable[binhigh].logk - PowerTable[binlow].logk;
  assert(dlogk > 0.0);

  const double u = (logk - PowerTable[binlow].logk) / dlogk;

  const double logD= (1 - u) * PowerTable[binlow].logD + u * 
                     PowerTable[binhigh].logD;

  const double Delta2 = pow(10.0, logD);

  double P = Norm * Delta2 / (4.0*M_PI*k*k*k);

  return P;
}



double sigma2_int(double k, void *param)
{
  double kr, kr3, kr2, w, x;

  kr = r_tophat * k;
  kr2 = kr * kr;
  kr3 = kr2 * kr;

  if(kr < 1e-8)
    return 0;

  w = 3 * (sin(kr) / kr3 - cos(kr) / kr2);
  x = 4 * M_PI * k * k * w * w * PowerSpec(k);

  return x;
}


double TopHatSigma2(double R)
{
  double result, abserr;
  gsl_integration_workspace *workspace;
  gsl_function F;

  workspace = gsl_integration_workspace_alloc(WORKSIZE);

  F.function = &sigma2_int;

  r_tophat = R;

  gsl_integration_qag(&F, 0, 500.0 * 1 / R,
          0, 1.0e-4, WORKSIZE, GSL_INTEG_GAUSS41, workspace, &result, &abserr);
  // high precision caused error
  gsl_integration_workspace_free(workspace);

  return result;

  // note: 500/R is here chosen as (effectively) infinity integration boundary
}

int main(int argc, char* argv[])
{
  if(argc < 2) {
    printf("rescale_power <P(k) file name> <sigma8>");
    return 0;
  }

  char* filename= argv[1];
  double sigma8= atof(argv[2]);

  read_power_table_camb(filename);

  const double R8 = 8.0; // 8 /h Mpc

  double sigma2 = TopHatSigma2(R8); 
  double fac= sigma8*sigma8/sigma2*(2.0*M_PI*M_PI);

  fprintf(stderr, "Input power spectrum sigma8 %f\n", sqrt(sigma2));

  for(int i=0; i<NPowerTable; i++) {
    double k= pow(10.0, PowerTable[i].logk); 
    double pk= fac*pow(10.0, PowerTable[i].logD)/(k*k*k);

    printf("%e %e\n", k, pk);
  }
    
  return 0;
}
