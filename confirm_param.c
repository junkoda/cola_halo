//
// Check necessary parameters are initialized and print the parameters
//
#include <stdio.h>
#include "parameters.h"
#include "msg.h"

void confirm_parameters(Parameters* const param, const char filename[])
{
  msg_printf(normal, "parameter file %s read.\n", filename);

  if(param->nc <= 0)
    msg_abort(11001, "nc is not given is the parameter file: %d\n", param->nc);
  else
    msg_printf(normal, "nc = %d\n", param->nc);

  msg_printf(normal, "pm_nc_factor = %d\n", param->pm_nc_factor);

  if(param->boxsize <= 0.0)
    msg_abort(11002, "boxsize is not given in the parameter file: %f\n",
	      param->boxsize);
  else
    msg_printf(normal, "boxsize = %lg\n", param->boxsize);

  msg_printf(normal, "a_final = %lg\n", param->a_final);

  if(param->ntimestep <= 0)
    msg_abort(11003, "ntimestep is not given in the parameter file: %d\n",
	      param->ntimestep);
  else
    msg_printf(normal, "ntimestep = %d\n", param->ntimestep);

  if(param->n_zout == 0)
    msg_abort(11004, "output_redshifts are not given in the parameter file.\n");
  else {
    msg_printf(normal, "output_redshift = {%lg", param->zout[0]);
    for(int i=1; i<param->n_zout; i++)
      msg_printf(normal, ", %lg", param->zout[i]);
    msg_printf(normal, "}\n");
  }

  msg_printf(normal, "random_seed = %d\n", param->random_seed);
  msg_printf(normal, "nrealization = %d\n", param->nrealization);

  if(param->omega_m <= 0.0)
    msg_abort(11005, "omega_m is not given in the parameter file.");
  else
    msg_printf(normal, "omega_m = %lg\n", param->omega_m);

  if(param->h <= 0.0)
    msg_abort(11006, "h is not given in the parameter file");
  else
    msg_printf(normal, "h = %lg\n", param->h);

  if(param->sigma8 <= 0.0)
    msg_abort(11007, "sigma8 is not given in the parameter file");
  else
    msg_printf(normal, "sigma8 (check) = %lg\n", param->sigma8);

  msg_printf(normal, "np_alloc_factor = %lg\n", param->np_alloc_factor);

  msg_printf(normal, "loglevel = %d\n", param->loglevel);

  if(param->power_spectrum_filename == 0)
    msg_abort(11008, "powerspectrum is not given in the parameter file\n");
  else if(param->strlen_power_spectrum_filename <= 1)
    msg_abort(11018, "Error: powerspectrum filename is empty\n");
  else
    msg_printf(normal, "powerspectrum = \"%s\"\n",
	       param->power_spectrum_filename);

  if(param->fof_filename == 0)
    msg_printf(normal, "No FoF halo finding. fof_filename is not given.\n");
  else if(param->strlen_fof_filename <= 1)
    msg_abort(11019, "Error: fof_filename is empty\n");
  else {
    msg_printf(normal, "fof_filename = \"%s\"\n", param->fof_filename);
    msg_printf(normal, "linking_factor = %lg\n", param->fof_linking_factor);
  }

  if(param->snapshot_filename == 0)
    msg_printf(normal, "No N-body particle output (snapshot)\n");
  else if(param->strlen_snapshot_filename <= 1)
    msg_abort(11020, "Error: snapshot filename is empty.\n");
  else {
    msg_printf(normal, "snapshot = \"%s\"\n", param->snapshot_filename);
    msg_printf(normal, "write_longid = %d\n", param->write_longid);
  }

  if(param->subsample_filename == 0)
    msg_printf(normal, "No subsampled particle output (subsample)\n");
  else if(param->strlen_subsample_filename <= 1)
    msg_abort(10021, "Error: subsample filename is empty.\n"); 
  else {
    msg_printf(normal, "subsample = \"%s\"\n", param->subsample_filename);
    msg_printf(normal, "subsample_factor = %lg\n", param->subsample_factor);
  }

  if(param->cgrid_filename == 0)
    msg_printf(normal, "No coarse density grid output (coarse_grid)\n");
  else if(param->strlen_cgrid_filename <= 1)
    msg_abort(10022, "coarse_grid filename is empty.\n");
  else {
    msg_printf(normal, "coarse_grid = \"%s\"\n", param->cgrid_filename);
    msg_printf(normal, "coarse_grid_nc = %d\n", param->cgrid_nc);
  }

  if(param->init_filename) {
    if(param->strlen_init_filename <= 1)
      msg_abort(10023, "init filename is empty.\n");

    msg_printf(normal, "inital= \"%s\"\n", param->init_filename);
    if(param->snapshot_filename == 0)
      msg_printf(normal, "write_longid= %d\n", param->write_longid);
  }
}

