 /*@@
   @file    InitialData.c
   @date
   @author  Brian Friesen
   @desc
            Initial data for spectral method radiative transfer
   @enddesc
 @@*/

#include <math.h>

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"

static const char *rcsid = "$Header$";

CCTK_FILEVERSION(SpectralRT_Galerkin_InitialData_c)

void Galerkin_InitialData(CCTK_ARGUMENTS);


 /*@@
   @routine    Galerkin_InitialData
   @date
   @author     Brian Friesen
   @desc
               Set up initial data for the spectral coefficients
   @enddesc
@@*/

void Galerkin_InitialData(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  int i;

  for (i=0; i<num_basis_funcs; i++)
  {
    sc[i] = 0.0;
  }
}
