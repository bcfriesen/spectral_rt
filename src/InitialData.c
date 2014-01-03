#include <math.h>

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"

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
