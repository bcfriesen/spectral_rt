#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"

void Calculate_Integration_Weights(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  int i;
  int status;
  int gsh[3];

  status = CCTK_GroupgshVN(cctkGH, 1, gsh, "spectral_method::x");

  if (CCTK_EQUALS(integration_type,"trapezoid"))
  {
    for (i=0; i<gsh[0]; i++)
    {
        int_weights[i] = 69.0;
    }
  }

}
