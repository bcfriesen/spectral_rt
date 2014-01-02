 /*@@
   @file      Startup.c
   @date
   @author    Brian Friesen
   @desc
              Register Galerkin method banner
   @enddesc
   @version $Header$
 @@*/

#include "cctk.h"

static const char *rcsid = "$Header$";

CCTK_FILEVERSION(SpectralRT_Galerkin_Startup_c)

int Galerkin_Startup(void);

 /*@@
   @routine    Galerkin_Startup
   @date
   @author     Brian Friesen
   @desc

   @enddesc
   @calls
   @calledby
   @history

   @endhistory

@@*/
int Galerkin_Startup(void)
{

   const char *banner = "SpectralRT: Calculate spectral coefficients for RTE";

   CCTK_RegisterBanner(banner);

   return 0;
}
