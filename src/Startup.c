#include "cctk.h"

int Galerkin_Startup(void)
{

   const char *banner = "SpectralRT: Calculate spectral coefficients for RTE";

   CCTK_RegisterBanner(banner);

   return 0;
}
