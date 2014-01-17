#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "util_Table.h"

#include "CoordBase.h"

/* Registers the 1-D Cartesian coordinate system. */
void Cart1D_RegisterCoordinates(CCTK_ARGUMENTS)
{
    DECLARE_CCTK_ARGUMENTS
    DECLARE_CCTK_PARAMETERS

    CCTK_INT ierr;
    CCTK_INT systemhandle;

    /* Set up a new 1-D coordinate system called 'cart1d'. */
    ierr = Coord_SystemRegister(cctkGH, 1, "cart1d");

    /* Get the handle for cart1d so we can add stuff to it. */
    systemhandle = Coord_SystemHandle(cctkGH, "cart1d");

    /* Register the variable 'x' as a coordinate along the first (and only)
     * dimension of this coordinate system. */
    ierr = Coord_CoordRegister(cctkGH, systemhandle, 1, "x");

    /* Set cart1d as the default coordinate system for all grid variables with
     * 1 dimension. This can be overridden by specifying the desired coordinate
     * system in the desired thorn. (See the CoordBase documentation for
     * details.) */
    systemhandle = Coord_SetDefaultSystem(cctkGH, "cart1d");
}
