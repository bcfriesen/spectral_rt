#include <math.h>

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"

void SetupCoords(CCTK_ARGUMENTS)
{
    DECLARE_CCTK_ARGUMENTS
    DECLARE_CCTK_PARAMETERS

    int i;
    int status;
    int gsh[1];
    double pi;
    CCTK_INT coord_handle, varindex, ierr;

    /* Get the array bounds for the grid functions. Since the 'grid1d'
     * implementation only provides a 1-D Cartesian grid, we'll only use
     * gsh[0]. */
    status = CCTK_GroupgshVN(cctkGH, 1, gsh, "grid1d::x");

    pi = acos(-1.0);

    for (i=0; i<gsh[0]; i++)
    {
        x[i] = cos(pi * (double)i / (double)gsh[0]);
    }

    coord_handle = Coord_CoordHandle(cctkGH, "x", "cart1d");
    varindex = CCTK_VarIndex("grid1d::x");
    ierr = Util_TableSetReal(coord_handle, -1.0, "PHYSICALMIN");
    ierr = Util_TableSetReal(coord_handle, -1.0, "COMPMIN");
    ierr = Util_TableSetReal(coord_handle, +1.0, "PHYSICALMAX");
    ierr = Util_TableSetReal(coord_handle, +1.0, "COMPMAX");
    ierr = Util_TableSetString(coord_handle, "nonuniform", "TYPE");
    ierr = Util_TableSetString(coord_handle, "no", "TIMEDEPENDENT");
    ierr = Util_TableSetString(coord_handle, "CCTK_REAL", "DATATYPE");
    ierr = Util_TableSetInt(coord_handle, varindex, "GAINDEX");

}
