#include <math.h>

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"

#include <GaussLobatto.h>

void SetupCoords(CCTK_ARGUMENTS)
{
    DECLARE_CCTK_ARGUMENTS
    DECLARE_CCTK_PARAMETERS

    int i;
    int status;
    int gsh[1];
    CCTK_INT coord_handle, varindex, ierr;

    /* Get the array bounds for the grid functions. Since this coordinate
     * implementation implementation only provides a 1-D Cartesian grid, we'll
     * only use gsh[0]. */
    status = CCTK_GroupgshVN(cctkGH, 1, gsh, "spectral_method::x");

    for (i=0; i<gsh[0]; i++)
    {
        x[i] = gauss_lobatto_point(i, gsh[0]);
    }

    coord_handle = Coord_CoordHandle(cctkGH, "x", "cart1d");
    varindex = CCTK_VarIndex("spectral_method::x");
    ierr = Util_TableSetReal(coord_handle, -1.0, "PHYSICALMIN");
    ierr = Util_TableSetReal(coord_handle, -1.0, "COMPMIN");
    ierr = Util_TableSetReal(coord_handle, +1.0, "PHYSICALMAX");
    ierr = Util_TableSetReal(coord_handle, +1.0, "COMPMAX");
    ierr = Util_TableSetString(coord_handle, "nonuniform", "TYPE");
    ierr = Util_TableSetString(coord_handle, "no", "TIMEDEPENDENT");
    ierr = Util_TableSetString(coord_handle, "CCTK_REAL", "DATATYPE");
    ierr = Util_TableSetInt(coord_handle, varindex, "GAINDEX");

}
