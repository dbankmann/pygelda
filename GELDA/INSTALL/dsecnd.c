/* dsecnd.f -- translated by f2c (version 19940714).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"
#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>

doublereal dsecnd_()
{
    /* System generated locals */
    doublereal ret_val;

    /* Local variables */
    int i;
    struct rusage usage;

/*  -- LAPACK auxiliary routine (version 2.0) -- */
/*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
/*     Courant Institute, Argonne National Lab, and Rice University */
/*     September 30, 1994 */

/*  Purpose */
/*  ======= */

/*  DSECND returns the user time for a process in seconds. */
/*  This version gets the time from the system function ETIME. */

/* ===================================================================== 
*/

/*     .. Local Scalars .. */
/*     .. */
/*     .. Local Arrays .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. Executable Statements .. */

    i = getrusage (RUSAGE_SELF, &usage);
    ret_val = usage.ru_utime.tv_sec;
    ret_val = ret_val + usage.ru_utime.tv_usec / 1000000.0;
    return ret_val;

/*     End of DSECND */

} /* dsecnd_ */




