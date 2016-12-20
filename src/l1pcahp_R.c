#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "type.h"

int allocateMemoryHp (ENTITYINFOptr entityinfo, PROBLEMINFOptr probleminfo);

int solveL1PCAHp (ENTITYINFOptr entityinfo, SOLVERINFOptr solverinfo, PROBLEMINFOptr probleminfo);

static void
   free_and_null (char **ptr);

void l1pcahp (double *points_XT, int *dataDim, double *threshold, double *initV, double *PCs)
{

  ENTITYINFO entityinfo;
  SOLVERINFO  solverinfo;
  PROBLEMINFO  probleminfo;

  probleminfo.status = 0;
  probleminfo.k = 0;
  int status = probleminfo.status;

  entityinfo.points_XT       = points_XT;
  entityinfo.numentities_n   = dataDim[1];
  entityinfo.numattributes_m = dataDim[0];
  probleminfo.initV = initV;


  probleminfo.numfactors = entityinfo.numattributes_m;

  solverinfo.model        = NULL;


  probleminfo.xx_obj      = 0;
  probleminfo.x_obj       = 0;
  probleminfo.currObj     = 0;
  probleminfo.point       = NULL;
  probleminfo.alpha       = NULL;
  probleminfo.numcols     = 0;
  probleminfo.numrows     = 0 ;
  probleminfo.obj         = NULL;
  probleminfo.PCs         = PCs;
  probleminfo.threshold   = *threshold;

  probleminfo.matbeg      = NULL;
  probleminfo.matval      = NULL;
  probleminfo.matind      = NULL;
  probleminfo.obj         = NULL;
  probleminfo.lb          = NULL;
  probleminfo.ub          = NULL;
  probleminfo.rhsL        = NULL;
  probleminfo.rhsU        = NULL;


  status = allocateMemoryHp (&entityinfo, &probleminfo); /* at the end of this file */
  if (status) {
    REprintf ("Unable to allocate memory\n");
    goto TERMINATE;
  }

  status = solveL1PCAHp (&entityinfo, &solverinfo, &probleminfo); /* in solveproblem.c*/

  if (status) {
    REprintf ("Unable to solve.  Terminating...; or done\n");
    goto TERMINATE;
  }

TERMINATE:

  free_and_null ((char **) &probleminfo.point);
  free_and_null ((char **) &probleminfo.rhsL);
  free_and_null ((char **) &probleminfo.rhsU);
  free_and_null ((char **) &probleminfo.matbeg);
  free_and_null ((char **) &probleminfo.matval);
  free_and_null ((char **) &probleminfo.matind);
  free_and_null ((char **) &probleminfo.obj);
  free_and_null ((char **) &probleminfo.lb);
  free_and_null ((char **) &probleminfo.ub);
  free_and_null ((char **) &probleminfo.alpha);
  if ( solverinfo.model != NULL ) {
    Clp_deleteModel (solverinfo.model);

    if ( status ) {
      REprintf ("Clp delete failed, error code %d.\n", status);
    }
  }
}
    
static void
free_and_null (char **ptr)
{
   if ( *ptr != NULL ) {
      free (*ptr);
      *ptr = NULL;
   }
} /* END free_and_null */  

int allocateMemoryHp (ENTITYINFOptr entityinfo, PROBLEMINFOptr probleminfo) {
  int numentities_n   = entityinfo->numentities_n;
  int numattributes_m = entityinfo->numattributes_m;
  
  int numcols = probleminfo->numcols;
  int numrows = probleminfo->numrows;
  int nzcnt   = probleminfo->nzcnt;

  numcols = numattributes_m + numentities_n;
  numrows = 2 * numentities_n ;
  nzcnt = 2 * numentities_n * (numattributes_m + 1);

  probleminfo->point   = (double *) malloc ((long unsigned int) numattributes_m * sizeof (double));


  /* allocate memory for columns */
  probleminfo->obj     = (double *) malloc ((long unsigned int) numcols * sizeof (double));
  probleminfo->lb      = (double *) malloc ((long unsigned int) numcols * sizeof (double));
  probleminfo->ub      = (double *) malloc ((long unsigned int) numcols * sizeof (double));
  probleminfo->alpha   = (double *) malloc ((long unsigned int) numcols * sizeof (double));

  /* allocate memory for constraints */

  probleminfo->rhsL    = (double *) malloc ((long unsigned int) numrows * sizeof (double));
  probleminfo->rhsU    = (double *) malloc ((long unsigned int) numrows * sizeof (double));
  probleminfo->matbeg = (int *)    malloc ((long unsigned int) (numcols+1) * sizeof (int));
  probleminfo->matind = (int *)    malloc ((long unsigned int)nzcnt * sizeof (int));
  probleminfo->matval = (double *) malloc ((long unsigned int)nzcnt * sizeof (double));

    
  return 0;
} /* end allocateMemory */
