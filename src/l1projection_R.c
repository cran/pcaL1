#include "Clp_C_Interface.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <Rmath.h>
#include "type.h"


static int allocateMemoryProj (ENTITYINFOptr entityinfo, PROBLEMINFOptr probleminfo); 
int solveL1Projection (ENTITYINFOptr entityinfo, SOLVERINFOptr solverinfo, PROBLEMINFOptr probleminfo);

static void
   free_and_null (char **ptr);

void l1projection (double *points_XT, int *dataDim, int *q, double *PCs, double *projPoints, double *alphas) 
{

  ENTITYINFO entityinfo;
  SOLVERINFO  solverinfo;
  PROBLEMINFO  probleminfo;

  probleminfo.status = 0;
  int status = probleminfo.status;
  

  solverinfo.model        = NULL;

  probleminfo.aind        = NULL;
  probleminfo.rhs         = NULL;
  probleminfo.matbeg      = NULL;
  probleminfo.matval      = NULL;
  probleminfo.matind      = NULL;
  probleminfo.obj         = NULL;
  probleminfo.lb          = NULL;
  probleminfo.ub          = NULL;
  /*probleminfo.colname     = NULL;*/
  probleminfo.projPoints  = projPoints;
  probleminfo.alphas      = alphas;
  probleminfo.projdim     = *q;/*desired number of PCs*/
 
  entityinfo.numentities_n   = dataDim[1];
  entityinfo.numattributes_m = dataDim[0];
  entityinfo.PCs          = PCs;
 
  entityinfo.points_XT = points_XT;/*transpose of data matrix*/

  status = allocateMemoryProj (&entityinfo, &probleminfo); /* at the end of this file */ 
  if (status) {
    REprintf ("Unable to allocate memory\n");
    goto TERMINATE;
  }

  status = solveL1Projection ( &entityinfo, &solverinfo, &probleminfo); /* in l1projection.c*/
  if (status) {
    REprintf ("Unable to solve.  Terminating...; or done\n");
    goto TERMINATE;
  }

  TERMINATE:
  
  for (probleminfo.i = 0; probleminfo.i < entityinfo.numentities_n; ++probleminfo.i) {
    free_and_null ((char **) &probleminfo.aind[probleminfo.i]);
  }
  free_and_null ((char **) &probleminfo.aind);
  free_and_null ((char **) &probleminfo.rhs);
  free_and_null ((char **) &probleminfo.matbeg);
  free_and_null ((char **) &probleminfo.matval);
  free_and_null ((char **) &probleminfo.matind);
  free_and_null ((char **) &probleminfo.obj);
  free_and_null ((char **) &probleminfo.lb);
  free_and_null ((char **) &probleminfo.ub);
  /*for (probleminfo.i = 0; probleminfo.i < entityinfo.numentities_n*probleminfo.projdim + 2*entityinfo.numentities_n*entityinfo.numattributes_m; ++probleminfo.i) {
    free_and_null ((char **) &probleminfo.colname[probleminfo.i]);
  }
  free_and_null ((char **) &probleminfo.colname);*/

  if ( solverinfo.model != NULL ) {
    Clp_deleteModel (solverinfo.model);
    if ( status ) {
      REprintf ("CPXfreeprob failed, error code %d.\n", status);
    }
  }
}
    
static void
free_and_null (char **ptr) {
  if ( *ptr != NULL ) {
     free (*ptr);
     *ptr = NULL;
  }
} /* END free_and_null */  

static int allocateMemoryProj (ENTITYINFOptr entityinfo, PROBLEMINFOptr probleminfo) {
  int numentities_n   = entityinfo->numentities_n;
  int numattributes_m = entityinfo->numattributes_m;
  
  int i = probleminfo->i;
  int projdim = probleminfo->projdim;
  int numcols = probleminfo->numcols;


  /* allocate memory for columns */
  probleminfo->aind   = (int **) malloc ((long unsigned int)numentities_n*sizeof(int *));
  for (i = 0; i < numentities_n; ++i) {
    probleminfo->aind[i]   = (int *) malloc ((long unsigned int)projdim*sizeof(int));
  }

  numcols = numentities_n*projdim + 2*numentities_n*numattributes_m;
  probleminfo->obj     = (double *) malloc ((long unsigned int)(numcols) * sizeof (double));
  probleminfo->lb      = (double *) malloc ((long unsigned int)(numcols) * sizeof (double));
  probleminfo->ub      = (double *) malloc ((long unsigned int)(numcols) * sizeof (double));
  /*probleminfo->colname = (char **)  malloc ((long unsigned int)(numcols) * sizeof (char *));
  for (i = 0; i < numcols; ++i) {
    probleminfo->colname[i] = (char *) malloc (20 * sizeof (char));
  }*/

  /* allocate memory for constraints */
  probleminfo->rhs    = (double *) malloc ((long unsigned int)numentities_n*numattributes_m*sizeof (double));
  probleminfo->matbeg = (int *) malloc ((long unsigned int)(numcols+1)* sizeof (int));
  probleminfo->matind = (int *) malloc ((long unsigned int)(numentities_n*projdim*numattributes_m + 2*numentities_n*numattributes_m) * sizeof (int));
  probleminfo->matval = (double *) malloc ((long unsigned int)(numentities_n*projdim*numattributes_m + 2*numentities_n*numattributes_m) * sizeof (double));
  return 0;
  
} /* end allocateMemory */
