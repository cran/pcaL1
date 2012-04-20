#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "type.h"


int allocateMemory2 (ENTITYINFOptr entityinfo, PROBLEMINFOptr probleminfo); 

int solveL1PCA (ENTITYINFOptr entityinfo, SOLVERINFOptr solverinfo, PROBLEMINFOptr probleminfo);

static void
   free_and_null (char **ptr);

void l1pca (double *points_XT, int *dataDim, int *q, double *tolerance, int *iterations, double *initV, double *PCs, double *Scores) 
{

  ENTITYINFO entityinfo;
  SOLVERINFO  solverinfo;
  PROBLEMINFO  probleminfo;

  int status = probleminfo.status;
  int i = probleminfo.i;
  status = 0;

  entityinfo.points_XT       = points_XT;
  entityinfo.numentities_n  = dataDim[1];
  entityinfo.numattributes_m= dataDim[0];

  solverinfo.modelU        = NULL;
  solverinfo.modelV        = NULL;


  probleminfo.q = *q; 
  probleminfo.initV = initV;
  probleminfo.tolerance  = *tolerance;
  probleminfo.iterations = *iterations;
  probleminfo.rhsL    = NULL;
  probleminfo.rhsU    = NULL;
  probleminfo.matbeg  = NULL;
  probleminfo.matval  = NULL;
  probleminfo.matind  = NULL;
  probleminfo.obj     = NULL;
  probleminfo.ub      = NULL;
  probleminfo.lb      = NULL;
  probleminfo.colname = NULL;   
  probleminfo.nu      = NULL;
  probleminfo.nv      = NULL;
  /*probleminfo.sigma   = NULL;*/
  probleminfo.Vtemp   = NULL;
  probleminfo.Utemp   = NULL;
  probleminfo.rowind  = NULL;
  probleminfo.xinda   = NULL;
  probleminfo.xindb   = NULL;
  probleminfo.vind   = NULL;
  /*probleminfo.Vsol    = NULL;
  probleminfo.Usol    = NULL;*/

  probleminfo.V = PCs;
  probleminfo.U = Scores;      

  status = allocateMemory2 (&entityinfo, &probleminfo); /* at the end of this file */ 
  if (status) {
    REprintf ("Unable to allocate memory\n");
    goto TERMINATE;
  }

  status = solveL1PCA (&entityinfo, &solverinfo, &probleminfo); /* in solveproblem.c*/
  if (status) {
    REprintf ("Unable to solve.  Terminating...; or done\n");
    goto TERMINATE;
  }
  REprintf ("\n");

TERMINATE:
  
  free_and_null ((char **) &probleminfo.rhsL);
  free_and_null ((char **) &probleminfo.rhsU);
  free_and_null ((char **) &probleminfo.matbeg);
  free_and_null ((char **) &probleminfo.matval);
  free_and_null ((char **) &probleminfo.matind);
  free_and_null ((char **) &probleminfo.obj);
  free_and_null ((char **) &probleminfo.lb);
  free_and_null ((char **) &probleminfo.ub);
  for (i = 0; i < 2*entityinfo.numentities_n*entityinfo.numattributes_m + entityinfo.numattributes_m*probleminfo.q; ++i) {
    free_and_null ((char **) &probleminfo.colname[i]);
  }
  free_and_null ((char **) &probleminfo.colname);
  free_and_null ((char **) &probleminfo.nu);
  free_and_null ((char **) &probleminfo.nv);
  /*free_and_null ((char **) &probleminfo.sigma);*/
  free_and_null ((char **) &probleminfo.Vtemp);
  free_and_null ((char **) &probleminfo.uind);
  /*free_and_null ((char **) &probleminfo.Vsol);*/
  free_and_null ((char **) &probleminfo.Utemp);
  for (i = 0; i < entityinfo.numentities_n; ++i) {
    free_and_null ((char **) &probleminfo.rowind[i]);
  }
  free_and_null ((char **) &probleminfo.rowind);
  free_and_null ((char **) &probleminfo.xinda);
  free_and_null ((char **) &probleminfo.xindb);
  for (i = 0; i < entityinfo.numattributes_m; ++i) {
    free_and_null ((char **) &probleminfo.vind[i]);
  }
  free_and_null ((char **) &probleminfo.vind);
  if ( solverinfo.modelU != NULL ) {
    Clp_deleteModel (solverinfo.modelU);
    if ( status ) {
      REprintf ("Clp delete failed, error code %d.\n", status);
    }
  }
  if ( solverinfo.modelV != NULL ) {
    Clp_deleteModel (solverinfo.modelV);
    if ( status ) {
      REprintf ("Clp delete failed, error code %d.\n", status);
    }
  }
  /*free_and_null ((char **) &probleminfo.Usol);*/
 /*return (status);*/
}
    
static void
free_and_null (char **ptr)
{
   if ( *ptr != NULL ) {
      free (*ptr);
      *ptr = NULL;
   }
} /* END free_and_null */  

int allocateMemory2 (ENTITYINFOptr entityinfo, PROBLEMINFOptr probleminfo) {
  int numentities_n   = entityinfo->numentities_n;
  int numattributes_m = entityinfo->numattributes_m;
  
  int numcols = probleminfo->numcols;
  int q       = probleminfo->q;
  int rcnt    = probleminfo->rcnt;
  int nzcnt   = probleminfo->nzcnt;
  
  int i = probleminfo->i;
  int j = probleminfo->j;


  /* allocate memory for columns */
  numcols = 2*numentities_n*numattributes_m + numattributes_m*q; /* number of columns when solving for V */
  probleminfo->obj     = (double *) malloc ((long unsigned int) numcols * sizeof (double));
  probleminfo->lb      = (double *) malloc ((long unsigned int) numcols * sizeof (double));
  probleminfo->ub      = (double *) malloc ((long unsigned int) numcols * sizeof (double));
  probleminfo->colname = (char **)  malloc ((long unsigned int) numcols * sizeof (char *));
  for (i = 0; i < numcols; ++i) {
    probleminfo->colname[i] = (char *) malloc (20 * sizeof (char));
  }

  /* allocate memory for constraints */
  rcnt = numentities_n*numattributes_m; /* number of rows when solving for V */
  nzcnt = 2*numentities_n*numattributes_m + numentities_n*numattributes_m*q; /* number of non-zeros when sovling for V */
  probleminfo->rhsL    = (double *) malloc ((long unsigned int) rcnt * sizeof (double));
  probleminfo->rhsU    = (double *) malloc ((long unsigned int) rcnt * sizeof (double));
  probleminfo->xinda   = (int *) malloc ((long unsigned int)rcnt * sizeof (int));
  probleminfo->xindb   = (int *) malloc ((long unsigned int)rcnt * sizeof (int));
  probleminfo->matbeg = (int *)    malloc ((long unsigned int)(numcols+1) * sizeof (int));
  probleminfo->matind = (int *)    malloc ((long unsigned int)nzcnt * sizeof (int));
  probleminfo->matval = (double *) malloc ((long unsigned int)nzcnt * sizeof (double));

  probleminfo->uind = (int *) malloc((long unsigned int)q*sizeof(int));
  
  probleminfo->Vtemp = (double *) malloc ((long unsigned int)numattributes_m*(long unsigned int) q*sizeof(double));
  probleminfo->Utemp = (double *) malloc ((long unsigned int)numentities_n*(long unsigned int)q*sizeof(double));
  probleminfo->nu    = (double *) malloc ((long unsigned int)q*sizeof(double));
  probleminfo->nv    = (double *) malloc ((long unsigned int)q*sizeof(double));
  /*probleminfo->sigma = (double *) malloc (q*sizeof(double));*/

  probleminfo->rowind = (int **) malloc((long unsigned int)numentities_n *sizeof(int *));
  for (i = 0; i < numentities_n; ++i) {
    probleminfo->rowind[i] = (int *) malloc((long unsigned int)numattributes_m * sizeof(int));
  }
  probleminfo->vind = (int **) malloc((long unsigned int)numattributes_m * sizeof(int *));
  for (j = 0; j < numattributes_m; ++j) {
    probleminfo->vind[j] = (int *) malloc((long unsigned int)q * sizeof(int));
  }
    
  return 0;
} /* end allocateMemory */
