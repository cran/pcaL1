#include "type.h"

static int allocateMemoryL1rs (ENTITYINFOptr entityinfo, PROBLEMINFOptr probleminfo); 
int solveSharpeL1rs (ENTITYINFOptr entityinfo, SOLVERINFOptr solverinfo, PROBLEMINFOptr probleminfo);

static void
   free_and_null (char **ptr);

void sharpel1rs(double *points_XT, int *dataDim, int *q, double *PCs, double *objectives) 
{
 
  ENTITYINFO entityinfo;
  SOLVERINFO  solverinfo;
  PROBLEMINFO  probleminfo;
  
  probleminfo.status = 0;
  int status        = probleminfo.status;

  solverinfo.model = NULL;

  probleminfo.obj       = NULL;
  probleminfo.lb        = NULL;
  probleminfo.ub        = NULL;
  probleminfo.rhs       = NULL;
  probleminfo.matbeg    = NULL;
  probleminfo.matind    = NULL;
  probleminfo.matval    = NULL;
  probleminfo.colname   = NULL;
  probleminfo.PCs       = PCs;
  /*probleminfo.getScores = *getScores;
  probleminfo.scores    = scores;*/
  probleminfo.objectives = objectives;
 
  entityinfo.numentities_n   = dataDim[1];
  entityinfo.numattributes_m = dataDim[0];
  
  entityinfo.points_XT = points_XT; /* transpose of data matrix */

  probleminfo.q = *q; /* desired number of PCs */

  status = allocateMemoryL1rs(&entityinfo, &probleminfo);
  if (status) {
    REprintf ("Unable to allocate memory\n");
    goto TERMINATE;
  }

  status = solveSharpeL1rs( &entityinfo, &solverinfo, &probleminfo); /* in l1line.c*/
  if (status) {
    REprintf ("Unable to solve.  Terminating...; or done\n");
    goto TERMINATE;
  }

  
TERMINATE:

  free_and_null ((char **) &probleminfo.ratios);     
  free_and_null ((char **) &probleminfo.weights);     
  free_and_null ((char **) &probleminfo.v);     
  free_and_null ((char **) &probleminfo.tosort);

} /* end l1line */

static void
free_and_null (char **ptr) {
  if ( *ptr != NULL ) {
     free (*ptr);
     *ptr = NULL;
  }
} /* END free_and_null */  

static int allocateMemoryL1rs (ENTITYINFOptr entityinfo, PROBLEMINFOptr probleminfo) {
  int numentities_n   = entityinfo->numentities_n;
  int numattributes_m = entityinfo->numattributes_m;

  probleminfo->ratios     = (double *) malloc (numentities_n*sizeof (double));
  probleminfo->tosort     = (double **) malloc (numentities_n*sizeof (double *));
  probleminfo->weights    = (double *) malloc (numentities_n*sizeof (double));
  probleminfo->v = (double *) malloc(numattributes_m*sizeof(double));

  return 0;
} /* end allocateMemoryL1rs */

