#include "Clp_C_Interface.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <Rmath.h>
#include "type.h"


int allocateMemory (ENTITYINFOptr entityinfo, PROBLEMINFOptr probleminfo); 
int solveL1PCAStar (ENTITYINFOptr entityinfo, SOLVERINFOptr solverinfo, PROBLEMINFOptr probleminfo);

static void
   free_and_null (char **ptr);

void l1pcastar (double *points_XT, int *dataDim, int *q, int *getScores, int *getProjPoints, double *PCs, double *scores, double *projPoints) 
{

  ENTITYINFO entityinfo;
  SOLVERINFO  solverinfo;
  PROBLEMINFO  probleminfo;

  probleminfo.status = 0;
  int status = probleminfo.status;
  

  solverinfo.model        = NULL;

  probleminfo.betaind   = NULL;
  probleminfo.eplusind  = NULL; 
  probleminfo.eminusind = NULL;
  probleminfo.beta      = NULL;
  probleminfo.bestdir   = NULL;
  probleminfo.xpluslambda_Z = NULL;
  probleminfo.xpluslambda_Z2 = NULL;
  probleminfo.S           = NULL;
  probleminfo.VT          = NULL;
  probleminfo.Vj         = NULL;
  probleminfo.preVj      = NULL;
  probleminfo.temppreVj  = NULL;
  probleminfo.VBeta       = NULL;
  probleminfo.preVBeta    = NULL;
  probleminfo.temppreVBeta= NULL;
  probleminfo.tempPC      = NULL;
  probleminfo.rhs         = NULL;
  probleminfo.matbeg      = NULL;
  probleminfo.matval      = NULL;
  probleminfo.matind      = NULL;
  probleminfo.obj         = NULL;
  probleminfo.lb          = NULL;
  probleminfo.ub          = NULL;
  probleminfo.colname     = NULL;
  probleminfo.a           = NULL;
  probleminfo.b           = PCs;
  probleminfo.projPoints  = projPoints;
  probleminfo.scores      = scores;
  probleminfo.work        = NULL;
  probleminfo.getScores     = *getScores;
  probleminfo.getProjPoints = *getProjPoints;
 
  entityinfo.numentities_n   = dataDim[1];
  entityinfo.numattributes_m = dataDim[0];

  if (entityinfo.numattributes_m <= entityinfo.numentities_n) {
    probleminfo.numfactors = entityinfo.numattributes_m;
  }
  else {
    probleminfo.numfactors = entityinfo.numentities_n;
  }
  
  probleminfo.q = *q;/*desired number of PCs*/
 
  entityinfo.points_XT = points_XT;/*transpose of data matrix*/

  if ((VERBOSITY) >= 2) {
    REprintf("n %d m %d\n", dataDim[1], dataDim[0]);
    REprintf("getScores %d getProjPoints %d\n", *getScores, *getProjPoints);
    REprintf("getScores %d getProjPoints %d\n", probleminfo.getScores, probleminfo.getProjPoints);
    REprintf("numfactors %d q %d\n", probleminfo.numfactors, probleminfo.q);
  }
 
  status = allocateMemory (&entityinfo, &probleminfo); /* at the end of this file */ 
  if (status) {
    REprintf ("Unable to allocate memory\n");
    goto TERMINATE;
  }

  status = solveL1PCAStar ( &entityinfo, &solverinfo, &probleminfo); /* in l1pcaStar.c*/
  if (status) {
    REprintf ("Unable to solve.  Terminating...; or done\n");
    goto TERMINATE;
  }

  REprintf ("\n");
  
TERMINATE:
  
  free_and_null ((char **) &probleminfo.betaind);
  free_and_null ((char **) &probleminfo.eplusind);
  free_and_null ((char **) &probleminfo.eminusind);
  free_and_null ((char **) &probleminfo.beta);
  free_and_null ((char **) &probleminfo.bestdir);
  free_and_null ((char **) &probleminfo.xpluslambda_Z);
  free_and_null ((char **) &probleminfo.xpluslambda_Z2);
  free_and_null ((char **) &probleminfo.S);
  free_and_null ((char **) &probleminfo.VT);
  free_and_null ((char **) &probleminfo.Vj);
  free_and_null ((char **) &probleminfo.preVj);
  free_and_null ((char **) &probleminfo.temppreVj);
  free_and_null ((char **) &probleminfo.VBeta);
  free_and_null ((char **) &probleminfo.preVBeta);
  free_and_null ((char **) &probleminfo.temppreVBeta);
  free_and_null ((char **) &probleminfo.tempPC);
  free_and_null ((char **) &probleminfo.rhs);
  free_and_null ((char **) &probleminfo.matbeg);
  free_and_null ((char **) &probleminfo.matval);
  free_and_null ((char **) &probleminfo.matind);
  free_and_null ((char **) &probleminfo.obj);
  free_and_null ((char **) &probleminfo.lb);
  free_and_null ((char **) &probleminfo.ub);
  for (probleminfo.i = 0; probleminfo.i < probleminfo.numfactors + 2 * entityinfo.numentities_n; ++probleminfo.i) {
    free_and_null ((char **) &probleminfo.colname[probleminfo.i]);
  }
  free_and_null ((char **) &probleminfo.colname);
  free_and_null ((char **) &probleminfo.a);
  free_and_null ((char **) &probleminfo.work);
  /*free_and_null ((char **) &probleminfo.projpoints);
  free_and_null ((char **) &probleminfo.scores);*/

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

int allocateMemory (ENTITYINFOptr entityinfo, PROBLEMINFOptr probleminfo) {
  int numentities_n   = entityinfo->numentities_n;
  int numattributes_m = entityinfo->numattributes_m;
  
  int i = probleminfo->i;
  int j = probleminfo->j;
  int numfactors = probleminfo->numfactors;

  /* allocate memory for best regression */
  probleminfo->beta      = (double *) malloc ((long unsigned int) numfactors * sizeof (double)); /* keep track of best linear regression */
  probleminfo->bestdir   = (int *) malloc ((long unsigned int)numfactors * sizeof (int)); /* to keep track of projections */

  /* allocate memory for columns */
  probleminfo->betaind   = (int *) malloc ((long unsigned int)numfactors * sizeof (int));
  probleminfo->eplusind  = (int *) malloc ((long unsigned int)numentities_n * sizeof (int));
  probleminfo->eminusind = (int *) malloc ((long unsigned int)numentities_n * sizeof (int));

  probleminfo->obj     = (double *) malloc ((long unsigned int)(numfactors + 2 * numentities_n) * sizeof (double));
  probleminfo->lb      = (double *) malloc ((long unsigned int)(numfactors + 2 * numentities_n) * sizeof (double));
  probleminfo->ub      = (double *) malloc ((long unsigned int)(numfactors + 2 * numentities_n) * sizeof (double));
  probleminfo->colname = (char **)  malloc ((long unsigned int)(numfactors + 2 * numentities_n) * sizeof (char *));
  for (i = 0; i < numfactors + 2 * numentities_n; ++i) {
    probleminfo->colname[i] = (char *) malloc (20 * sizeof (char));
  }

  /* allocate memory for constraints */
  probleminfo->rhs    = (double *) malloc ((long unsigned int)numentities_n * sizeof (double));
  probleminfo->matbeg = (int *) malloc ((long unsigned int)(numfactors + 2*numentities_n + 1)* sizeof (int));
  probleminfo->matind = (int *)    malloc ((long unsigned int)(numfactors * numentities_n + 2 * numentities_n) * sizeof (int));
  probleminfo->matval = (double *) malloc ((long unsigned int)(numfactors * numentities_n + 2 * numentities_n) * sizeof (double));
  
  /* allocate memory for chgbds-no more changing one bound at a time */

  /* allocate memory for projected points */
  probleminfo->xpluslambda_Z  = (double *) malloc ((long unsigned int)numentities_n * (long unsigned int)numfactors * sizeof (double));
  probleminfo->xpluslambda_Z2 = (double *) malloc ((long unsigned int)numentities_n * (long unsigned int)numfactors * sizeof (double));
  probleminfo->lwork        = 9 *(numentities_n + numattributes_m) * (NBMAX);
  probleminfo->work         = (double *) malloc ((long unsigned int)probleminfo->lwork * sizeof (double));
  probleminfo->S            = (double *) malloc ((long unsigned int)numattributes_m * sizeof (double));
  probleminfo->Vj          = (double *) malloc ((long unsigned int)numattributes_m * (long unsigned int)numattributes_m * sizeof (double));
  probleminfo->VT           = (double *) malloc ((long unsigned int)numentities_n * (long unsigned int)numentities_n * sizeof (double));
  probleminfo->a            = (double *) malloc ((long unsigned int)numfactors * sizeof (double));
  probleminfo->VBeta        = (double *) malloc ((long unsigned int)numattributes_m * (long unsigned int)numattributes_m * sizeof (double));
  probleminfo->preVj       = (double *) malloc ((long unsigned int)numattributes_m * (long unsigned int)numattributes_m * sizeof (double));
  probleminfo->preVBeta     = (double *) malloc ((long unsigned int)numattributes_m * (long unsigned int)numattributes_m * sizeof (double));
  for (i = 0; i < numattributes_m; ++i) {
    for (j = 0; j < numattributes_m; ++j) {
      if (i != j) {
        probleminfo->preVj[i * numattributes_m + j] = 0.0;
        probleminfo->preVBeta[i * numattributes_m + j] = 0.0;
      }
      else {
        probleminfo->preVj[i * numattributes_m + j] = 1.0;
        probleminfo->preVBeta[i * numattributes_m + j] = 1.0;
      }
    }
  }
  probleminfo->temppreVj   = (double *) malloc ((long unsigned int)numattributes_m * (long unsigned int)numattributes_m * sizeof (double));
  probleminfo->temppreVBeta = (double *) malloc ((long unsigned int)numattributes_m *(long unsigned int) numattributes_m * sizeof (double));
  probleminfo->tempPC       = (double *) malloc ((long unsigned int)numattributes_m * sizeof (double));
  return 0;
} /* end allocateMemory */
