/*#include<config.h>*/
#include<stdlib.h>
#include<stdio.h>
#include<string.h>
#include<Rmath.h>
#include <time.h>
#include "type.h"

static int allocateMemoryPcaLp (ENTITYINFOptr entityinfo, PROBLEMINFOptr probleminfo);
int solvePcaLp(ENTITYINFOptr entityinfo, PROBLEMINFOptr probleminfo);

static void
  free_and_null (char **ptr);
void pcalp (double *points_XT, int *dataDim, int *q, double *p, double *PCs, int *initMethod, int *solMethod, double *initV, double *epsilon, double *lratio)
{
  ENTITYINFO entityinfo;
  PROBLEMINFO probleminfo;
  probleminfo.status = 0;
  int status = probleminfo.status;

  probleminfo.polarity = NULL;
  probleminfo.wT       = NULL;
  probleminfo.wTOld    = NULL;  
  probleminfo.wT       = NULL;
  probleminfo.work     = NULL;
  probleminfo.S        = NULL;
  probleminfo.points_XT_temp = NULL;

  probleminfo.PCs      = PCs;
  entityinfo.numentities_n = (int) dataDim[1];
  entityinfo.numattributes_m = (int) dataDim[0];

  probleminfo.q = *q;/*desired number of PCs*/
  probleminfo.epsilon = *epsilon;
  probleminfo.lratio = *lratio;
  probleminfo.p = *p;
  probleminfo.initMethod = *initMethod;
  probleminfo.solMethod = *solMethod;
  probleminfo.wT = initV;
  
  entityinfo.points_XT = points_XT; /* transpose of data matrix */

  GetRNGstate();

  status = allocateMemoryPcaLp (&entityinfo, &probleminfo);
  if (status) {
    REprintf("Unable to allocate memory\n");
    goto TERMINATE;
  }

  
  status = solvePcaLp(&entityinfo, &probleminfo);
  if(status) {
    REprintf("Unable to solve. Terminating...; or done\n");
    goto TERMINATE;
  }

  PutRNGstate();

  REprintf("\n");
TERMINATE:

  free_and_null ((char **) &probleminfo.polarity);
  free_and_null ((char **) &probleminfo.wT);
  free_and_null ((char **) &probleminfo.wTOld);
  free_and_null ((char **) &probleminfo.work);
  free_and_null ((char **) &probleminfo.S);
  free_and_null ((char **) &probleminfo.points_XT_temp);
}

static void
free_and_null (char **ptr)
{
  if( *ptr != NULL){
    free (*ptr);
    *ptr = NULL;
  }
}

static int allocateMemoryPcaLp (ENTITYINFOptr entityinfo, PROBLEMINFOptr probleminfo) {
  int numentities_n   = entityinfo->numentities_n;
  int numattributes_m = entityinfo->numattributes_m;

  probleminfo->polarity = (double *) malloc ((long unsigned int) numentities_n * sizeof(double));
  probleminfo->wT       = (double *) malloc ((long unsigned int) numattributes_m * sizeof(double));
  probleminfo->wTOld    = (double *) malloc ((long unsigned int) numattributes_m * sizeof(double));
  probleminfo->Dw       = (double *) malloc ((long unsigned int) numattributes_m * sizeof(double));
  probleminfo->lwork = 9*(numattributes_m + numentities_n) * (NBMAX);
  probleminfo->work = (double *) malloc((long unsigned int) probleminfo->lwork * sizeof(double));
  probleminfo->S = (double *) malloc ((long unsigned int) numattributes_m * sizeof(double));
  probleminfo->points_XT_temp = (double *) malloc ((long unsigned int) numentities_n * (long unsigned int) numattributes_m * sizeof(double));
/*  probleminfo->PCs      = (double *) malloc (probleminfo->q * numentities_n * sizeof(double));*/
  
  return 0;
}
