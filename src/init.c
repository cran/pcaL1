#include "type.h"
#include <stdlib.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Visibility.h>


void l1pcahp (double *points_XT, int *dataDim, double *threshold, double *initV, double *PCs);
void l1pcastar (double *points_XT, int *dataDim, int *q, double *PCs);
void l1pca (double *points_XT, int *dataDim, int *q, double *tolerance, int *iterations, double *initV, double *PCs, double *Scores);
void l1projection (double *points_XT, int *dataDim, int *q, double *PCs, double *projPoints, double *alphas);
void pcal1 (double *points_XT, int *dataDim, int *q, double *PCs, int *initMethod, double *initV);
void pcalp (double *points_XT, int *dataDim, int *q, double *p, double *PCs, int *initMethod, int *solMethod, double *initV, double *epsilon, double *lratio);
void sharpel1pca (double *points_XT, int *dataDim, int *q, double *PCs, double *objectives);

#define C_DEF(name, n) {#name, (DL_FUNC) &name, n}

static const R_CMethodDef cMethods[] = {
  C_DEF(l1pcahp, 5),
  C_DEF(l1pcastar, 4),
  C_DEF(l1pca, 8),
  C_DEF(l1projection, 6),
  C_DEF(pcal1, 6),
  C_DEF(pcalp, 10),
  C_DEF(sharpel1pca, 5),
  {NULL, NULL, 0}
};

void attribute_visible R_init_pcaL1(DllInfo *info)
{
  R_registerRoutines(info, cMethods, NULL, NULL, NULL);
  R_useDynamicSymbols(info, FALSE);
  R_forceSymbols(info, TRUE); 
}


