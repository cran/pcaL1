#include "Clp_C_Interface.h"
#include <stdlib.h>
#include <math.h>
#include "type.h"

int solveL1Projection (ENTITYINFOptr entityinfo, SOLVERINFOptr solverinfo, PROBLEMINFOptr probleminfo);

static int setupCLPProj (SOLVERINFOptr solverinfo); 
static int loadClpProblemProj (ENTITYINFOptr entityinfo, SOLVERINFOptr solverinfo, PROBLEMINFOptr probleminfo);
static int optimizeProj (SOLVERINFOptr solverinfo, PROBLEMINFOptr probleminfo); 
static int getAlphas (ENTITYINFOptr entityinfo, SOLVERINFOptr solverinfo, PROBLEMINFOptr probleminfo); 
static int getProjectedPointsProj (ENTITYINFOptr entityinfo, PROBLEMINFOptr probleminfo);
static void dgemmProj (char transa, char transb, int m, int n, int k, double alpha, double *A, int lda, double *B, int ldb, double beta, double *C, int ldc); /* multiply A *B = C */

int solveL1Projection (ENTITYINFOptr entityinfo, SOLVERINFOptr solverinfo, PROBLEMINFOptr probleminfo) {
  int status = probleminfo->status;
  int numentities_n     = entityinfo->numentities_n;

  status = setupCLPProj (solverinfo);  /*initialize CLP */
  if (status) {
    REprintf ("Error setting up CLP\n");
    return 1;
  }
 
  for (probleminfo->i=0; probleminfo->i < numentities_n; ++probleminfo->i) {
    if ((VERBOSITY) > 0) {
      REprintf("%d ", probleminfo->i);
    }
      
    status = loadClpProblemProj(entityinfo, solverinfo, probleminfo); /*set up columns/variables */
    if (status) {
      REprintf ("Error with initial load of problem\n");
      return 1;
    }
 
    status = optimizeProj(solverinfo, probleminfo); /* solve with Clp, get objective value */
    if (status) {
      REprintf ("Error solving l1 projection problem\n");
      return 1;
    }
    
    status = getAlphas(entityinfo, solverinfo, probleminfo); /* store scores */
    if (status) {
      REprintf ("Error storing alphas\n");
      return 1;
    }
  }

  status = getProjectedPointsProj (entityinfo, probleminfo);/* get projected points */
  if (status) {
    REprintf ("Error getting projected points\n");
    return 1;
  }

  return 0;
} /*end solveL1Projection */

static int setupCLPProj (SOLVERINFOptr solverinfo) {

  solverinfo->model=Clp_newModel();

  Clp_setLogLevel(solverinfo->model, 0);
  return 0;
} /* end setupCLP */

static int loadClpProblemProj (ENTITYINFOptr entityinfo, SOLVERINFOptr solverinfo, PROBLEMINFOptr probleminfo) {
  int numcols = probleminfo->numcols;
  int i       = probleminfo->i;
  int j       = probleminfo->j;
  int k       = probleminfo->k;
  double *obj = probleminfo->obj;
  /*char   **colname = probleminfo->colname;*/
  int    rcnt    = probleminfo->rcnt;
  int    nzcnt   = probleminfo->nzcnt;
  double *rhs    = probleminfo->rhs;
  int    *matind = probleminfo->matind;
  double *matval = probleminfo->matval;
  int    *matbeg = probleminfo->matbeg;
  int    projdim = probleminfo->projdim;
  double    *lb  = probleminfo->lb;
  double    *ub  = probleminfo->ub;

  double *points_XT   = entityinfo->points_XT;
  int numattributes_m = entityinfo->numattributes_m;
  double *PCs         = entityinfo->PCs;

  rcnt = 0;
  for (j = 0; j < numattributes_m; ++j) {
    rhs[rcnt] = points_XT[i*numattributes_m + j];
    ++rcnt;
  }
  

  nzcnt   = 0;
  numcols = 0;
 
  for (k = 0; k < projdim; ++k) {
    matbeg[numcols] = nzcnt;
    probleminfo->aind[k] = numcols;
    obj[numcols] = 0.0;
    lb[numcols] = -(DBL_MAX);
    ub[numcols] = DBL_MAX;
    /*sprintf (colname[numcols], "alpha_%d_%d",i,k);*/
    for (j = 0; j < numattributes_m; ++j) {
      matind[nzcnt] = j;
      matval[nzcnt] = PCs[k*numattributes_m + j];  
      ++nzcnt;
    }
    ++numcols;
  }
  
      
  for (j = 0; j < numattributes_m; ++j) {
    matbeg[numcols] = nzcnt;
   /* probleminfo->lplus[i][j] = numcols;*/
    obj[numcols] = 1.0;
    lb[numcols] = 0.0;
    ub[numcols] = DBL_MAX;
    /*sprintf (colname[numcols], "lplus_%d_%d",i,j);*/
    matind[nzcnt] = j;
    matval[nzcnt] = 1.0;
    ++nzcnt;
    ++numcols;
  }
  
  for (j = 0; j < numattributes_m; ++j) {
    matbeg[numcols] = nzcnt;
    /*probleminfo->lminus[i][j] = numcols;*/
    obj[numcols] = 1.0;
    lb[numcols] = 0.0;
    ub[numcols] = DBL_MAX;
    /*sprintf (colname[numcols], "lminus_%d_%d",i,j);*/
    matind[nzcnt] = j;
    matval[nzcnt] = -1.0;
    ++nzcnt;
    ++numcols;
  }

  matbeg[numcols] = nzcnt;

  Clp_loadProblem(solverinfo->model, numcols, rcnt, matbeg, matind, matval, lb, ub, obj, rhs, rhs);
  
  return 0;
} /* end loadClpProblem */


static int optimizeProj(SOLVERINFOptr solverinfo, PROBLEMINFOptr probleminfo) {

  int status  = probleminfo->status;
  int solstat = probleminfo->solstat;

  status=Clp_dual(solverinfo->model, 0);
  if (status) {
    /*REprintf ("error solving dual simplex\n");
    return 1;*/
  }

  solstat=Clp_status(solverinfo->model);
  if (solstat != 0) {
    /*REprintf ("solstat %d\n", solstat);
    return 1;*/
  }

  probleminfo->objective = Clp_getObjValue(solverinfo->model);
  if ((VERBOSITY) >= 4) {
    REprintf ("objective value %f\n", probleminfo->objective);
  }

  return 0;
} /* end optimize */

static int getAlphas (ENTITYINFOptr entityinfo, SOLVERINFOptr solverinfo, PROBLEMINFOptr probleminfo) {

  int i       = probleminfo->i;
  int k       = probleminfo->k;
  int projdim = probleminfo->projdim;
  int *aind   = probleminfo->aind;
  const double *projSolution = probleminfo->projSolution;
  
  int numentities_n  = entityinfo->numentities_n;

  projSolution = (const double *) Clp_getColSolution (solverinfo->model);

  for (k = 0; k < projdim; ++k) {
    probleminfo->alphas[numentities_n*k+i] = projSolution[aind[k]];
  }
  

  if ((VERBOSITY) >= 4) {
    for (k = 0; k < projdim; ++k) {
      REprintf("%f ", probleminfo->alphas[numentities_n*k+i]);
    }
    REprintf("\n");
  }
  return 0;
} /* end getAlphas */

static int getProjectedPointsProj(ENTITYINFOptr entityinfo, PROBLEMINFOptr probleminfo) {
  
  int numattributes_m = entityinfo->numattributes_m;
  int numentities_n   = entityinfo->numentities_n;
  double *PCs         = entityinfo->PCs;

  int   projdim       = probleminfo->projdim;
  double *alphas      = probleminfo->alphas;

  dgemmProj('N', 'T', numattributes_m, numentities_n, projdim, 1.0, PCs, numattributes_m, alphas, numentities_n, 0.0, probleminfo->projPoints, numattributes_m);/*get projected points in original coordinates by multiplying preVj by xpluslambda_Z, projected points are columns */

  return 0;
} /* end getProjectedPoints */


static void dgemmProj(char transa, char transb, int m, int n, int k, double alpha, double *A, int lda, double *B, int ldb, double beta, double *C, int ldc) { /* multiply A *B = C */
  extern void dgemm_ (const char *transap, const char *transbp, const int *mp, const int *np, const int *kp, double *alphap, double *A, const int *ldap, double *B, const int *ldbp, const double *betap, double *C, const int *ldcp); 
  dgemm_ (&transa, &transb, &m, &n, &k, &alpha, A, &lda, B, &ldb, &beta, C, &ldc);
} /* end dgemm, multiply A*B */

