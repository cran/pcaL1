#include "type.h"

int solveL1PCAStar (ENTITYINFOptr entityinfo, SOLVERINFOptr solverinfo, PROBLEMINFOptr probleminfo);

static int initialSVD (ENTITYINFOptr entityinfo, PROBLEMINFOptr probleminfo);
static int setupCLP (SOLVERINFOptr solverinfo); 
static int loadClpProblem (ENTITYINFOptr entityinfo, SOLVERINFOptr solverinfo, PROBLEMINFOptr probleminfo);
static int changeBounds (SOLVERINFOptr solverinfo, PROBLEMINFOptr probleminfo, int l);
static int optimize (SOLVERINFOptr solverinfo, PROBLEMINFOptr probleminfo); 
static int getBeta (SOLVERINFOptr solverinfo, PROBLEMINFOptr probleminfo); 
static int getProjectedPoints (ENTITYINFOptr entityinfo, PROBLEMINFOptr probleminfo);
static int dgesvd (char jobu, char jobvt, int m, int n, double *A, int lda, double *S, double *VT, int ldu, double *Umat, int ldvt, double *work, int lwork); /* SVD */
static void dgemm (char transa, char transb, int m, int n, int k, double alpha, double *A, int lda, double *B, int ldb, double mybeta, double *C, int ldc); /* multiply A *B = C */
static void dgemv (char trans, int m, int n, double alpha, double *A, int lda, double *x, int incx, double mybeta, double *y, int incy);  /* multiply Ax = y */

int solveL1PCAStar (ENTITYINFOptr entityinfo, SOLVERINFOptr solverinfo, PROBLEMINFOptr probleminfo) {
  int status = probleminfo->status;
  int l      = probleminfo->l;
  double minobjective = probleminfo->minobjective;
  int numprojdim      = probleminfo->numprojdim;
  int numfactors      = probleminfo->numfactors;

  int numentities_n     = entityinfo->numentities_n;
  int numattributes_m   = entityinfo->numattributes_m;
  numprojdim = 0;

  status = setupCLP (solverinfo);  /*initialize CLP */
  if (status) {
    REprintf ("Error setting up CLP\n");
    return 1;
  }

  /* if numattributes_m > numentities_n, project into numentities_n space, then set numattributes_m = numentities_n and continue */
  if (numattributes_m > numentities_n) {
    status = initialSVD (entityinfo, probleminfo);
    if (status) {
      REprintf ("Error with initial SVD\n");
      return 1;
    }
  }
 
  for (probleminfo->projdim = numfactors - 1; probleminfo->projdim > 0; --probleminfo->projdim) {  /* the main loop */
    minobjective = (OBJ_INIT); /* objective of the best l1 regression */
    Clp_setDualObjectiveLimit(solverinfo->model, minobjective);

    status = loadClpProblem(entityinfo, solverinfo, probleminfo); /*set up columns/variables */
    if (status) {
      REprintf ("Error with initial load of problem\n");
      return 1;
    }
 
    REprintf ("%d ", probleminfo->projdim);
    if ((VERBOSITY) >= 1) {
      REprintf ("projdim %d\n", probleminfo->projdim);
    }

    for (l = 0; l <= probleminfo->projdim; ++l) {  /* solve l1 regression for each direction, lth attribute is response */
      if ((VERBOSITY) >= 1) {
        REprintf ("l %d minobjective %f\n", l, minobjective);
      }

      status = changeBounds (solverinfo, probleminfo, l); /* set beta_l = -1, others unbounded  */
      if (status) {
        REprintf ("Error changing beta bounds\n");
        return 1;
      }
   
      status = optimize (solverinfo, probleminfo); /* solve with Clp, get objective value */
      if (status) {
        REprintf ("Error solving l1 regression subproblem\n");
        return 1;
      }
  
      if (minobjective > probleminfo->objective) { /* determine if objective is best so far */
        probleminfo->bestdir[probleminfo->projdim] = l;
        minobjective = probleminfo->objective;
        probleminfo->minobjective = minobjective;
        if ((VERBOSITY) >= 4) {
          REprintf ("new bestdir %d new best objective %f\n", l, probleminfo->objective);
        }
         
        status = getBeta(solverinfo, probleminfo); /* store coefficients of regression */
        if (status) {
          REprintf ("Error storing beta's\n");
          return 1;
        }
        
        Clp_setDualObjectiveLimit(solverinfo->model, minobjective);
      }

      if (minobjective < (EPSILON)) { /* if error is 0, go to next dimension */
        l = probleminfo->projdim + 1;
      }
    } /* end loop on l1 regressions */
    if ((VERBOSITY) >= 1) {
      REprintf ("avg error %f\n", minobjective/((double) entityinfo->numentities_n)); 
    }
    if (minobjective/((double) entityinfo->numentities_n) < 1.0) {
      ++numprojdim;
    }

    if ((VERBOSITY) >= 4) {
      REprintf ("bestdir %d best objective %f\n", probleminfo->bestdir[probleminfo->projdim], minobjective);
    }
    
    status = getProjectedPoints (entityinfo, probleminfo);/* get projected points */
    if (status) {
      REprintf ("Error getting projected points, or done\n");
      return 1;
    }
  } 

  
  return 0;
} /*end solveL1PCAStar */

int setupCLP (SOLVERINFOptr solverinfo) {

  solverinfo->model=Clp_newModel();

  Clp_setLogLevel(solverinfo->model, 0);
  return 0;
} /* end setupCLP */

static int initialSVD (ENTITYINFOptr entityinfo, PROBLEMINFOptr probleminfo) {
  int    i = probleminfo->i;
  int    j = probleminfo->j;
  int    status = probleminfo->status;
  double *work = probleminfo->work;
  int    lwork = probleminfo->lwork;
  double *S    = probleminfo->S;
  double *VT   = probleminfo->VT;
  

  int numattributes_m = entityinfo->numattributes_m;
  int numentities_n   = entityinfo->numentities_n;
  double *points_XT = entityinfo->points_XT;

  if ((VERBOSITY) >= 7) {
    REprintf ("points\n");
    for (i = 0; i < numattributes_m; ++i) {
      for (j = 0; j < numentities_n; ++j) {
        REprintf ("%f ", points_XT[numattributes_m * j + i]);
      }
      REprintf ("\n");
    }
  }

  status = dgesvd ( 'A', 'A', numattributes_m, numentities_n, points_XT, numattributes_m, S, probleminfo->preVj, numattributes_m, VT, numentities_n, work, lwork); /* get A_q = preVj, points gets destroyed here; The columns of preVj define the new subspace */
  if (status) { /* preV^jT is the first projdim columns of U, derived by doing an SVD on points_XT (which has dimension numattributes_m X numentities_n).  It's like doing PCA on the transpose of the data matrix; U*S has the scores */
    REprintf ("Error getting SVD, error %d\n", status);
    return 1;
  }

  for (i = 0; i < numentities_n; ++i) {
    for (j = 0; j < numentities_n; ++j) {
      entityinfo->points_XT[numentities_n * j + i] =  probleminfo->VT[numentities_n * j + i] * S[i];
    }
  }
  
 return 0;

} /* end initial SVD */


static int loadClpProblem (ENTITYINFOptr entityinfo, SOLVERINFOptr solverinfo, PROBLEMINFOptr probleminfo) {
  int numcols = probleminfo->numcols;
  int i       = probleminfo->i;
  int j       = probleminfo->j;
  double *obj = probleminfo->obj;
  char   **colname = probleminfo->colname;
  int    rcnt    = probleminfo->rcnt;
  int    nzcnt   = probleminfo->nzcnt;
  double *rhs    = probleminfo->rhs;
  int    *matind = probleminfo->matind;
  double *matval = probleminfo->matval;
  int    *matbeg = probleminfo->matbeg;
  int    projdim   = probleminfo->projdim;
  double *points_XT = entityinfo->points_XT;

  int numentities_n     = entityinfo->numentities_n;


  rcnt = numentities_n;
  for (i = 0; i < numentities_n; ++i) {
    rhs[i] = 0.0;
  }

  nzcnt   = 0;
  numcols = 0;
  for (j = 0; j < (projdim + 1); ++j) {
    matbeg[numcols] = nzcnt;
    probleminfo->betaind[j] = numcols;
    obj[numcols] = 0.0;
    probleminfo->lb[numcols]=-(DBL_MAX);
    probleminfo->ub[numcols]=DBL_MAX;
    sprintf (colname[numcols], "beta_%d", j);
    for (i = 0; i < numentities_n; ++i) {
      if (points_XT[i*(projdim + 1) + j]!=0.0) {
        matind[nzcnt]=i;
        matval[nzcnt]=points_XT[i * (projdim + 1) + j];
        ++nzcnt;
      }
    }
    ++numcols;
  }
  for (i = 0; i < numentities_n; ++i) {
    matbeg[numcols] = nzcnt;
    probleminfo->eplusind[i] = numcols;
    obj[numcols] = 1.0;
    probleminfo->lb[numcols]  = 0.0;
    probleminfo->ub[numcols]=(DBL_MAX);
    sprintf (colname[numcols], "eplus_%d", i);
    matind[nzcnt] = i;
    matval[nzcnt] = 1.0;
    ++nzcnt;
    ++numcols;
  }
  for (i = 0; i < numentities_n; ++i) {
    matbeg[numcols] = nzcnt;
    probleminfo->eminusind[i] = numcols;
    obj[numcols] = 1.0;
    probleminfo->lb[numcols]  = 0.0;
    probleminfo->ub[numcols]  =(DBL_MAX);/* (CPX_INFBOUND);*/
    sprintf (colname[numcols], "eminus_%d", i);
    matind[nzcnt] = i;
    matval[nzcnt] = -1.0;
    ++nzcnt;
    ++numcols;
  }
  matbeg[numcols] = nzcnt;

  Clp_loadProblem(solverinfo->model, numcols, rcnt, matbeg, matind, matval, probleminfo->lb, probleminfo->ub, obj, rhs, rhs);

  
  return 0;
} /* end loadClpProblem */


static int changeBounds (SOLVERINFOptr solverinfo, PROBLEMINFOptr probleminfo, int l) {
  int j        = probleminfo->j;
  int projdim  = probleminfo->projdim;
  int *betaind = probleminfo->betaind;

  for (j = 0; j < (projdim + 1); ++j) {
    probleminfo->lb[betaind[j]]=-(DBL_MAX);
    probleminfo->ub[betaind[j]]=DBL_MAX;
  }
  probleminfo->lb[betaind[l]] = -1.0;
  probleminfo->ub[betaind[l]] = -1.0;

  Clp_chgColumnLower(solverinfo->model, probleminfo->lb);
  Clp_chgColumnUpper(solverinfo->model, probleminfo->ub);

  /*Clp_writeMps(solverinfo->model, "test.mps");*/
  return 0;
} /* end changeBounds */
      
   
static int optimize (SOLVERINFOptr solverinfo, PROBLEMINFOptr probleminfo) {

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

static int getBeta (SOLVERINFOptr solverinfo, PROBLEMINFOptr probleminfo) {

  int i       = probleminfo->i;
  int numfactors = probleminfo->numfactors;

  probleminfo->currBeta = (const double *) Clp_getColSolution (solverinfo->model);
 
  for (i = 0; i < numfactors; ++i) {
    probleminfo->mybeta[i] = probleminfo->currBeta[i];
  }

  if ((VERBOSITY) >= 4) {
    for (i = 0; i < numfactors; ++i) {
      REprintf("beta[%d] %f\n", i, probleminfo->mybeta[i]);
    }
  }
  return 0;
} /* end getBeta */

static int getProjectedPoints (ENTITYINFOptr entityinfo, PROBLEMINFOptr probleminfo) {
  
  int numattributes_m = entityinfo->numattributes_m;
  int numentities_n   = entityinfo->numentities_n;
  double *points_XT    = entityinfo->points_XT;

  int i = probleminfo->i;
  int j = probleminfo->j;
  int   status         = probleminfo->status;
  int   projdim        = probleminfo->projdim;
  int   *bestdir       = probleminfo->bestdir;
  double *xpluslambda_Z  = probleminfo->xpluslambda_Z;
  double *xpluslambda_Z2 = probleminfo->xpluslambda_Z2;
  const double *mybeta   = probleminfo->mybeta;
  double *work         = probleminfo->work;
  int    lwork         = probleminfo->lwork;
  double    *S         = probleminfo->S;
  double    *VT        = probleminfo->VT;
  double  *a           = probleminfo->a;
  double *temppreVj   = probleminfo->temppreVj;
  double *temppreVBeta = probleminfo->temppreVBeta;
  double *tempPC       = probleminfo->tempPC;
  /*int    q             = probleminfo->q;*/

  if ((VERBOSITY) >= 4) {
    REprintf ("projdim %d bestdir %d \n", projdim, bestdir[projdim]);
    for (j = 0; j <= projdim; ++j) {
      REprintf ("beta %f\n", mybeta[j]);
    }
    for (i = 0; i < numentities_n; ++i) {
      for (j = 0; j <= projdim; ++j) {
        REprintf ("%f ", points_XT[i*(projdim+1)+j]);
      }
      REprintf ("\n");
    }
  }
  
  for (i = 0; i < numentities_n; ++i) {   /* get projected points, in terms of coordinates of current subspace */
    xpluslambda_Z[i * (projdim + 1)+ bestdir[projdim]] = 0.0;
    for (j = 0; j <= projdim; ++j) {
      if (j != bestdir[projdim]) {
        xpluslambda_Z[i * (projdim + 1) + bestdir[projdim]] += mybeta[j] * points_XT[i * (projdim + 1) + j];
        xpluslambda_Z[i * (projdim + 1) + j] = points_XT[i * (projdim + 1) + j];
      }
    }
  }
  for (i = 0; i < numentities_n; ++i) {
    for (j = 0; j <= projdim; ++j) {
      xpluslambda_Z2[i * (projdim + 1) + j] = xpluslambda_Z[i * (projdim + 1) + j];
      if ((VERBOSITY) >= 4) {
        REprintf ("%f ", xpluslambda_Z2[i*(projdim+1)+j]);
      }
    }
    if ((VERBOSITY) >= 4) {
      REprintf ("\n");
    }
  }
  if ((VERBOSITY) < 4) {
    status = dgesvd ( 'A', 'N', projdim + 1, numentities_n, xpluslambda_Z, projdim + 1, S, probleminfo->Vj, projdim + 1, VT, numentities_n, work, lwork); /* get V^jT, xpluslambda_Z gets destroyed here; The columns of Vj define the new subspace */
    if (status) { /* V^jT is the first projdim columns of U, derived by doing an SVD on xpluslambda_Z (which has dimension projdim + 1 X numentities_n).  It's like doing PCA on the transpose of the data matrix */
      REprintf ("Error getting SVD, error %d\n", status);
      return 1;
    }
  }
  else {
    status = dgesvd ( 'A', 'S', projdim + 1, numentities_n, xpluslambda_Z, projdim + 1, S, probleminfo->Vj, projdim + 1, VT, numentities_n, work, lwork); /* get V^jT, xpluslambda_Z gets destroyed here; The columns of Vj define the new subspace */
    if (status) { /* V^jT is the first projdim columns of U, derived by doing an SVD on xpluslambda_Z (which has dimension projdim + 1 X numentities_n).  It's like doing PCA on the transpose of the data matrix */
      REprintf ("Error getting SVD, error %d\n", status);
      return 1;
    }
  }
  
  for (i = 0; i < projdim + 1; ++i) {
    a[i] = probleminfo->Vj[projdim * (projdim + 1) + i];/*saving normalized beta*/
  }
  for (i = 0; i < numattributes_m; ++i) {
    tempPC[i] = 0.0;
  }
  if ((VERBOSITY) >= 4) {
    for (i = 0; i < numattributes_m * (projdim + 1); ++i) {
      REprintf ("preVj[%d] %f\n", i, probleminfo->preVj[i]); /* product of all previous Vj */
    }
  }
  /*if (probleminfo->getProjPoints == 1) {
    if (projdim == q) {*/ /*print projected points in original coordinates*/
      /*for (i = 0; i < numattributes_m * numentities_n; ++i) {
        probleminfo->projPoints[i] = 0.0;
      }
      dgemm('N', 'N', numattributes_m, numentities_n, projdim + 1, 1.0, probleminfo->preVj, numattributes_m, xpluslambda_Z2, projdim + 1, 1.0, probleminfo->projPoints, numattributes_m);*//*get projected points in original coordinates by multiplying preVj by xpluslambda_Z, projected points are columns */
/*    }
  }*/

  dgemv ('N', numattributes_m, projdim + 1, 1.0, probleminfo->preVj, numattributes_m, a, 1, 1.0, tempPC, 1); /*multiply old Vj's by a, to get projdim^th PC */
  for (j = 0; j <numattributes_m; ++j) {
    probleminfo->b[numattributes_m * projdim + j] = tempPC[j];
  }
  /*new VBeta times previous VBeta */
  for (i = 0; i < numattributes_m * (projdim + 1); ++i) {
    temppreVBeta[i] = 0.0;
  }
  for (i = 0; i < projdim; ++i) {
    for (j = 0; j < projdim + 1; ++j) {
      if (j != bestdir[projdim]) {
        probleminfo->VBeta[projdim * j + i] = probleminfo->Vj[(projdim + 1) * i + j] + probleminfo->Vj[(projdim + 1) * i + bestdir[projdim]] * mybeta[j];
      }
      else {
        probleminfo->VBeta[projdim * j + i] = 0.0;
      } 
    }
  }

  dgemm ('N', 'N', projdim, numattributes_m, projdim + 1, 1.0, probleminfo->VBeta, projdim, probleminfo->preVBeta, projdim + 1, 0.0, temppreVBeta, projdim); /* get preVBeta - projdim by numattributes_m  */

  for (i = 0; i < numattributes_m * projdim; ++i) {
    probleminfo->preVBeta[i] = temppreVBeta[i];
  }
  /*previous Vj's times Vj */
  for (i = 0; i < numattributes_m * (projdim + 1); ++i) {
    temppreVj[i] = 0.0;
  }
  dgemm ('N', 'N', numattributes_m, projdim, projdim + 1, 1.0, probleminfo->preVj, numattributes_m, probleminfo->Vj, projdim + 1, 1.0, temppreVj, numattributes_m);/* get preVj - numattributes_m by projdim*/

  for (i = 0; i < numattributes_m * projdim; ++i) {
    probleminfo->preVj[i] = temppreVj[i];
  }
  if (projdim == 1) { /* get 1st PC */
    for (i = 0; i < numattributes_m; ++i) {
      probleminfo->b[i] = probleminfo->preVj[i];
    }
  }

  if ((VERBOSITY) >= 4) {
    REprintf ("V_j^T\n");
    for (i = 0; i < projdim + 1; ++i) {
      for (j = 0; j < projdim; ++j) {
	REprintf ("%f ", probleminfo->Vj[j * (projdim + 1) + i]);
      }
      REprintf ("\n");
    }
  }

  for (i = 0; i < projdim * numentities_n; ++i) {
    entityinfo->points_XT[i] = 0.0;
  }

  dgemm ('T', 'N', projdim, numentities_n, projdim + 1, 1.0, probleminfo->Vj, projdim + 1, xpluslambda_Z2, projdim + 1, 1.0, entityinfo->points_XT, projdim); /* get new points in new coordinates by multiplying Vj (transpose of Vj) by xpluslambda_Z */
  
  if ((VERBOSITY) >= 4) {
    REprintf ("x^j-1\n");
    for (i = 0; i < numentities_n; ++i) {
      for (j = 0; j < projdim; ++j) {
        REprintf ("%f ", entityinfo->points_XT[i * projdim + j]);
      }
      REprintf ("\n");
    }
  }
  /* get scores- projected points in terms of new coordinates */
  /*if ((probleminfo->getScores == 1) && (projdim == q)) {
    for (i = 0; i < numentities_n; ++i) {
      for (j = 0; j < projdim; ++j) {
        probleminfo->scores[i * projdim + j] = entityinfo->points_XT[i * projdim + j];
     }
    }
  }*/
  return 0;
} /* end getProjectedPoints */

static int dgesvd (char jobu, char jobvt, int m, int n, double *A, int lda, double *S, double *VT, int ldu, double *Umat, int ldvt, double *work, int lwork) { /* SVD */
  extern void dgesvd_(const char *jobup, const char *jobvtp, const int *mp, const int *np, double *A, int *ldap, double *S, double *U, const int *ldup, double *Umat, int *ldvtp, double *work, int *lworkp, int *infop);
 
  int info;
  dgesvd_ (&jobu, &jobvt, &m, &n, A, &lda, S, VT, &ldu, Umat, &ldvt, work, &lwork, &info);
  return info;
} /* end dgesvd, SVD */


static void dgemm (char transa, char transb, int m, int n, int k, double alpha, double *A, int lda, double *B, int ldb, double mybeta, double *C, int ldc) { /* multiply A *B = C */
  extern void dgemm_ (const char *transap, const char *transbp, const int *mp, const int *np, const int *kp, double *alphap, double *A, const int *ldap, double *B, const int *ldbp, const double *betap, double *C, const int *ldcp); 
  dgemm_ (&transa, &transb, &m, &n, &k, &alpha, A, &lda, B, &ldb, &mybeta, C, &ldc);
} /* end dgemm, multiply A*B */

static void dgemv (char trans, int m, int n, double alpha, double *A, int lda, double *x, int incx, double mybeta, double *y, int incy) { /* multiply Ax = y */
  extern void dgemv_ (const char *transp, const int *mp, const int *np, double *alphap, double *A, const int *ldap, double *x, const int *incxp, const double *betap, double *y, const int *incyp); 
  dgemv_ (&trans, &m, &n, &alpha, A, &lda, x, &incx, &mybeta, y, &incy);
} /* end dgemv, multiply Ax */

