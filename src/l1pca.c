#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "type.h"

static int setupCLPforL1PCA (SOLVERINFOptr solverinfo); 
static int solveforU (ENTITYINFOptr entityinfo, SOLVERINFOptr solverinfo, PROBLEMINFOptr probleminfo);
static int solveforV (ENTITYINFOptr entityinfo, SOLVERINFOptr solverinfo, PROBLEMINFOptr probleminfo);
static int normalize(ENTITYINFOptr entityinfo, PROBLEMINFOptr probleminfo);

int solveL1PCA (ENTITYINFOptr entityinfo, SOLVERINFOptr solverinfo, PROBLEMINFOptr probleminfo);

int solveL1PCA (ENTITYINFOptr entityinfo, SOLVERINFOptr solverinfo, PROBLEMINFOptr probleminfo) {
  int    status        = probleminfo->status;
  
  status = setupCLPforL1PCA (solverinfo);  /*initialize CLP */
  if (status) {
    REprintf ("Error setting up CLP\n");
    return 1;
  }
  
  probleminfo->maxdifference = (DBL_MAX);

  probleminfo->iter = probleminfo->iterations;
  while ((probleminfo->iter > 0) && (probleminfo->maxdifference > probleminfo->tolerance)) {
    REprintf("%d ", probleminfo->iterations - probleminfo->iter + 1);
     
    probleminfo->maxdifference = 0.0;
    status = solveforU(entityinfo, solverinfo, probleminfo); /* solve for U, matrix of scores */
    if (status) {
      REprintf ("Error solving for U\n");
      return 1;
    }
      /* solve for V, rotation matrix */  
    status = solveforV(entityinfo, solverinfo, probleminfo); /* solve for V, rotation matrix */
    if (status) {
      REprintf ("Error solving for V\n");
      return 1;
    }
    status = normalize(entityinfo, probleminfo);
    if(status){
      REprintf("Error normalizing\n");
      return 1;
    }
    --probleminfo->iter;
  }
 
  return 0;
} /*end solveproblem */

int setupCLPforL1PCA (SOLVERINFOptr solverinfo) {
  solverinfo->modelU=Clp_newModel(); 
  solverinfo->modelV=Clp_newModel();

  Clp_setLogLevel(solverinfo->modelU, 0);
  Clp_setLogLevel(solverinfo->modelV, 0);
  return 0;
} /* end setupCLP */


int solveforU (ENTITYINFOptr entityinfo, SOLVERINFOptr solverinfo, PROBLEMINFOptr probleminfo) {
  int status       = probleminfo->status;
  int numcols      = probleminfo->numcols;
  int i            = probleminfo->i;
  int j            = probleminfo->j;
  int k            = probleminfo->k;
  int q            = probleminfo->q;
  double *obj      = probleminfo->obj;
  double *lb       = probleminfo->lb;
  double *ub       = probleminfo->ub;
  int    rcnt    = probleminfo->rcnt;
  int    nzcnt   = probleminfo->nzcnt;
  double *rhsL   = probleminfo->rhsL;
  double *rhsU   = probleminfo->rhsU;
  int    *matind = probleminfo->matind;
  double *matval = probleminfo->matval;
  int    *matbeg = probleminfo->matbeg;
  int    *xinda = probleminfo->xinda;
  int    *xindb = probleminfo->xindb;
  double *initV = probleminfo->initV;
  double *V     = probleminfo->V;
  int    solstat = probleminfo->solstat;
  int    *uind  = probleminfo->uind;
  const double *Usol = probleminfo->Usol;
  double *points_XT = entityinfo->points_XT;
  int numentities_n = entityinfo->numentities_n;
  int numattributes_m = entityinfo->numattributes_m;

  /* add rows */
  for (i = 0; i < numentities_n; ++i) {
    rcnt = 0;
    for (j = 0; j < numattributes_m; ++j) {
      xinda[j] = rcnt;
      rhsL[rcnt] = -(DBL_MAX);
      rhsU[rcnt] = points_XT[i*numattributes_m + j];
      ++rcnt;
      rhsL[rcnt] = points_XT[i*numattributes_m + j];
      rhsU[rcnt] = (DBL_MAX);
      xindb[j] = rcnt;
      ++rcnt;
    }

    nzcnt   = 0;
    numcols = 0;
    /* add u_i variables */
    for (k = 0; k < q; ++k) {
      matbeg[numcols] = nzcnt;
      obj[numcols] = 0.0;
      lb[numcols]  = -(DBL_MAX);
      ub[numcols]  = DBL_MAX;
      uind[k] = numcols;
      for (j = 0; j < numattributes_m; ++j) {
        matind[nzcnt] = xinda[j];
        if(probleminfo->iter != probleminfo->iterations){
          matval[nzcnt] = V[numattributes_m*k +j];
          ++nzcnt;
          matind[nzcnt] = xindb[j];
          matval[nzcnt] = V[numattributes_m*k +j];
        }
        else{
          matval[nzcnt] = initV[numattributes_m*k +j];
          ++nzcnt;
          matind[nzcnt] = xindb[j];
          matval[nzcnt] = initV[numattributes_m*k +j];
        }
        ++nzcnt;
      }
      ++numcols;
    }
    /* add delta_j's */
    for (j = 0; j < numattributes_m; ++j) {
      matbeg[numcols] = nzcnt;
      obj[numcols] = 1.0;
      lb[numcols]  = 0.0;
      ub[numcols]  = DBL_MAX;
      matind[nzcnt] = xinda[j];
      matval[nzcnt] = -1.0;
      ++nzcnt;
      matind[nzcnt] = xindb[j];
      matval[nzcnt] = 1.0;
      ++nzcnt;
      ++numcols;
    }
    matbeg[numcols] = nzcnt;
    
    Clp_loadProblem(solverinfo->modelU, (int) numcols, rcnt, matbeg, matind, matval, lb, ub, obj, rhsL, rhsU);

    status = Clp_dual(solverinfo->modelU, 0);
    if (status) {
      REprintf ("error solving dual simplex for U\n");
      return 1;
    }

    solstat=Clp_status(solverinfo->modelU);
    if (solstat != 0) {
      /*REprintf ("solstat %d\n", solstat);
      return 1;*/
    }
    Usol = (const double *) Clp_getColSolution (solverinfo->modelU);
    for (k = 0; k < q; ++k) {
      probleminfo->Utemp[k*numentities_n+i] = Usol[uind[k]];
    }
  }
  return 0;
} /* end solveforU */


int solveforV (ENTITYINFOptr entityinfo, SOLVERINFOptr solverinfo, PROBLEMINFOptr probleminfo) {
  int status       = probleminfo->status;
  int numcols      = probleminfo->numcols;
  int i            = probleminfo->i;
  int j            = probleminfo->j;
  int k            = probleminfo->k;
  int q            = probleminfo->q;
  double *obj      = probleminfo->obj;
  double *lb       = probleminfo->lb;
  double *ub       = probleminfo->ub;
  int    rcnt    = probleminfo->rcnt;
  int    nzcnt   = probleminfo->nzcnt;
  double *rhsL   = probleminfo->rhsL;
  double *rhsU   = probleminfo->rhsU;
  int    *matind = probleminfo->matind;
  double *matval = probleminfo->matval;
  int    *matbeg = probleminfo->matbeg;
  int    **rowind = probleminfo->rowind;
  int    solstat = probleminfo->solstat;
  int    **vind  = probleminfo->vind;
  const double *Vsol = probleminfo->Vsol;
  double *Utemp = probleminfo->Utemp;
  double *points_XT    = entityinfo->points_XT;
  int numentities_n   = entityinfo->numentities_n;
  int numattributes_m = entityinfo->numattributes_m;
 /* add rows */
  rcnt = 0;
  for (i = 0; i < numentities_n; ++i) {
    for (j = 0; j < numattributes_m; ++j) {
      rowind[i][j] = rcnt;
      rhsL[rcnt] = 0.0-points_XT[i*numattributes_m + j];
      rhsU[rcnt] = 0.0-points_XT[i*numattributes_m +j];
      ++rcnt;
    }
  }
  nzcnt = 0;
  numcols = 0;
  /* add lambda variables */
  for (i = 0; i < numentities_n; ++i) {
    for (j = 0; j < numattributes_m; ++j) {
      matbeg[numcols] = nzcnt;
      obj[numcols] = 1.0;
      lb[numcols] = 0.0;
      ub[numcols] = DBL_MAX;
      matind[nzcnt] = rowind[i][j];
      matval[nzcnt] = 1.0;
      ++nzcnt;
      ++numcols;
      matbeg[numcols] = nzcnt;
      obj[numcols] = 1.0;
      lb[numcols] = 0.0;
      ub[numcols] = DBL_MAX;
      matind[nzcnt] = rowind[i][j];
      matval[nzcnt] = -1.0;
      ++nzcnt;
      ++numcols;
    }
  }
  /* add V variables */
  for (j = 0; j < numattributes_m; ++j) {
    for (k = 0; k < q; ++k) {
      matbeg[numcols] = nzcnt;
      vind[j][k] = numcols;
      lb[numcols] = -(DBL_MAX);
      ub[numcols] = DBL_MAX;
      obj[numcols] = 0.0;
      for (i = 0; i < numentities_n; ++i) {
	matind[nzcnt] = rowind[i][j];
        matval[nzcnt] = 0.0-Utemp[k*numentities_n + i];
        ++nzcnt;
      }
      ++numcols;
    }
  }
  matbeg[numcols] = nzcnt;

  Clp_loadProblem(solverinfo->modelV, (int) numcols, rcnt, matbeg, matind, matval, lb, ub, obj, rhsL, rhsU);
  status = Clp_dual(solverinfo->modelV, 0);
  if (status) {
    REprintf ("error solving dual simplex for V\n");
    return 1;
  }

  solstat=Clp_status(solverinfo->modelV);
  if (solstat != 0) {
    REprintf ("solstat %d\n", solstat);
    return 1;
  }
  Vsol = (const double *) Clp_getColSolution (solverinfo->modelV);
  for (j = 0; j < numattributes_m; ++j) {
    for (k = 0; k < q; ++k) {
      probleminfo->Vtemp[numattributes_m*k+j] = Vsol[vind[j][k]];
    }
  }
 
  return 0;
} /* end solveforV */

int normalize(ENTITYINFOptr entityinfo, PROBLEMINFOptr probleminfo){
  
  double *nv    = probleminfo->nv;
  int i         = probleminfo->i;
  int j         = probleminfo->j;
  int k         = probleminfo->k;
  int numattributes_m = entityinfo->numattributes_m;
  int numentities_n   = entityinfo->numentities_n;
  double udiff = probleminfo->udiff;
  int iter     = probleminfo->iter;
  int iterations     = probleminfo->iterations;
  double vdiff  = probleminfo->vdiff;
  double *Vtemp = probleminfo->Vtemp;
  double *Utemp = probleminfo->Utemp;
  int q = probleminfo->q;

  for(i=0; i < q; ++i){/*calculate nv */
    nv[i] = 0.0;
    for(j=0; j < numattributes_m; ++j){
      nv[i] += Vtemp[i*numattributes_m+j] * Vtemp[i*numattributes_m+j];
    }
  }
  
  for(i=0; i < q; ++i){
    for(j=0; j < numattributes_m; ++j){
      if (nv[i] > (EPSILON)) {
        Vtemp[i*numattributes_m+j] = Vtemp[i*numattributes_m+j] * (1/sqrt(nv[i]));
      }
    }
    /* calculating maxdifference after V*/
    vdiff = 0.0;
    for (j = 0; j < numattributes_m; ++j) {
      for (k = 0; k < q; ++k) {
        vdiff = fabs(Vtemp[numattributes_m*k+j] - probleminfo->V[numattributes_m*k+j]);
        if (vdiff > probleminfo->maxdifference) {
          probleminfo->maxdifference = vdiff;
        }
        probleminfo->V[numattributes_m*k+j] = Vtemp[numattributes_m*k+j];
      }
    }
  }

/*calculating mxdifference after U*/
  udiff = 0.0;
  for( i = 0; i < numentities_n; ++i){
    for (k = 0; k < q; ++k) {
      if (iter != iterations) {
	udiff = fabs(Utemp[k*numentities_n+i] - probleminfo->U[k*numentities_n+i]);
        if (udiff > probleminfo->maxdifference) {
          probleminfo->maxdifference = udiff;
        }
      }
      probleminfo->U[k*numentities_n+i] = Utemp[k*numentities_n+i];
    }
  }


  return 0;
}
      



