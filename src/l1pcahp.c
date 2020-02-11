#include "type.h"

static int setupCLPforL1PCAHp (SOLVERINFOptr solverinfo); 

int solveL1PCAHp (ENTITYINFOptr entityinfo, SOLVERINFOptr solverinfo, PROBLEMINFOptr probleminfo);

int solveL1PCAHp (ENTITYINFOptr entityinfo, SOLVERINFOptr solverinfo, PROBLEMINFOptr probleminfo) {

  int numentities_n     = entityinfo->numentities_n;
  int numattributes_m   = entityinfo->numattributes_m;
  double *points_XT     = entityinfo->points_XT;

  int status            = probleminfo->status;
  double xx_obj         = probleminfo->xx_obj;
  double x_obj          = probleminfo->x_obj;
  double currObj        = probleminfo->currObj;
  double *point         = probleminfo->point;
  double *initV         = probleminfo->initV;
  double *PCs           = probleminfo->PCs;
  const double *alpha   = probleminfo->alpha;
  int numcols           = probleminfo->numcols;
  int numrows           = probleminfo->numrows;
  int nzcnt             = probleminfo->nzcnt;
  double *obj           = probleminfo->obj;
  double *lb            = probleminfo->lb;
  double *ub            = probleminfo->ub;
  int *matbeg           = probleminfo->matbeg;
  int *matind           = probleminfo->matind;
  double *matval        = probleminfo->matval;
  double *rhsL          = probleminfo->rhsL;
  double *rhsU          = probleminfo->rhsU;
  int i                 = probleminfo->i;
  int j                 = probleminfo->j;
  int k                 = probleminfo->k;
  int l                 = probleminfo->l;

  k = 1;

  status = setupCLPforL1PCAHp (solverinfo);  /*initialize CLP */
  if (status) {
    REprintf ("Error setting up CLP\n");
    return 1;
  }

  numcols = numattributes_m + numentities_n;
  numrows = 2 * numentities_n;
  nzcnt = 2 * numentities_n * (numattributes_m + 1);

  /* Set objective function */
  for (i = 0; i < numentities_n; ++i) {
    obj[i] = 1.0;
  }
  for (; i < numattributes_m + numentities_n; ++i) {
    obj[i] = 0.0;
  }

  /* Set columns lower bound */
  for (i = 0; i < numentities_n; ++i) {
    lb[i] = 0.0;
  }
  for (; i < numcols; ++i) {
    lb[i] = -(DBL_MAX);
  }

  /* Set columns upper bound */
  for (i = 0; i < numentities_n + numattributes_m; ++i) {
    ub[i] = (DBL_MAX);
  }

  /* Set rows lower bound */
  for (i = 0; i < numrows; ++i) {
    rhsL[i] = 0.0;
  }

  /* Set rows upper bound */
  for (i = 0; i < 2 * numentities_n; ++i) {
    rhsU[i] = (DBL_MAX);
  }
  for (i = 2*numentities_n; i < numrows; ++i) {
    rhsU[i] = 0.0;
  }

  /* Set matbeg */
  j = 0;
  for (i = 0; i < numentities_n; ++i) {
    matbeg[i] = j;
    j += 2;
  }
  for (i = numentities_n; i < numcols; ++i) {
    matbeg[i] = j;
    j += numrows;
  }
  matbeg[numcols] = nzcnt;

  /* Set matind */
  j = 0;
  for (i = 0; i < numentities_n * 2; ++i) {
    matind[i] = j;
    matind[++i] = j++ + numentities_n;
  }
  for (l = 0; l < numattributes_m; l++)
  {
    for (j = 0; j < numrows; j++)
      matind[i++] = j;
  }

  /* Set matval */
  j = 0;
  for (i = 0; i < numentities_n * 2; ++i) {
    matval[i] = 1;
  }
  for (l = 0;l < numattributes_m;l++)
  {
    for (j = 0; j < numentities_n; j++)
    {
      matval[i++] = points_XT[l +j*numattributes_m];
    }
    for (j = 0; j < numentities_n; j++)
    {
      matval[i++] = - points_XT[l +j*numattributes_m];
    }
  }


  Clp_loadProblem(solverinfo->model, numcols, numrows, matbeg, matind, matval, lb, ub, obj, rhsL, rhsU);
  Clp_setObjSense(solverinfo->model,1); 


  while (k <= numattributes_m)
  {
    REprintf ("%d ",numattributes_m + 1 - k);
    xx_obj = (DBL_MAX);
    x_obj = xx_obj/2;
    currObj = 0;

    for (i = 0; i < numattributes_m; i++)
    {
      point[i] = initV[i+ numattributes_m * (k-1)];
    }
    while ((xx_obj - x_obj > probleminfo->threshold) && (x_obj > probleminfo->threshold))
    {
      matbeg[0] = 0;
      matbeg[1] = numattributes_m;
      rhsL[0] = 1.0;
      for (i = 0; i < numattributes_m; i++)
        matind[i] = numentities_n + i;
      for (i = 0; i < numattributes_m; i++)
        matval[i] = point[i];


      Clp_addRows(solverinfo->model,1,rhsL,rhsL,matbeg,matind,matval);

      Clp_chgObjCoefficients(solverinfo->model, obj);
      Clp_dual(solverinfo->model, 0);

      status=Clp_status(solverinfo->model);
      if (status != 0) {
        REprintf ("solstat %d\n", status);
        /*return 1;*/
      }

      currObj = Clp_getObjValue(solverinfo->model);

      if (currObj > x_obj)
      {
        x_obj = currObj;
      }
      xx_obj = x_obj;
      x_obj = currObj;
      alpha = (const double *) Clp_getColSolution (solverinfo->model);

      //remove the last row
      matbeg[0] = Clp_numberRows(solverinfo->model) - 1;
      int rem = 1;
      Clp_deleteRows(solverinfo->model,rem,matbeg);

      double qr = 0;
      for (i = 0; i < numattributes_m; i++)
      {
        point[i] = alpha[i+numentities_n];
        qr+= point[i]*point[i];
      }
      qr = sqrt(qr);
      for (i = 0; i < numattributes_m; i++)
      {
        point[i] = point[i]/qr;
      }

    } //end while

    matbeg[0] = 0;
    matbeg[1] = numattributes_m;
    rhsL[0] = 0.0;
    for (i = 0; i < numattributes_m; i++)
      matind[i] = numentities_n + i;
    for (i = 0; i < numattributes_m; i++)
      matval[i] = point[i];


    Clp_addRows(solverinfo->model,1,rhsL,rhsL,matbeg,matind,matval);

    for (i = 0; i < numattributes_m; i++)
    {
      PCs[i + numattributes_m *(k-1)] = point[i];
    }
    k++;
  } 
  REprintf("\n");
  return 0;
} /*end solveL1PCAHp */

static int setupCLPforL1PCAHp (SOLVERINFOptr solverinfo) {

  solverinfo->model=Clp_newModel();

  Clp_setLogLevel(solverinfo->model, 0);
  return 0;
} /* end setupCLP */




