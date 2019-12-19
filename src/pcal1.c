#include <stdlib.h>
#include <math.h>
#include "type.h"


static int initialize(ENTITYINFOptr entityinfo, PROBLEMINFOptr probleminfo);
static int polarityChk(ENTITYINFOptr entityinfo, PROBLEMINFOptr probleminfo);
static int ChkConvergence(PROBLEMINFOptr probleminfo, ENTITYINFOptr entityinfo);
static int dgesvd_pcal1 (char jobu, char jobvt, int m, int n, double *A, int lda, double *S, double *VT, int ldu, double *Umat, int ldvt, double *work, int lwork); /* SVD */

int solvePcaL1(ENTITYINFOptr entityinfo, PROBLEMINFOptr probleminfo)
{
 
 int numattributes_m = entityinfo->numattributes_m;
 int numentities_n   = entityinfo->numentities_n;

 int status = probleminfo->status;
 int q = probleminfo->q;
 int i = probleminfo->i;
 int j = probleminfo->j;
 int k = probleminfo->k;
 int l = probleminfo->l;
 double innerprod = probleminfo->innerprod;
 int initMethod = probleminfo->initMethod;
 /*FILE *projFile;*/

 int position;

 for (k = 0; k < q; ++k) {
   REprintf("%d ", k+1);
   if (k!=0) { /*compute new points_XT for Xj*/ 
    
     /*projFile = fopen("projPointsl1pca.txt", "a");
    
     for (l = 0; l < numattributes_m; ++l) {
       fprintf(projFile, "%f ", probleminfo->wT[l]);
     }
     fprintf(projFile, "\n");
     for (i = 0; i < numentities_n ; ++i) {
       for (l = 0; l < numattributes_m; ++l) {
         fprintf(projFile, "%f ", entityinfo->points_XT[numattributes_m*i+l]);
       }
       fprintf(projFile, "\n");
     }*/
     for (i = 0; i < numentities_n ; ++i) {
       innerprod = 0.0;
       
       for (l = 0; l < numattributes_m; ++l) {
         innerprod += probleminfo->wT[l] * entityinfo->points_XT[numattributes_m*i+l];
       }
       for (j = 0; j < numattributes_m; ++j) {
         position = numattributes_m*i + j;
         entityinfo->points_XT[position] = entityinfo->points_XT[position] - probleminfo->wT[j] * innerprod;
         /*fprintf(projFile,"%f ", entityinfo->points_XT[position]);*/
       }
       /*fprintf(projFile,"%f ", innerprod);
       fprintf(projFile,"\n");*/
     }
     /*fflush(projFile);
     fclose(projFile);*/
   }
   if ((initMethod < 3) || (k != 0)) { /* if user-supplied vector (initMethod = 3), use that on first iteration */
     status = initialize(entityinfo, probleminfo);
     if (status) {
       REprintf("Unable to initialize \n");
       return 1;
     }
   }
    
   probleminfo->convergent = 0;
   if (k < (numattributes_m - 1)) {
     while (probleminfo->convergent == 0) {
         status = polarityChk(entityinfo, probleminfo);
         if(status) {
           REprintf("Polarity Check failed \n");
           return 1;
         }
         status = ChkConvergence(probleminfo, entityinfo);
         if(status) {
           REprintf("Convergence Check failed\n");
           return 1;
         }
     }
   }
   for (j = 0; j < numattributes_m ; ++j) {
     probleminfo->PCs[numattributes_m * k + j] = probleminfo->wT[j];/*save PCs*/
   }
   /*projFile = fopen("projPointsl1pca.txt", "a");
   
   for (l = 0; l < numattributes_m; ++l) {
     fprintf(projFile, "%f ", probleminfo->wT[l]);
   }
   fprintf(projFile, "\n");
   fclose(projFile);*/
 }
 return 0;
} /* end solvePcaL1 */
static int initialize(ENTITYINFOptr entityinfo, PROBLEMINFOptr probleminfo) {
  int numentities_n   = entityinfo->numentities_n;
  int numattributes_m = entityinfo->numattributes_m;
  double *points_XT    = entityinfo->points_XT;
  int initMethod   = probleminfo->initMethod;

  int i = probleminfo->i;
  int j = probleminfo->j;
  
  double *work = probleminfo->work;
  int lwork = probleminfo->lwork;
  double x = probleminfo->x;
  double xSum = probleminfo->xSum;
  double maxSum = probleminfo->maxSum;
  int argmax = probleminfo->argmax;
  double *S = probleminfo->S;
  int status = probleminfo->status;
  double *points_XT_temp = probleminfo->points_XT_temp;
  double normalizer = probleminfo->normalizer;


  if (initMethod == 0) { /* use L2 PCA first component */
    /* get first PC */
    for (i = 0; i < numentities_n; ++i) {
      for (j = 0; j < numattributes_m; ++j) {
        points_XT_temp[numattributes_m*i+j] = points_XT[numattributes_m*i+j];
      }
    }
   
    status = dgesvd_pcal1 ('O', 'N', numattributes_m, numentities_n, points_XT, numattributes_m, S, NULL, 1, NULL, 1, work, lwork);  /* the pointer to points_XT gets destroyed here */
    if (status) {
      REprintf("dgesvd failed\n");
      return 1;
    }

    for (j = 0; j < numattributes_m; ++j) {
      probleminfo->wT[j] = points_XT[j];
    }

    for (i = 0; i < numentities_n; ++i) {
      for (j = 0; j < numattributes_m; ++j) {
        entityinfo->points_XT[numattributes_m *i + j] = points_XT_temp[numattributes_m*i+j];
      }
    }
    return 0; /* need to add options */
  }
  else if (initMethod == 1) { /* find entity with largest norm */
    maxSum=0.0;
 
    for (i=0; i < numentities_n; ++i) {
      xSum = 0.0;
      for (j=0; j<numattributes_m; ++j) {
        x = points_XT[numattributes_m * i + j];/*get the point*/
        xSum=xSum+x*x;/*create sum*/
      }
      xSum=sqrt(xSum);
      if(maxSum < xSum) {/*compare to the max*/
        maxSum=xSum;/*save the max sum*/
        argmax = i;/*save the observations*/
      }
    }
    for (j = 0; j < numattributes_m; ++j) {
      probleminfo->wT[j]=points_XT[numattributes_m * argmax + j]/maxSum;
    } 
    return 0;
  }
  else if (initMethod == 2) { /* random vector */
    xSum = 0.0;
    for (j = 0; j < numattributes_m; ++j) {
      probleminfo->wT[j]= unif_rand() - 0.5;
      xSum += probleminfo->wT[j]*probleminfo->wT[j];
    }
    normalizer = sqrt(xSum);
    for (j = 0; j < numattributes_m; ++j) {
      probleminfo->wT[j]=probleminfo->wT[j]/normalizer;
    }
  
    return 0;
  } 
  return 0;
}
static int polarityChk(ENTITYINFOptr entityinfo, PROBLEMINFOptr probleminfo) {
  double *polarity = probleminfo->polarity;
  int i         = probleminfo->i;
  int j         = probleminfo->j;
  double wTSum = probleminfo->wTSum;
  double dotProd = probleminfo->dotProd;
  double normalizer = probleminfo->normalizer;
  double normx = probleminfo->normx;

  int numattributes_m = entityinfo->numattributes_m;
  int numentities_n   = entityinfo->numentities_n;
  double *points_XT    = entityinfo->points_XT;

  dotProd   = 0.0;
  wTSum = 0.0;

  probleminfo->dotConv=0;

  for (i = 0; i < numentities_n ; ++i) {
    dotProd = 0.0;
    normx = 0.0;
    for (j = 0; j < numattributes_m; ++j) {
      dotProd += probleminfo->wT[j]*points_XT[numattributes_m * i + j];/*obtaining dot product of wT and ith observation */
      normx += points_XT[numattributes_m*i+j]*points_XT[numattributes_m * i + j];/*obtaining norm of XT */
    }
    if (dotProd < 0.0)
      polarity[i] = -1.0;
    else {
      if ((dotProd < (EPSILON)) & (normx > (EPSILON))) {
        probleminfo->dotConv=1;/*saving for convergence check*/
      }
      polarity[i] = 1.0;
    }
  }
  for(i = 0; i < numattributes_m ; ++i) { 
    probleminfo->wTOld[i] = probleminfo->wT[i];/*saving w(T-1)*/
    probleminfo->wT[i] = 0.0;
  }
  
  for (i = 0; i < numattributes_m; ++i) {
    for (j = 0; j < numentities_n ; ++j) {
     probleminfo->wT[i] += polarity[j]*points_XT[numattributes_m * j + i];/*new wT*/
    } 
  }
  for (j = 0; j < numattributes_m; ++j) {
    wTSum = wTSum + probleminfo->wT[j] * probleminfo->wT[j];
  }
  normalizer = sqrt(wTSum);
  
  for (j = 0; j < numattributes_m; ++j) {
    probleminfo->wT[j] = probleminfo->wT[j]/normalizer;
  }
  return 0;
}

static int ChkConvergence(PROBLEMINFOptr probleminfo, ENTITYINFOptr entityinfo) {
  double *wTOld = probleminfo->wTOld;
  int j = probleminfo->j;
  int dotConv = probleminfo->dotConv;

  double wTSum = probleminfo->wTSum;
  double Normalizer = probleminfo->Normalizer;

  int numattributes_m = entityinfo->numattributes_m;

  Normalizer = 0.0;
  wTSum = 0.0;
  
  probleminfo->convergent = 1;

  for (j = 0; j < numattributes_m ; ++j) {
    if (probleminfo->wT[j] - wTOld[j] > (EPSILON)) {
      probleminfo->convergent = 0;
      return 0;
    }
  }
  if (dotConv == 1) {/*if wT(Xi)=0, then add a small non-zero random vector and do Polarity Check again*/
    probleminfo->convergent = 0;
    for(j = 0; j < numattributes_m; ++j){
      probleminfo->wT[j] = probleminfo->wT[j] + (unif_rand() - 0.5)/100.0;
      REprintf("%f\n", probleminfo->wT[j]);
      wTSum = wTSum + probleminfo->wT[j] * probleminfo->wT[j];
    }
    Normalizer = sqrt (wTSum);
    for (j = 0; j < numattributes_m ; ++j){
      probleminfo->wT[j] = probleminfo->wT[j]/Normalizer;
    }
    dotConv = 0;/*reset dotConv for next iteration*/
  }
  
  return 0;
}

static int dgesvd_pcal1 (char jobu, char jobvt, int m, int n, double *A, int lda, double *S, double *VT, int ldu, double *Umat, int ldvt, double *work, int lwork) { /* SVD */
  extern void dgesvd_(const char *jobup, const char *jobvtp, const int *mp, const int *np, double *A, int *ldap, double *S, double *U, const int *ldup, double *Umat, int *ldvtp, double *work, int *lworkp, int *infop);
 
  int info;
  dgesvd_ (&jobu, &jobvt, &m, &n, A, &lda, S, VT, &ldu, Umat, &ldvt, work, &lwork, &info);
  return info;
} /* end dgesvd, SVD */


