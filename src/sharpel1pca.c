#include "type.h"

static int cmp(const void *x, const void *y);
int solveSharpeL1PCA (ENTITYINFOptr entityinfo, SOLVERINFOptr solverinfo, PROBLEMINFOptr probleminfo);

int solveSharpeL1PCA (ENTITYINFOptr entityinfo, SOLVERINFOptr solverinfo, PROBLEMINFOptr probleminfo) {
  int numattributes_m = entityinfo->numattributes_m;
  int numentities_n   = entityinfo->numentities_n;
 
  int q = probleminfo->q;
  int i = probleminfo->i;
  int j = probleminfo->j;
  int k = probleminfo->k;
  int l = probleminfo->l;
  int l1 = probleminfo->l1;
  int l2 = probleminfo->l2;

  double minobjective = probleminfo->minobjective;
  double objective    = probleminfo->objective;
  double *ratios      = probleminfo->ratios;
  double **tosort     = probleminfo->tosort;
  double *weights     = probleminfo->weights;
  double medWeight    = probleminfo->medWeight;
  double sumWeights   = probleminfo->sumWeights;
  double *v           = probleminfo->v;
  double innerprod    = probleminfo->innerprod;
  double normv        = probleminfo->normv;
  int    origIndex    = probleminfo->origIndex;
  int    lstar        = probleminfo->lstar;
  /*int    getScores    = probleminfo->getScores;*/
  FILE *projFile; 

  int index = probleminfo->index;

  for (k = 0; k < q; ++k) { 
    minobjective = (OBJ_INIT); /* objective for the best coordinate to fix */
    for (l = 0; l < numattributes_m; ++l) {
      if (VERBOSITY > 0) {
        REprintf ("l %d ", l);
      }

      /* find v */
      for (j = 0; j < numattributes_m; ++j) {
	if (j != l) {
	  sumWeights = 0.0;
	  index = 0;
          for (i = 0; i < numentities_n; ++i) {
	    if (entityinfo->points_XT[numattributes_m*i+l] != 0.0) {
	      ratios[index] = entityinfo->points_XT[numattributes_m*i+j]/entityinfo->points_XT[numattributes_m*i+l]; /* store ratios */
              tosort[index]=&(ratios[index]); /* sort the pointers to the ratios */
	      sumWeights += fabs(entityinfo->points_XT[numattributes_m*i+l]);
	      weights[index] = fabs(entityinfo->points_XT[numattributes_m*i+l]);
	      index += 1;
	    }
	  }
	  /* get weighted median */
	  qsort(tosort, index,sizeof(double *),cmp); 

	  medWeight = 0.0;
	  for (i = 0; i < index; ++i) {
	    origIndex = tosort[i] - ratios;
	    medWeight += weights[origIndex];
	    if (medWeight > 0.5*sumWeights) {
	      v[j] = *tosort[i];
	      i = index;
	    }
	  }
	}
	else {
	  v[j] = 1.0;
	}
      }

      /* get objective function value */
      objective = 0.0;
      for (i = 0; i < numentities_n; ++i) {
	for (j = 0; j < numattributes_m; ++j) {
	  if (j != l) {
	    objective += fabs(entityinfo->points_XT[numattributes_m*i+j] - entityinfo->points_XT[numattributes_m*i+l]*v[j]);
	  }
	}
      }
      if (VERBOSITY > 1) {
	REprintf ("objective %f\n", objective);
      }

      /* check if best */
      if (objective < minobjective) {
        minobjective = objective;
       
        lstar = l; 

	normv = 0.0;
        for (j = 0; j < numattributes_m; ++j) {
	  probleminfo->PCs[numattributes_m*k+j] = v[j];
	  normv += v[j]*v[j];
        }
      }
    }
    probleminfo->objectives[k] = minobjective;
    if (VERBOSITY > 1) {
      REprintf("k %d lstar %d minobjective %f\n", k, lstar, minobjective);
    }

    /* normalize v */
    for (j = 0; j < numattributes_m; ++j) {
      probleminfo->PCs[numattributes_m*k+j] = probleminfo->PCs[numattributes_m*k+j]/sqrt(normv);
      v[j] = probleminfo->PCs[numattributes_m*k+j]; /* for use in subtracting out below */
    }

    /* project out elements of previous vectors */
    for (j = 0; j < numattributes_m; ++j) {
      for (l2 = 0; l2 < k; ++l2) {
        for (l1 = 0; l1 < numattributes_m; ++l1) {
            probleminfo->PCs[numattributes_m*k+j] -= probleminfo->PCs[numattributes_m*l2+j]*probleminfo->PCs[numattributes_m*l2+l1]*v[l1];
        }
      }
    }
    /* renormalize */
    normv = 0.0;
    for (j = 0; j < numattributes_m; ++j) {
      normv += probleminfo->PCs[numattributes_m*k+j]*probleminfo->PCs[numattributes_m*k+j];
    }
    for (j = 0; j < numattributes_m; ++j) {
      probleminfo->PCs[numattributes_m*k+j] = probleminfo->PCs[numattributes_m*k+j]/sqrt(normv);
    }
    
    /* get scores */
    /*if (getScores == 1) {
      if (VERBOSITY > 2) {
        REprintf("getting scores");
      }
      for (i = 0; i < numentities_n; ++i) {
        probleminfo->scores[numentities_n*k + i] =  entityinfo->points_XT[numattributes_m*i + lstar]/probleminfo->PCs[numattributes_m*k+lstar];
      }
    }

    if (VERBOSITY > 2) {
      REprintf("project data");
    }*/

    /* project data into orthogonal space */
    if (VERBOSITY > 2) {
      projFile = fopen("projPoints.txt", "a");
      for (l = 0; l < numattributes_m; ++l) {
        fprintf(projFile, "%f ", probleminfo->PCs[numattributes_m*k+l]);
      }
      for (i = 0; i < numentities_n ; ++i) {
        for (l = 0; l < numattributes_m; ++l) {
          fprintf(projFile, "%f ", entityinfo->points_XT[numattributes_m*i+l]);
        }
        fprintf(projFile, "\n");
      }
      fflush(projFile);
    }
    for (i = 0; i < numentities_n ; ++i) {
      innerprod = 0.0;
      for (l = 0; l < numattributes_m; ++l) {
        innerprod += probleminfo->PCs[numattributes_m*k+l] * entityinfo->points_XT[numattributes_m*i+l];
      }
      for (j = 0; j < numattributes_m; ++j) {
        entityinfo->points_XT[numattributes_m*i+j] = entityinfo->points_XT[numattributes_m*i+j]-probleminfo->PCs[numattributes_m*k+j]*innerprod;
    
        if (VERBOSITY > 2) {
          fprintf(projFile, "%f ", entityinfo->points_XT[numattributes_m*i+j]);
        }
      }
      if (VERBOSITY > 2) {
        fprintf(projFile, "%f ", innerprod);
        fprintf(projFile, "\n");
      }
    }
    if(VERBOSITY > 2) {
      fflush(projFile);
      fclose(projFile);
    } 
  }

  return 0;
} /*end solveproblem */


static int cmp(const void *x, const void *y) {
  const double **xx = (const double **)x;
  const double **yy = (const double **)y;
  if (**xx < **yy) return -1;
  if (**xx > **yy) return 1;
  return 0;

  /*double xx = *(double*)x, yy=*(double*)y;
  if (xx < yy) return -1;
  if (xx > yy) return 1;
  return 0;*/
}


