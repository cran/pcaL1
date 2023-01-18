#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "type.h"

int solveSparsEl(ENTITYINFOptr entityinfo, PROBLEMINFOptr probleminfo);

int sparsEl_pointer_cmp(const void *x, const void *y);
int sparsEl_cmp(const void *x, const void *y);

int solveSparsEl(ENTITYINFOptr entityinfo, PROBLEMINFOptr probleminfo) {

  /*FILE *rotationfile = ioinfo->rotationfile;*/

  int numattributes_m = entityinfo->numattributes_m;
  int numentities_n   = entityinfo->numentities_n;
 
  int q = probleminfo->q;
  int i = probleminfo->i;
  int j = probleminfo->j;
  int k = probleminfo->k;
  int l = probleminfo->l;
  /*int l1 = probleminfo->l1;
  int l2 = probleminfo->l2;*/

  double *ratios      = probleminfo->ratios;
  double **tosort     = probleminfo->tosort;
  double *weights     = probleminfo->weights;
  /*double medWeight    = probleminfo->medWeight;*/
  double sumWeights   = probleminfo->sumWeights;
  double normv        = probleminfo->normv;
  int    origIndex    = probleminfo->origIndex;

  int index = probleminfo->index;

  double sum_below        = probleminfo->sum_below;
  double sum_above        = probleminfo->sum_above;
  double lambdaU          = probleminfo->lambdaU;
  double lambdaL          = probleminfo->lambdaL;

  double *lambdas         = probleminfo->lambdas;
  int    num_lambdas      = probleminfo->num_lambdas;
  double lambda_max       = probleminfo->lambda_max;
  int    **num_lambdas_lj = probleminfo->num_lambdas_lj;
  double ***v_lj          = probleminfo->v_lj;
  double ***lambdas_lj    = probleminfo->lambdas_lj;
  double ****lambdas_lj_sort = probleminfo->lambdas_lj_sort;
  int    **curr_lambda_lj = probleminfo->curr_lambda_lj;
  /*int    t                = probleminfo->t;*/
  int    num_distinct_lambdas = probleminfo->num_distinct_lambdas;
  int    max_memory       = probleminfo->max_memory;
  int    **max_memory_lj  = probleminfo->max_memory_lj;
  int    max_memory_lambdas = probleminfo->max_memory_lambdas;
  double *zs    = probleminfo->zs;
  double *vs    = probleminfo->vs;
  double maxL   = probleminfo->maxL;
  double minU   = probleminfo->minU;
  int lambda_curr = probleminfo->lambda_curr;
  int lambda_next = probleminfo->lambda_next;
  int check_lambda = probleminfo->check_lambda;
  int v_index = probleminfo->v_index;

  double minobjective = probleminfo->minobjective;
  double objective    = probleminfo->objective;
  double *v           = probleminfo->v;
  double innerprod    = probleminfo->innerprod; 
  int    lstar        = probleminfo->lstar;
  FILE *projFile; 
  int l1 = probleminfo->l1;
  int l2 = probleminfo->l2;

  if (probleminfo->lambdas_out[0] >= 0.0) {  /* one lambda specified */
    for (k = 0; k < q; ++k) { 
      minobjective = (OBJ_INIT); /* objective for the best coordinate to preserve */
      if (VERBOSITY > 0) {
        fprintf(stdout, "k %d\n", k);
        fflush(stdout);
      }

      num_lambdas = 1;
      for (l = 0; l < numattributes_m; ++l) { /* for each fixed coordinate l */
        if (VERBOSITY > 0) {
          fprintf (stdout, "l %d ", l);
          fflush (stdout);
        }
        /* find v */
        v[l] = 1.0;
        for (j = 0; j < numattributes_m; ++j) {
          if (j != l) {
            sumWeights = 0.0;
            index = 0;
            for (i = 0; i < numentities_n; ++i) {
              if (entityinfo->points_XT[numattributes_m*i+l] != 0.0) {
                ratios[index] = entityinfo->points_XT[numattributes_m*i+j]/entityinfo->points_XT[numattributes_m*i+l]; /* store ratios */
                tosort[index]=&(ratios[index]); /* sort the pointers to the ratios */
                weights[index] = fabs(entityinfo->points_XT[numattributes_m*i+l]);
                sumWeights += fabs(entityinfo->points_XT[numattributes_m*i+l]);
                index += 1;
              }
            }
            /* sort */
            qsort(tosort,index,sizeof(double *),sparsEl_pointer_cmp); 

            sum_below=0.0; /* keep track of sum below and sum above; updating involves one addition/subtraction each */
            sum_above=sumWeights;
            v[j] = 0.0;
            for (i = 0; i < index; ++i) {
              origIndex = tosort[i] - ratios;
              sum_above -= weights[origIndex];
             
              /* (lambdaL, lambdaU) is the interval where the current ratio is the best for this j/l combination (l is jhat) */
              lambdaU = (sum_above - sum_below) + weights[origIndex];
              if (*tosort[i] < 0.0) { /* here's where we check the sign of the ratio */
                lambdaU = (sum_below - sum_above) + weights[origIndex];
              }
              lambdaL = lambdaU - 2*weights[origIndex];
              if (lambdaU > 0.0) {  /* if lambdaU <= 0.0, then current point is never the best ratio */
                lambdaL = lambdaU - 2*weights[origIndex];
                if ((probleminfo->lambdas_out[0] >= lambdaL) && (probleminfo->lambdas_out[0] <= lambdaU)) { /* check if given lambda is in interval */
                  v[j] = *tosort[i];
                  i=index;
                }
              }
              sum_below += weights[origIndex];
            }
          }
        }
        objective = 0.0;
        for (j = 0; j < numattributes_m; ++j) {
          if (j != l) {
            objective += probleminfo->lambdas_out[0]*fabs(v[j]);
            for (i = 0; i < numentities_n; ++i) {
              objective += fabs(entityinfo->points_XT[numattributes_m*i+j] - entityinfo->points_XT[numattributes_m*i+l]*v[j]);
            }
          }
          else {
            objective += probleminfo->lambdas_out[0];
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
      /* project data into orthogonal space */
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
    }
  }

  if (probleminfo->lambdas_out[0] < 0.0) { /* enumerate all lambdas and find global opt */
    k=0;
    if (VERBOSITY > 0) {
      fprintf(stdout, "k %d\n", k);
      fflush(stdout);
    }

    lambdas[0] = 0.0;
    num_lambdas = 1;
    for (l = 0; l < numattributes_m; ++l) { /* for each fixed coordinate l */
      if (VERBOSITY > 0) {
        fprintf (stdout, "l %d ", l);
        fflush (stdout);
      }
      /* find v */
      for (j = 0; j < numattributes_m; ++j) {
        lambda_max = 0.0; /* keep track of largest upper bound.  Set v_j = 0 beyond this value */
        num_lambdas_lj[l][j] = 0;
        if (j != l) {
          sumWeights = 0.0;
          index = 0;
          for (i = 0; i < numentities_n; ++i) {
            if (entityinfo->points_XT[numattributes_m*i+l] != 0.0) {
              ratios[index] = entityinfo->points_XT[numattributes_m*i+j]/entityinfo->points_XT[numattributes_m*i+l]; /* store ratios */
              tosort[index]=&(ratios[index]); /* sort the pointers to the ratios */
              weights[index] = fabs(entityinfo->points_XT[numattributes_m*i+l]);
              sumWeights += fabs(entityinfo->points_XT[numattributes_m*i+l]);
              index += 1;
            }
          }
          /* sort */
          qsort(tosort,index,sizeof(double *),sparsEl_pointer_cmp); 

          sum_below=0.0; /* keep track of sum below and sum above; updating involves one addition/subtraction each */
          sum_above=sumWeights;
          for (i = 0; i < index; ++i) {
            v_lj[l][l][num_lambdas_lj[l][j]] = 1.0;
            origIndex = tosort[i] - ratios;
            sum_above -= weights[origIndex];
            /* (lambdaL, lambdaU) is the interval where the current ratio is the best for this j/l combination (l is jhat) */
            lambdaU = (sum_above - sum_below) + weights[origIndex];
            if (*tosort[i] < 0.0) {
              lambdaU = (sum_below - sum_above) + weights[origIndex];
            }
            lambdaL = lambdaU - 2*weights[origIndex];
            /*if ((lambdaU > 0.1) && (lambdaU < 0.2)) {
              fprintf(stdout, "jhat %d (lambdaL, lambdaU) (%f,%f)\n", l, lambdaU - 2*weights[origIndex], lambdaU);
            }*/
            if (lambdaU > 0.0) {  /* if lambdaU <= 0.0, then current point is never the best ratio */
              lambdaL = lambdaU - 2*weights[origIndex];
              lambdas[num_lambdas] = 0.0;  /* if the interval (lambdaL, lambdaU) straddles 0, then we will record the value for lambda = 0 */
              if (lambdaU > lambda_max) {
        	lambda_max = lambdaU;
              }
              v_lj[l][j][num_lambdas_lj[l][j]] = *tosort[i];
              lambdas[num_lambdas] = 0.0; 
              lambdas_lj[l][j][num_lambdas_lj[l][j]] = 0.0;
              lambdas_lj_sort[l][j][num_lambdas_lj[l][j]] = &(lambdas_lj[l][j][num_lambdas_lj[l][j]]);
              if (lambdaL > 0.0) {
                lambdas[num_lambdas] = lambdaL; 
        	lambdas_lj[l][j][num_lambdas_lj[l][j]] = lambdaL;
                lambdas_lj_sort[l][j][num_lambdas_lj[l][j]] = &(lambdas_lj[l][j][num_lambdas_lj[l][j]]);
              }

              if (num_lambdas > 1) {
                if (lambdas[num_lambdas] != lambdas[num_lambdas-1]) { /* no need to store duplicate lambdas */
                  num_lambdas += 1;
                }
              }
              else {
                num_lambdas += 1;
              }

              if (num_lambdas_lj[l][j] > 1) {
                if (lambdas_lj[l][j][num_lambdas_lj[l][j]] != lambdas_lj[l][j][num_lambdas_lj[l][j]-1]) { /* no need to store duplicate lambdas */
                  num_lambdas_lj[l][j] += 1;
                }
              } 
              else {
                num_lambdas_lj[l][j] += 1;
              }
              if (num_lambdas >= max_memory_lambdas) {
                REprintf("Needs more memory 1");
                /*lambdas = (double *) realloc(lambdas, (max_memory_lambdas+numentities_n)*sizeof(double));
                max_memory_lambdas += numentities_n;*/
              }
              if (num_lambdas_lj[l][j] >= max_memory_lj[l][j]) { /* add room for n more solutions */
                REprintf("Needs more memory 2");
                /*lambdas_lj[l][j] = (double *) realloc(lambdas_lj[l][j], (max_memory_lj[l][j]+numentities_n)*sizeof(double));
                lambdas_lj_sort[l][j] = (double **) realloc(lambdas_lj_sort[l][j], (max_memory_lj[l][j]+numentities_n)*sizeof(double *));
                v_lj[l][j] = (double *) realloc(v_lj[l][j], (max_memory_lj[l][j]+numentities_n)*sizeof(double));
                max_memory_lj[l][j] += numentities_n;*/
              }

              /*fprintf(stdout, "v %f\n", *tosort[i]);
              fprintf(stdout, "numlambdas %d\n", num_lambdas_lj[l][j]);*/
            }
            sum_below += weights[origIndex];
        }
          /* keep track of lambda_max solution, even if lambda_max = 0 */
          lambdas_lj[l][j][num_lambdas_lj[l][j]] = lambda_max; /* beyond highest upper bound, set v_j = 0 */
          lambdas_lj_sort[l][j][num_lambdas_lj[l][j]] = &(lambdas_lj[l][j][num_lambdas_lj[l][j]]);
          v_lj[l][j][num_lambdas_lj[l][j]] = 0.0;
          num_lambdas_lj[l][j] += 1;
          lambdas[num_lambdas] = lambda_max;
          num_lambdas += 1;
          if (num_lambdas >= max_memory_lambdas) {
            REprintf("Needs more memory 3");
            /*lambdas = (double *) realloc(lambdas, (max_memory_lambdas+numentities_n)*sizeof(double));
            max_memory_lambdas += numentities_n;*/
            /*fprintf(stderr, "%d\n", max_memory_lambdas);*/
          }
          if (num_lambdas_lj[l][j] >= max_memory_lj[l][j]) { /* add room for n more solutions */
            REprintf("Needs more memory 4");
            /*lambdas_lj[l][j] = (double *) realloc(lambdas_lj[l][j], (max_memory_lj[l][j]+numentities_n)*sizeof(double));
            lambdas_lj_sort[l][j] = (double **) realloc(lambdas_lj[l][j], (max_memory_lj[l][j]+numentities_n)*sizeof(double *));
            v_lj[l][j] = (double *) realloc(v_lj[l][j], (max_memory_lj[l][j]+numentities_n)*sizeof(double));
            max_memory_lj[l][j] += numentities_n;*/
          }

          if (num_lambdas_lj[l][j] > 1) {
            qsort(lambdas_lj_sort[l][j],num_lambdas_lj[l][j],sizeof(double *),sparsEl_pointer_cmp); /* sort lambdas for each l, j combo */
          }
        }
      }
    }
    
    qsort(lambdas,num_lambdas,sizeof(double *),sparsEl_cmp); 

    num_distinct_lambdas = 0;

    for (l = 0; l < numattributes_m; ++l) {
      for (j = 0; j < numattributes_m; ++j) {
        curr_lambda_lj[l][j] = 0;
      }
    }

    lambda_curr = 0;
    lambda_next = 1;
    while (lambda_next < num_lambdas) {
      /*fprintf(stderr, "lambda curr %d next %d %d\n", lambda_curr, lambda_next, num_lambdas);
      fprintf(stdout, "lambda curr %f \n", lambdas[lambda_curr]);
      fprintf(stdout, "lambda next %f \n", lambdas[lambda_next]);
      fflush(stderr);*/
      while (fabs(lambdas[lambda_curr] - lambdas[lambda_next]) < 0.00001) {
        lambda_next += 1;
        if (lambda_next == num_lambdas) {
          break;
        }
      }
      /*fprintf(stdout, "lambda curr %f next %f\n", lambdas[lambda_curr], lambdas[lambda_next]);*/
      for (l = 0; l < numattributes_m; ++l) {
        zs[l] = 0.0; /*obj function  for lambda*/
        vs[l] = 0.0; /* L1 magnitude of v */
        for (j = 0; j < numattributes_m; ++j) {
          v_index = 0;
          if (num_lambdas_lj[l][j] > 1) {
            v_index = lambdas_lj_sort[l][j][curr_lambda_lj[l][j]] - lambdas_lj[l][j];
          }
          /*fprintf(stdout, "num lambdas %d curr_lambda %d %f \n", num_lambdas_lj[l][j], curr_lambda_lj[l][j], lambdas[lambda_curr]); */
          if (curr_lambda_lj[l][j] < num_lambdas_lj[l][j]-1) {
            /*fprintf(stdout, "num lambdas %d curr_lambda %d %f %f \n", num_lambdas_lj[l][j], curr_lambda_lj[l][j], *lambdas_lj_sort[l][j][curr_lambda_lj[l][j]+1], lambdas[lambda_curr]); */
            if (fabs(*lambdas_lj_sort[l][j][curr_lambda_lj[l][j]+1] - lambdas[lambda_curr]) < 0.00001) {
              curr_lambda_lj[l][j] += 1;
              /*fprintf(stdout, "lambda %f\n", *lambdas_lj_sort[l][j][curr_lambda_lj[l][j]]);*/
              v_index = lambdas_lj_sort[l][j][curr_lambda_lj[l][j]] - lambdas_lj[l][j];
              /*if ((l == 9) && (j == 0)) {
                fprintf(stdout, "curr_lambda %f lambda_lj %d %f %ld\n",lambdas[lambda_curr], curr_lambda_lj[l][j],  *lambdas_lj_sort[l][j][curr_lambda_lj[l][j]], lambdas_lj_sort[l][j][curr_lambda_lj[l][j]] - lambdas_lj[l][j]);
              }*/
            }
          }
          for (i=0; i<numentities_n; ++i) {
            /*if ((l == 9) && (j == 0)) {
              fprintf(stdout, "z %d %ld\n", curr_lambda_lj[l][j], lambdas_lj_sort[l][j][curr_lambda_lj[l][j]]- lambdas_lj[l][j]);
            }*/
            zs[l] += fabs(entityinfo->points_XT[numattributes_m*i+j] - entityinfo->points_XT[numattributes_m*i+l]*v_lj[l][j][v_index]); 
            /*fprintf(stdout, "l %d j %d i %d obj %f v_lj %f \n", l, j, i, fabs(entityinfo->points_XT[numattributes_m*i+j] - entityinfo->points_XT[numattributes_m*i+l]*v_lj[l][j][curr_lambda_lj[l][j]]), v_lj[l][j][curr_lambda_lj[l][j]]);*/
            
          }
          /*fprintf(stdout, "v l %d j %d: %f vindex %d\n", l,j, v_lj[l][j][v_index], v_index);*/
          zs[l] += fabs(lambdas[lambda_curr]*v_lj[l][j][v_index]);
          vs[l] += fabs(v_lj[l][j][v_index]);
        }
      }
      for (l = 0; l < numattributes_m; ++l) {
        maxL = 0.0;
        minU = 10000000.0;
        check_lambda = 1;
        for (j = 0; j < numattributes_m; ++j) {
          if (l != j) {
            if (fabs(vs[l] - vs[j]) < 0.00001) {
              if (zs[l] > zs[j]) {
                check_lambda = 0;
                j = numattributes_m + 1;
              }
            }
            else if ((vs[l] < vs[j]) && (maxL < (zs[l]-zs[j])/(vs[j]-vs[l]))) {
              maxL = (zs[l]-zs[j])/(vs[j]-vs[l]);
            }
            else if ((vs[l] > vs[j]) && (minU > (zs[j]-zs[l])/(vs[l]-vs[j]))) {
              minU = (zs[j]-zs[l])/(vs[l]-vs[j]);
            }
          }
        }
        /*fprintf(stderr, "lambdas[t] %f maxL %f minU %f lambdas[t+1] %f \n", lambdas[lambda_curr], maxL, minU, lambdas[lambda_next]);*/
        if ((0.0 < maxL) && (maxL < minU) && (lambdas[lambda_curr] + maxL <= lambdas[lambda_next]) && (check_lambda == 1)) {
          normv = 0.0;
          for (j = 0; j < numattributes_m; ++j) { /* calculate the squared norm of the new line */
            if (num_lambdas_lj[l][j] > 1) {
              normv += v_lj[l][j][lambdas_lj_sort[l][j][curr_lambda_lj[l][j]]- lambdas_lj[l][j]]*v_lj[l][j][lambdas_lj_sort[l][j][curr_lambda_lj[l][j]]- lambdas_lj[l][j]];
            }
            else {
              normv += v_lj[l][j][0]*v_lj[l][j][0];
            }
          }

          for (j = 0; j < numattributes_m; ++j) { /* store the normalized line */
            /*fprintf(stdout, "v_lj %f curr_lambda %f", v_lj[l][j][curr_lambda_lj[l][j]], lambdas_lj[l][j][curr_lambda_lj[l][j]]); */
            if (num_lambdas_lj[l][j] > 1) {
              probleminfo->PCs[numattributes_m*q*num_distinct_lambdas+numattributes_m*k+j] = v_lj[l][j][lambdas_lj_sort[l][j][curr_lambda_lj[l][j]]- lambdas_lj[l][j]]/sqrt(normv); 
            }
            else {
              probleminfo->PCs[numattributes_m*q*num_distinct_lambdas+numattributes_m*k+j] = v_lj[l][j][0]/sqrt(normv); 
            }
          }
          /*fprintf(stdout, " %d \n", new_lambda);*/
          /*fprintf(stdout, "lambda %f\n", lambdas[lambda_curr] + maxL);*/
          probleminfo->lambdas_out[num_distinct_lambdas] = lambdas[lambda_curr] + maxL;
          num_distinct_lambdas += 1;
          if (num_distinct_lambdas >= max_memory) { /* add space for n more solutions */
            REprintf("Needs more memory 5");
            /*probleminfo->PCs = (double *) realloc(probleminfo->PCs, (numattributes_m*q*max_memory+numattributes_m*q*numentities_n)*sizeof(double));
            probleminfo->lambdas_out = (double *) realloc(probleminfo->lambdas_out, (max_memory+numentities_n)*sizeof(double));
            max_memory += numentities_n;*/
          }
        } 
        else if ((maxL <= 0.0) && (minU >= 0.0) && (check_lambda == 1)) {
          normv = 0.0;
          for (j = 0; j < numattributes_m; ++j) { /* calculate the squared norm of the new line */
            if (num_lambdas_lj[l][j] > 1) {
              normv += v_lj[l][j][lambdas_lj_sort[l][j][curr_lambda_lj[l][j]] - lambdas_lj[l][j]]*v_lj[l][j][lambdas_lj_sort[l][j][curr_lambda_lj[l][j]] - lambdas_lj[l][j]];
            }
            else {
              normv += v_lj[l][j][0]*v_lj[l][j][0];
            }
          }

          for (j = 0; j < numattributes_m; ++j) { /* store the normalized line */
            if (num_lambdas_lj[l][j] > 1) {
              /*fprintf(stdout, "v_lj %f curr_lambda %f ", v_lj[l][j][curr_lambda_lj[l][j]], *lambdas_lj_sort[l][j][curr_lambda_lj[l][j]]); */
              probleminfo->PCs[numattributes_m*q*num_distinct_lambdas+numattributes_m*k+j] = v_lj[l][j][lambdas_lj_sort[l][j][curr_lambda_lj[l][j]]-lambdas_lj[l][j]]/sqrt(normv); 
            }
            else {
              probleminfo->PCs[numattributes_m*q*num_distinct_lambdas+numattributes_m*k+j] = v_lj[l][j][0]/sqrt(normv); 
            }
          }
          /*fprintf(stdout, " %d \n", new_lambda);*/
          probleminfo->lambdas_out[num_distinct_lambdas] = lambdas[lambda_curr];
          num_distinct_lambdas += 1;
          if (num_distinct_lambdas >= max_memory) { /* add space for n more solutions */
            REprintf("Needs more memory 6");
            /*probleminfo->PCs = (double *) realloc(probleminfo->PCs, (numattributes_m*q*max_memory+numattributes_m*q*numentities_n)*sizeof(double));
            probleminfo->lambdas_out = (double *) realloc(probleminfo->lambdas_out, (max_memory+numentities_n)*sizeof(double));
            max_memory += numentities_n;*/
          }
        }
      }
      lambda_curr = lambda_next;  
    }

    probleminfo->num_distinct_lambdas = num_distinct_lambdas;


    /*if ((VERBOSITY) >= 1) {

      for (t = 0; t < num_distinct_lambdas; ++t) {
        fprintf (rotationfile, "%f ", lambdas_out[t]);
        for (j = 0; j < numattributes_m; ++j) {
          fprintf (rotationfile, "%f ", probleminfo->PCs[numattributes_m*q*t+numattributes_m*k+j]);
        }
        fprintf (rotationfile, "\n");
      }
    }*/

    /*for (j = 0; j < numattributes_m; ++j) {
      for (l2 = 0; l2 < k; ++l2) {
        for (l1 = 0; l1 < numattributes_m; ++l1) {
            probleminfo->PCs[numattributes_m*k+j] -= probleminfo->PCs[numattributes_m*l2+j]*probleminfo->PCs[numattributes_m*l2+l1]*v[l1];
        }
      }
    }*/
    /* renormalize */
    /*normv = 0.0;
    for (j = 0; j < numattributes_m; ++j) {
      normv += probleminfo->PCs[numattributes_m*k+j]*probleminfo->PCs[numattributes_m*k+j];
    }
    for (j = 0; j < numattributes_m; ++j) {
      probleminfo->PCs[numattributes_m*k+j] = probleminfo->PCs[numattributes_m*k+j]/sqrt(normv);
    }*/

    /* project data into orthogonal space */
    /*for (i = 0; i < numentities_n ; ++i) {
      innerprod = 0.0;
      for (l = 0; l < numattributes_m; ++l) {
        innerprod += probleminfo->PCs[numattributes_m*k+l] * entityinfo->points_XT[numattributes_m*i+l];
      }
      for (j = 0; j < numattributes_m; ++j) {
        entityinfo->points_XT[numattributes_m*i+j] = entityinfo->points_XT[numattributes_m*i+j]-probleminfo->PCs[numattributes_m*k+j]*innerprod;
      }
    }


  if ((VERBOSITY) >= 1) {
    for (j = 0; j < numattributes_m; ++j) {
      for (k = 0; k < q; ++k) {
        fprintf (rotationfile, "%f ", probleminfo->PCs[numattributes_m*k+j]);
      }
      fprintf (rotationfile, "\n");
    }
  }*/
  }
  free(lambdas);
  /*}*/

  return 0;
} /*end solveproblem */

int sparsEl_pointer_cmp(const void *x, const void *y) {
  const double **xx = (const double **)x;
  const double **yy = (const double **)y;
  if (**xx < **yy) return -1;
  if (**xx > **yy) return 1;
  return 0;

}
int sparsEl_cmp(const void *x, const void *y) {
  double xx = *(double*)x, yy=*(double*)y;
  if (xx < yy) return -1;
  if (xx > yy) return 1;
  return 0;
}

