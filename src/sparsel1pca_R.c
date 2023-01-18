#include "type.h"

static int allocateMemoryL1Line (ENTITYINFOptr entityinfo, PROBLEMINFOptr probleminfo); 
int solveSparsEl(ENTITYINFOptr entityinfo, PROBLEMINFOptr probleminfo);

static void
   free_and_null (char **ptr);

void sparsel1pca (double *points_XT, int *dataDim, int *q, double *PCs, double *objectives, double *lambdas_out) 
{
 
  ENTITYINFO entityinfo;
  SOLVERINFO  solverinfo;
  PROBLEMINFO  probleminfo;
  
  probleminfo.status = 0;
  int status        = probleminfo.status;

  solverinfo.model = NULL;

  probleminfo.obj       = NULL;
  probleminfo.lb        = NULL;
  probleminfo.ub        = NULL;
  probleminfo.rhs       = NULL;
  probleminfo.matbeg    = NULL;
  probleminfo.matind    = NULL;
  probleminfo.matval    = NULL;
  probleminfo.colname   = NULL;
  probleminfo.PCs       = PCs;
  /*probleminfo.getScores = *getScores;
  probleminfo.scores    = scores;*/
  probleminfo.objectives = objectives;

  probleminfo.lambdas        = NULL;
  probleminfo.num_lambdas_lj = NULL;
  probleminfo.v_lj           = NULL;
  probleminfo.lambdas_lj     = NULL;
  probleminfo.curr_lambda_lj = NULL;
  probleminfo.lambdas_out    = lambdas_out;
  probleminfo.max_memory_lj  = NULL;
  probleminfo.vs             = NULL;
  probleminfo.zs             = NULL;
  probleminfo.max_memory     = (MAXLAMBDAS); /* number of unique solutions that we can store */

  
 
  entityinfo.numentities_n   = dataDim[1];
  entityinfo.numattributes_m = dataDim[0];
  
  entityinfo.points_XT = points_XT; /* transpose of data matrix */

  probleminfo.q = *q; /* desired number of PCs */

  status = allocateMemoryL1Line(&entityinfo, &probleminfo);
  if (status) {
    REprintf ("Unable to allocate memory\n");
    goto TERMINATE;
  }

  status = solveSparsEl( &entityinfo, &probleminfo); /* in sparsel1pca.c.c*/
  if (status) {
    REprintf ("Unable to solve.  Terminating...; or done\n");
    goto TERMINATE;
  }
  
  
TERMINATE:

  free_and_null ((char **) &probleminfo.ratios);     
  free_and_null ((char **) &probleminfo.weights);     
  free_and_null ((char **) &probleminfo.v);     
  free_and_null ((char **) &probleminfo.tosort);

  if (lambdas_out[0] < 0.0) {
    for (probleminfo.l = 0; probleminfo.l < entityinfo.numattributes_m; ++probleminfo.l){
      free_and_null ((char **) &probleminfo.num_lambdas_lj[probleminfo.l]);
      free_and_null ((char **) &probleminfo.curr_lambda_lj[probleminfo.l]);
      free_and_null ((char **) &probleminfo.max_memory_lj[probleminfo.l]);
      for (probleminfo.j = 0; probleminfo.j < entityinfo.numattributes_m; ++probleminfo.j){
        free_and_null ((char **) &probleminfo.lambdas_lj[probleminfo.l][probleminfo.j]);
        free_and_null ((char **) &probleminfo.v_lj[probleminfo.l][probleminfo.j]);
      }
    }
    for (probleminfo.l = 0; probleminfo.l < entityinfo.numattributes_m; ++probleminfo.l){
      free_and_null ((char **) &probleminfo.lambdas_lj[probleminfo.l]);
      free_and_null ((char **) &probleminfo.v_lj[probleminfo.l]);

    }
    free_and_null ((char **) &probleminfo.v_lj);
    free_and_null ((char **) &probleminfo.lambdas_lj);
    free_and_null ((char **) &probleminfo.max_memory_lj);
    free_and_null ((char **) &probleminfo.curr_lambda_lj);
    free_and_null ((char **) &probleminfo.num_lambdas_lj);
    free_and_null ((char **) &probleminfo.v_lj);
    free_and_null ((char **) &probleminfo.lambdas_lj);
    free_and_null ((char **) &probleminfo.curr_lambda_lj);
    free_and_null ((char **) &probleminfo.max_memory_lj);
    free_and_null ((char **) &probleminfo.vs);
    free_and_null ((char **) &probleminfo.zs);
  }

} /* end l1line */

static void
free_and_null (char **ptr) {
  if ( *ptr != NULL ) {
     free (*ptr);
     *ptr = NULL;
  }
} /* END free_and_null */  

static int allocateMemoryL1Line (ENTITYINFOptr entityinfo, PROBLEMINFOptr probleminfo) {
  int numentities_n   = entityinfo->numentities_n;
  int numattributes_m = entityinfo->numattributes_m;
  int j               = probleminfo->j;
  int l               = probleminfo->l;

  probleminfo->ratios     = (double *) malloc (numentities_n*sizeof (double));
  probleminfo->tosort     = (double **) malloc (numentities_n*sizeof (double *));
  probleminfo->weights    = (double *) malloc (numentities_n*sizeof (double));
  probleminfo->v = (double *) malloc(numattributes_m*sizeof(double));

  if (probleminfo->lambdas_out[0] < 0.0) { 
    probleminfo->max_memory_lj = (int **) malloc(numattributes_m*sizeof(int *)); 
    probleminfo->max_memory_lambdas = probleminfo->max_memory; 

    /*probleminfo->PCs = (double *) Realloc(probleminfo->PCs, numattributes_m*q*probleminfo->max_memory, double);
    probleminfo->lambdas_out = (double *) Realloc(probleminfo->lambdas_out, probleminfo->max_memory, double);*/ /* to reallocate memory, use .Call instead of .C */

    probleminfo->lambdas = (double *) malloc(probleminfo->max_memory*sizeof(double));
    probleminfo->num_lambdas_lj = (int **) malloc(numattributes_m*sizeof(int *));
    probleminfo->lambdas_lj = (double ***) malloc(numattributes_m*sizeof(double **));
    probleminfo->lambdas_lj_sort = (double ****) malloc(numattributes_m*sizeof(double ***));
    probleminfo->curr_lambda_lj = (int **) malloc(numattributes_m*sizeof(int *));
    probleminfo->v_lj = (double ***) malloc(numattributes_m*sizeof(double));
    probleminfo->zs = (double *) malloc(numattributes_m*sizeof(double));
    probleminfo->vs = (double *) malloc(numattributes_m*sizeof(double));
    for (l = 0; l < numattributes_m; ++l) {
      probleminfo->num_lambdas_lj[l] = (int *) malloc(numattributes_m*sizeof(int));
      probleminfo->lambdas_lj[l] = (double **) malloc(numattributes_m*sizeof(double *));
      probleminfo->lambdas_lj_sort[l] = (double ***) malloc(numattributes_m*sizeof(double **));
      probleminfo->curr_lambda_lj[l] = (int *) malloc(numattributes_m*sizeof(int));
      probleminfo->v_lj[l] = (double **) malloc(numattributes_m*sizeof(double *));
      probleminfo->max_memory_lj[l] = (int *) malloc(numattributes_m*sizeof(int));
      for (j = 0; j < numattributes_m; ++j) { /* allocate memory for n solutions initially */
        probleminfo->max_memory_lj[l][j] = probleminfo->max_memory;
        probleminfo->v_lj[l][j] = (double *) malloc(probleminfo->max_memory*sizeof(double));
        probleminfo->lambdas_lj[l][j] = (double *) malloc(probleminfo->max_memory*sizeof(double));
        probleminfo->lambdas_lj_sort[l][j] = (double **) malloc(probleminfo->max_memory*sizeof(double *));
      }
    }
  }

  return 0;
} /* end allocateMemoryL1Line */

