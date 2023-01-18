#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>
#include <assert.h>
#include <float.h>
#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <Clp_C_Interface.h>
#include <ClpConfig.h>
#define PATHLENGTH 500 /* length of path name describing location of data file*/
#define OBJ_INIT DBL_MAX /* initial incumbent value */
#define VERBOSITY 0  /* pcal1: 3 - print L2 estimate of w at each iteration.  l1pcastar: 1- print projdim, rotation matrix, average error, l 3- projected points in terms of original coords, 4-screen output on, all objective values, best betas, projected points at each iteration, orthogonal direction, columns of rotation matrix at each iteration, initial SVD info, 7 - write points to screen*/
#define NBMAX 64 /* maximum possible block size, see dgeqrf.f in LAPACK */
#define EPSILON 0.0000001 /* check if objective is 0, check if wT == wTold in pcal1, check if norm of x is 0 in pcal1 */
#define MAXLAMBDAS 1000 /* maximum number of lambdas to store */

/* input data */
struct entityinfo {
  int numentities_n;/*# of rows/points in your data*/
  int numattributes_m;/*# of columns*/
  double *points_XT;/*actual data points in column major format (i.e., transposed) - in a single array*/
  double *PCs; /* loadings for solving l1projection */
};
typedef struct entityinfo ENTITYINFO, *ENTITYINFOptr;

/* main CPLEX variables (others are in cutcallbackinfo) */
struct solverinfo {
/*  *CPXENVptr env;
  CPXLPptr  lp;*/
  struct Clp_Simplex  *model;
  int       status;
  char      *errmsg;
  struct Clp_Simplex *modelU;
  struct Clp_Simplex *modelV;
};
typedef struct solverinfo SOLVERINFO, *SOLVERINFOptr;


/* problem info */
struct probleminfo {
  /* variables shared by two or more methods */
  int rcnt;/*# of constraints in LP*/
  int nzcnt;/*number of coefficients to change*/
  int *matbeg;/*starting index of column*/
  int *matind;/*list of row-indices for each value in matval*/
  double *matval;/*co-efficients for the linear constraints*/
  int    q;/*projected dimension*/
  int    solstat;/*solution status -0-optimal*/      
  double *b;/*matrix of PC's*/
  int    numcols;/*# of columns for constraint matrix*/
  int    i;
  int    j;
  int    k;
  int    status;/*checking for the success of function*/
  double *obj; /*coefficients of objective function*/
  double *lb;/*lower bound for variables*/
  double *ub;/*upper bound for variables*/
  double *work;/*used in Fortran routine 'dgesvd'*/
  int    lwork;/*used in Fortran routine 'dgesvd'*/
  double *S;/*eigenvalues in SVD*/

  /* variables for L1-PCA-star */
  double *rhs;/*vector of RHS-0 for all of them*/
  int    *betaind;/*index of betas*/    
  int    *eplusind;/*indices for e+*/   
  int    *eminusind;  /*indices for e-*/
  double objective;/*stores the objective CLP returns*/
  int    projdim;/*index for the loop-dimension we're projecting into (k-1)*/     
  int    *bestdir;/*stores the best direction for every dimension*/      
  double *mybeta;/*betas-coefficients for datapoints in constraint*/
  const double *currBeta; /* pointer to solution for beta */ 
  double minobjective;/*stores the best objective for k linear regression*/
  int    l;/*variable for loop on regression -j (article)*/
  char   **colname;/*variable names*/
  double *xpluslambda_Z;/*Z-projected points on the subspace in terms of original co-ordinates*/
  double *xpluslambda_Z2;/*copy of Z*/
  double *Vj;/*V matrix*/
  double *VT;/*U matrix*/
  double *a;/*normalized beta*/
  double *preVj;/*product of V matrix*/
  double *temppreVj;/*used for obtaining the product of V matrices*/
  double *VBeta;/*product of VjT with Identity matrix (which has a row of betas for best direction dimension)*/
  double *preVBeta;/*transformation for finding the new score of a point*/
  double *temppreVBeta;
  double *tempPC;/**/
  double *projPoints;/*projected points in the subspace in terms of original co-ordinates*/
  int    numprojdim; /* number of projection directions with "L1 variation" greater than 1 */
  int    numfactors;/*min (#of attributes, # of points)*/
  double *scores; /*stores scores- projected points in terms of new coordinates*/ 
  int getScores; /*0 for no scores; default is 0*/
  int getProjPoints; /* 0 for no projected point calculations; default is 0*/

  /*variables for PCA-L1 */

  int initMethod; /* initialization method 0 - use L2-PCA first component, 1 - use direction of entity with largest norm, 2 - random */
  double *polarity;/*polarities for each column in data matrix */
  int dotConv;/*stores if w.x(i)=0*/
  double *wT;/*Tth w*/
  double *wTOld;/*(T-1)th w*/
  double *Xj; /*points for the Xj matrix -See Kwak*/
  double innerprod; /*w^Tx */
  int    convergent; /* is pcaL1 convergent */
  double wTSum; /* for finding norm-squared of w */
  double Normalizer; /* sqrt of wTSum */
  double dotProd; /* wT . x */
  double normx; /* x . x */
  double normalizer; /* for normalizing new direction to unit length  */
  int    argmax; /* storing the entity with the largest norm */
  double x; /* for finding norm of points */
  double maxSum; /* for finding norm of points */
  double xSum; /* norm squared of point */
  double *points_XT_temp; /* temporary storage of points */
  double *PCs;/*PCs for PCA-L1*/

  /*for L1-PCA  */
  int initialize; /*whether to use random initialization or not*/
  double *initV; /*initial V */ 
  double *V; /* rotation matrix */
  double *Vtemp;/*temp V matrix*/
  double *nv; /*normalized diagonal vector for V*/
  double *U; /* scores matrix */
  double *Utemp; /*temp U matrix*/
  double *nu; /*normalized diagonal vector for U*/
  double *sigma;
  double tolerance; /*user defined epsilon*/
  int iterations; /*user defined no. of iterations for algorithm*/
  int    iter; /* iteration number */
  double maxdifference; /* max difference between elements of U and V between iterations*/
  double *rhsL; /* lower bound for constraints */
  double *rhsU; /* upper bound for constraints */
  int    *xinda; /* index of constraints when solving for U */
  int    *xindb; /* index of constraints when solving for U */
  int    *uind; /* index of variables when solving for U */
  const double *Usol; /* solution to problem when solving for U, including delta's */
  int    **rowind; /* index of constraints when solving for V */
  int    **vind; /* index for variables when solving for V */
  const double *Vsol; /* solution to problem when solving for V, including lambda's */
  double udiff; /* checking for convergence of U */
  double vdiff; /* checking for convergence of V */

  /*for SharpEl  */
  double *ratios; /* ratios x_j/x_j* */
  double **tosort; /* pointers to the ratios to be sorted */
  double *weights; /* weights for points */
  double medWeight; /* weighted median */
  double sumWeights; /* the sum of the weights */
  double *v;  /* optimal line */
  double normv;  /* norm of v */
  int origIndex; /* to keep track of original index of sorted values */ 
  int index; /* for keeping track of nonzero and defined ratios */
  int lstar; /* direction to keep fixed for best-fit line */
  double *objectives; /* stores objective function values for each projection dimension */
  int normalizeV; /* Whether to normalize the components */
  int l1;
  int l2;
  
  /*for L1-PCAHp  */
  double threshold; /*user define threshold*/
  double currObj; /*current objective function*/
  double xx_obj; /*previous objective function*/
  double x_obj; /*previous objective function*/
  double *point; /*most recent projection we ar trying to improve*/
  const double *alpha; /*most recent projection we ar trying to improve*/
  int numrows; /* also used in l1projection */

  /* for l1projection */
  int *aind; /* index of alpha variables */
  const double *projSolution; /* LP solution */ 
  double *alphas;/* scores */

  /*for PCA-Lp  */
  int solMethod; /*whether to use Gradient Ascendent or Lagrangian */
  double lratio; /*learning ratio*/
  double p; /*p-norm that will be applied */
  double epsilon; /*user define tolerance for the convergence check*/
  double *Dw; /*gradient */
  double pNormValue; /* wT . x */

  /* for sparsEl */
  double *lambdas; /* breakpoints for regularization parameter */
  int num_lambdas; /* number of breakpoints */
  double lambda_max; /* largest breakpoint; beyond this value, v_j = 0 */
  int **num_lambdas_lj; /* number of breakpoints for fixed coord l and coord j */
  double sum_below; /* keep track of sum of weights for weighted median */
  double sum_above;/* keep track of sum of weights for weighted median */
  double lambdaU; /* upper bound of interval where a ratio is optimal for v_j */
  double lambdaL;/* lower bound of interval where a ratio is optimal for v_j */
  double ***v_lj; /* solutions for fixed coord l and coord j and for each breakpoint */
  double ***lambdas_lj; /* breakpoints for fixed coord l and coord j */
  double ****lambdas_lj_sort; /* ordered breakpoints */
  int **curr_lambda_lj; /* for iterating through lambdas when comparing solutions across fixed coord l */
  double *lambdas_out;  /* final breakpoints for output; given lambda on input */
  int num_distinct_lambdas; /* number of breakpoints */
  int max_memory; /* keep track of memory for breakpoints and solutions for output; upper bound on number of breakpoints */
  int **max_memory_lj;  /* keep track of memory for breakpoints and solutions for each fixed coord l and coord j */
  int max_memory_lambdas; /* keep track of memory for breakpoints and solutions from Algorithm 1 */ 
  int sparse; /* whether to solve sparse method */ 
  double *vs; /* keep track of solutions for each fixed coord l */
  double *zs; /* keep track of objective function values for each coord l */
  double maxL; /* keep track of lower bound on interval for penalty where solution is best */
  double minU;/* keep track of lower bound on interval where solution is best */
  int lambda_curr; /* keep track of index of current penalty parameter */
  int lambda_next; /* keep track of index of next penalty parameter (avoid duplicates) */
  int check_lambda; /* keep track of whether to check the interval for this fixed coord l */
  int v_index; /* keep track of which v to select for a fixed coordinate l and coordinate j for a given penalty */
  int t; /* index for distinct lambdas */

};
typedef struct probleminfo PROBLEMINFO, *PROBLEMINFOptr;
