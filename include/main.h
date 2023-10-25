#ifndef MAIN_H
#define MAIN_H

#include <argp.h>
#include <float.h>
#include <math.h>
#include <signal.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "ts_lib.h"

#define UNDET_SYS -1
#define INFEASIBLE_SYS 0
#define INFEASIBLE_SYS_MSG "Schedulable system!"

#define FEASIBLE_SYS 1

#define NO_SOL 0

/* Structure reserved for the result of the CP-KERN algorithm*/

typedef struct {
  int feasible;   /* 1 if the system is feasible, 0 otherwise */
  double opt_sol; /* if feasible, stores the optimal solution */
} cp_kern_res_t;

/* Kernel parameters
typedef struct {
  int n;
  int *alpha;
  int beta;
  int *a, *b;
} kern_param_t;*/

#endif /* MAIN_H */

/*
 * Set the initial values to zero WITHOUT allocating
 */
void cp_kern_res_set_zero(cp_kern_res_t *res);

/* Initialize task set with values from input file */
void ts_initialize();

/* Initialize constants in kernel: alpha, beta, A, B*/
void cp_kern_init_params();

/* Allocate memory for arrays used in the kernel*/
void cp_kern_alloc(int n);

/* Free memory allocated for arrays used in the kernel*/
void cp_kern_free();

/* Computes CP-KERN procedure and stores results in global structure
 * cp_kern_res_t*/
void edf_branching();

/* Returns upper bound L, calculated as follows:
  - if the total utilization is greater than 1, then the hyperperiod is the
  upper bound
  - otherwise, the upper bound is L_b as described in the paper*/
double calculate_upper_bound();

/* Print arrays alpha, x, y and values p, q, r*/
void print_kernel_params(int n, double p, double q, double r);

/* Merge the two sorted parts (before i and after) of pi mantaining a non
 * increasing order of y[pi[i]]*/
void merge_pi(int n, int i, int *pi1, int *pi2, int *pi);
