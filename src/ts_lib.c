#include "ts_lib.h"

#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*
 * Helper macros to compact a bit the code. They are relative to a
 * data structure of type ts_t, with name cur_ts
 */
#define N (cur_ts->num)
#define T(i) (cur_ts->per[(i)])
#define D(i) (cur_ts->dl[(i)])
#define C(i) (cur_ts->wcet[(i)])
#define DD(i) (cur_ts->dl[(i)] - cur_ts->jtr[(i)])
#define J(i) (cur_ts->jtr[(i)])
#define HAS_JTR (cur_ts->has_jtr)
#define H (cur_ts->h_per)
#define EPS (cur_ts->eps)

/*
 * It computes the least common multiple H of all periods (often
 * called hyperperiod). The algorithm guarantees exactness by at least
 * EPS: there is certainly a multiple of any period in the interval
 * [H, H*(1+EPS)). If (cur_ts->h_per_tol == 0) then the hyperperiod is
 * exact.
 */
static void compute_hyperperiod(ts_t *cur_ts) {
  unsigned long *jobs;
  double *multi;
  int id_max = 0, id_min = 0, i;

  /* Number of jobs per task within (current) hyperperiod  */
  jobs = malloc(sizeof(*jobs) * N);
  for (i = 0; i < N; i++)
    jobs[i] = 1;
  /* Array of jobs[i]*T(i) */
  multi = malloc(sizeof(*multi) * N);
  memcpy(multi, cur_ts->per, sizeof(*multi) * N);
  /* searching for min and max */
  for (i = 1; i < N; i++) {
    if (T(i) > T(id_max))
      id_max = i;
    if (T(i) < T(id_min))
      id_min = i;
  }
  while (multi[id_max] / multi[id_min] - 1 > EPS) {
    if (++jobs[id_min] >= ULONG_MAX) {
      fprintf(stderr, "%s:%d: Exceeded precision\n", __FILE__, __LINE__);
      H = 0;
      return;
    }
    multi[id_min] = jobs[id_min] * T(id_min);
    if (multi[id_min] > multi[id_max])
      id_max = id_min;
    for (i = 0; i < N; i++) {
      if (multi[i] < multi[id_min])
        id_min = i;
    }
  }
  H = multi[id_min];
  cur_ts->h_per_tol = multi[id_max] / multi[id_min] - 1;
  free(jobs);
  free(multi);
}

/*
 * Compute the maximum deadline and minimum DD (deadline minus period and
 * jitter)
 */
static void store_max_dl(ts_t *cur_ts) {
  int i;

  /* Computing max deadline */
  cur_ts->max_d = D(0);
  cur_ts->min_dd = D(0) - J(0);
  for (i = 1; i < N; i++) {
    if (D(i) > cur_ts->max_d)
      cur_ts->max_d = D(i);
    if (DD(i) < cur_ts->min_dd)
      cur_ts->min_dd = DD(i);
  }
}

void ts_set_zero(ts_t *cur_ts) {
  cur_ts->num = 0;
  cur_ts->per = NULL;
  cur_ts->dl = NULL;
  cur_ts->jtr = NULL;
  cur_ts->wcet = NULL;
}

void ts_realloc(ts_t *cur_ts) {
  /* allocate the data structure */
  cur_ts->per = (double *)realloc(cur_ts->per, sizeof(double) * N);
  cur_ts->dl = (double *)realloc(cur_ts->dl, sizeof(double) * N);
  cur_ts->jtr = (double *)realloc(cur_ts->jtr, sizeof(double) * N);
  cur_ts->wcet = (double *)realloc(cur_ts->wcet, sizeof(double) * N);
}

void ts_read_alloc(ts_t *cur_ts, char *input_filename) {
  int i, num;
  double cur_T, cur_D, cur_J, cur_C;
  int cur_phasing;

  FILE *file_ptr = fopen(input_filename, "r");
  if (NULL == file_ptr) {
    fprintf(stderr, "%s:%d: Input file can't be opened.\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }
  /* initially set no phasing */
  cur_phasing = 0;

  /* read the number of tasks */
  fscanf(file_ptr, "%d", &num);
  if (N != num) {
    N = num;
    /* Now alloc proper amount */
    ts_realloc(cur_ts);
  } else {
    /*Nothing, already allocated */
  }

  fscanf(file_ptr, "%lf", &(EPS));

  for (i = 0; i < N; i++) {
    /* read the task set parameters */
    fscanf(file_ptr, "%lf\t%lf\t%lf\t%lf", &cur_T, &cur_D, &cur_J, &cur_C);
    T(i) = cur_T;
    D(i) = cur_D;
    cur_J = cur_J - floor(cur_J / cur_T) * cur_T;
    J(i) = cur_J;
    C(i) = cur_C;
  }

  cur_ts->constrained_dls = 1;

  /* Closing input file once the task set has been read*/
  fclose(file_ptr);

  compute_hyperperiod(cur_ts);
  store_max_dl(cur_ts);

  /*Checking if it is a constrained-deadline system*/
  for (i = 0; i < N; i++) {
    if (D(i) > T(i)) {
      cur_ts->constrained_dls = 0;
    }
  }
}

void ts_print(const ts_t *cur_ts) {
  int i;

  printf("\nTask system with %s release jitter\n",
         HAS_JTR ? "NON-ZERO" : "ZERO");
  printf(" i: (      Ti,       Di,       Ji,       Ci)\n");
  for (i = 0; i < N; i++) {
    printf("%2i: (%8.1f, %8f, %8f, %8f)\n", i, T(i), D(i), J(i), C(i));
  }
  printf("Hyperperiod: %f, computed with tolerance %g (if == 0, then exact)\n",
         H, cur_ts->h_per_tol);
  fflush(stdout);
}

void ts_free(ts_t *cur_ts) {
  free(cur_ts->per);
  free(cur_ts->dl);
  free(cur_ts->jtr);
  free(cur_ts->wcet);
}

void ts_sort_tasks(ts_t *cur_ts) {
  int i, prev_i;
  double tmp, tmp_per, tmp_dl, tmp_jtr, tmp_wcet;
  double *key_arr = malloc(sizeof(double) * cur_ts->num);

  /* Copying the task set to arrays */
  for (i = 0; i < cur_ts->num; i++) {
    key_arr[i] = cur_ts->dl[i] - cur_ts->jtr[i] - cur_ts->per[i];
  }

  /* Insertion sort (preferable for small arrays) */
  for (i = 1; i < cur_ts->num; i++) {
    tmp = key_arr[i];
    prev_i = i - 1;

    tmp_per = cur_ts->per[i];
    tmp_dl = cur_ts->dl[i];
    tmp_jtr = cur_ts->jtr[i];
    tmp_wcet = cur_ts->wcet[i];

    while (prev_i >= 0 && key_arr[prev_i] > tmp) {
      key_arr[prev_i + 1] = key_arr[prev_i];
      cur_ts->per[prev_i + 1] = cur_ts->per[prev_i];
      cur_ts->dl[prev_i + 1] = cur_ts->dl[prev_i];
      cur_ts->jtr[prev_i + 1] = cur_ts->jtr[prev_i];
      cur_ts->wcet[prev_i + 1] = cur_ts->wcet[prev_i];

      --prev_i;
    }

    key_arr[prev_i + 1] = tmp;
    cur_ts->per[prev_i + 1] = tmp_per;
    cur_ts->dl[prev_i + 1] = tmp_dl;
    cur_ts->jtr[prev_i + 1] = tmp_jtr;
    cur_ts->wcet[prev_i + 1] = tmp_wcet;
  }

  free(key_arr);
}
