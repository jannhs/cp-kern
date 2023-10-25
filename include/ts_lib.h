#ifndef _TS_LIB_H
#define _TS_LIB_H
#include <stdio.h>

/* TASK SET DATA STRUCTURE */

/*
 * Task set  related data. It  is preferred to store  homogeneous data
 * (periods,  deadline,...) in  unique vector,  rather than  storing a
 * vector  of  "task  structures",  where each  element  contains  the
 * parameters. This is for efficiency.
 */
typedef struct {
  int num;          /* number of tasks */
  double *per;      /* periods */
  double *dl;       /* rel. deadline (even > period) */
  double *jtr;      /* (J) max jitter experienced by task */
  double *wcet;     /* worst case execution time */
  int has_jtr;      /* 0 if system has zero release jitter*/
  double h_per;     /* hyperperiod = lcm of periods */
  double h_per_tol; /* h_per tolerance. If 0, exact */
  double max_d;     /* max deadline */
  double min_dd; /* D̂_min (where D̂ is obtained subtracting period and jitter to
                    the corresponding task's deadline) */
  int constrained_dls; /* 1 if all deadlines are constrained, 0 otherwise */
  double eps;          /* (t_1-t_0)/t_1 < eps => t_1, t_0 are same */
} ts_t;

/*
 * Set the initial values to zero WITHOUT allocating
 */
void ts_set_zero(ts_t *cur_ts);

/*
 * Allocate arrays for task set. The size is taken from
 * cur_ts->number. If cur_ts->number == 0, then the pointers MUST be
 * NULL
 */
void ts_realloc(ts_t *cur_ts);

/*
 * Initialize the task set using data read from stdin
 */
void ts_read_alloc(ts_t *cur_ts, char *input_filename);

/*
 * Print the data of the task set
 */
void ts_print(const ts_t *cur_ts);

/*
 * Free all the dynamically allocated memory
 */
void ts_free(ts_t *cur_ts);

/*
 * Sort tasks in a non-decreasing order using the key i = D_i - J_i - T_i
 */
void ts_sort_tasks(ts_t *cur_ts);

#endif /* _TS_LIB_H */
