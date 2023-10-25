#include "main.h"

char *input_file = "input.txt";
ts_t cur_ts;       /* Current task set*/
cp_kern_res_t res; /* CP-KERN algorithm's results*/

/* KERNEL PARAMETERS*/
double *alpha;
double beta;
double *A, *B;
double L;
double *x, *y;
int *pi, *pi1, *pi2;

int cp_kern(int k);

int main(int argc, char *argv[]) {
  ts_initialize();
  cp_kern_res_set_zero(&res);
  cp_kern_init_params();

  /*edf_branching();*/
  cp_kern_fp();

  if (res.feasible == FEASIBLE_SYS) {
    printf("Unschedulable system: missed deadline at %f\n", res.opt_sol);
  } else if (res.feasible == UNDET_SYS) {
    printf("Undetermined system!\n");
  } else {
    printf("Schedulable system!\n");
  }

  cp_kern_free();
}

double calculate_upper_bound() {
  int i;
  double L = 0.0;
  double first_candidate, second_candidate;
  double max_diff = -DBL_MAX;
  double sum_numerator = 0.0;
  double temp, total_U = 0.0;

  /* Total amount of utilization*/
  for (i = 0; i < cur_ts.num; i++) {
    printf("U_%d = %f / %f = %f\n", i, cur_ts.wcet[i], cur_ts.per[i],
           cur_ts.wcet[i] / cur_ts.per[i]);
    total_U += cur_ts.wcet[i] / cur_ts.per[i];
  }

  printf("Total utilization: %f\n", total_U);

  printf("Upper bound = max {");

  if (total_U < 1) {

    for (i = 0; i < cur_ts.num; i++) {
      temp = cur_ts.dl[i] - (cur_ts.jtr[i] + cur_ts.per[i]);
      if (max_diff < temp) {
        max_diff = temp;
      }
    }

    first_candidate = max_diff;
    printf("%f,", first_candidate);

    for (i = 0; i < cur_ts.num; i++) {
      sum_numerator += (cur_ts.per[i] - cur_ts.dl[i] + cur_ts.jtr[i]) *
                       (cur_ts.wcet[i] / cur_ts.per[i]);
    }

    /*Implementation as described in algorithm 1*/
    /*second_candidate = sum_numerator / (1 - total_U);*/

    /*Implementation as described in rtsched/sched_test/edf.py */
    second_candidate = floor(sum_numerator / (1 - total_U) - 1);

    printf("%f}", second_candidate);

    if (first_candidate < second_candidate) {
      L = second_candidate;
    } else {
      L = first_candidate;
    }
  } else { /* Upper bound will be the hyper-period of the system*/
    L = cur_ts.h_per;
  }
  printf(" -> %f\n", L);
  return L;
}

void ts_initialize() {
  ts_set_zero(&cur_ts);
  ts_read_alloc(&cur_ts, input_file);
  printf("INITIAL TASK SET:");
  ts_print(&cur_ts);
  printf("\n--- Sorting task set ---\n");

  ts_sort_tasks(&cur_ts);
  printf("SORTED TASK SET:");
  ts_print(&cur_ts);
}

#define DD(i) (cur_ts.dl[i] - cur_ts.jtr[i])
#define DD_min (cur_ts.min_dd)
#define T(i) (cur_ts.per[i])
#define J(i) (cur_ts.jtr[i])

void cp_kern_init_params() {
  int i;
  alpha = malloc(cur_ts.num * sizeof(double));
  A = malloc(cur_ts.num * sizeof(double));
  B = malloc(cur_ts.num * sizeof(double));
}

void edf_branching() {
  int i, h;
  int k;         /* index of the task being considered */
  int found_sol; /* 1 if an optimal solution exists, 0 otherwise*/

  /* INIITIALIZING KERNEL PARAMETERS */

  int p, q;

  L = calculate_upper_bound();

  /* alpha_i = D_i - J_i -  T_i */
  for (i = 0; i < cur_ts.num; i++) {
    alpha[i] = cur_ts.dl[i] - cur_ts.jtr[i] - cur_ts.per[i];
  }

  beta = 1; /* or zero? */

  printf("D̂_min: %f\n", cur_ts.min_dd);

  if (DD_min > L) {
    printf("D̂_min > L\n");
    res.feasible = INFEASIBLE_SYS;
    return;
  }

  if (cur_ts.constrained_dls) {
    printf("Constrained deadlines! Only one interval to check.\n");
    p = q = cur_ts.num - 1;
  } else {
    /*
     * Computing p and q to restrict search space:
     * - p is the first index such that D̂_p - T_p >= D̂_min
     * - q is the last index such that D̂_q - T_q < L (in this instance L is
     *the hyperperiod)
     **/
    for (i = 0; i < cur_ts.num && DD(i + 1) - T(i + 1) < cur_ts.min_dd; i++) {
    }
    p = i;
    for (i = cur_ts.num - 1; i >= 0 && DD(i) - T(i) >= L; i--) {
    }
    q = i;
  }

  printf(
      "\n--- Initializing intervals in EDF-BRANCHING over [p,q]=[%d,%d] ---\n",
      p, q);

  /*Print a e b*/
  for (h = q; h >= p; h--) {
    if (h == cur_ts.num - 1) {
      B[h] = L;

    } else {
      B[h] = DD(h + 1) - T(h + 1);
    }
    if (DD_min >= alpha[h]) {
      A[h] = DD_min;
    } else {
      A[h] = alpha[h];
    }
    printf("[a_%d,b_%d]=[%f,%f]\n", h, h, A[h], B[h]);
  }

  k = q;
  while (k >= p) {
    cp_kern_alloc(k + 1);
    found_sol = cp_kern_edf(k); /*MUST RETURN NEGATION OF OPTIMAL VALUE*/
    if (found_sol) {            /*optimal solution found*/
      return;
    }
    k--;
  }
}

/* Kernel's parameters for FP schedulability test:
 * - n = n
 * - alpha =  J
 * - beta = 0
 * - a = 1
 * - b = D̂_n - T_n
 * */
int cp_kern_fp() {
  double a, b;
  int i = 0, k = 0, h = 0;
  int h1, h2;
  int j, iter = 0;
  int i_start; /*smallest point of the domain of f*/
  int n = cur_ts.num;
  double p, q, r, d;

  cp_kern_alloc(n);

  a = 1;
  b = DD(n - 1) - J(n - 1);

  printf("\n --- Initializing parameters in CP-KERN over "
         "[a,b]=[%f,%f]----\n",
         a, b);

  if (a > b) {
    res.feasible = INFEASIBLE_SYS;
    return INFEASIBLE_SYS;
  }

  /*None tasks*/
  if (n == 1) {
    res.feasible = FEASIBLE_SYS;
    res.opt_sol = a;
    return FEASIBLE_SYS;
  }

  /*Initializing p,q and r*/
  p = r = beta;
  q = 1;

  /* Initializing x, y and pi*/
  for (h = 0; h < n; h++) {
    x[h] = ceil((a + alpha[h]) / cur_ts.per[h]);
    if (x[h] == -0.0) {
      x[h] = 0.0;
    }
    printf("x_%d = ceil((%.2f + %.2f) / %.2f) = %.2f\n", i, a, alpha[h],
           cur_ts.per[h], x[h]);

    y[h] = cur_ts.per[h] * x[h] - alpha[h];
    printf("y_% = %.2f * %.2f - %.2f = %.2f\n", i, cur_ts.per[h], x[h],
           alpha[h], y[h]);
    pi[h] = h; /* pi stores the indexes of y so that for any i <= j y[pi[i]]
                      <= y[pi[j]], right now is the identity permutation*/
    p += (cur_ts.wcet[h] * alpha[h]) / cur_ts.per[h];
    q -= (cur_ts.wcet[h] / cur_ts.per[h]);
    r += cur_ts.wcet[h] * x[h];
  }

  print_kernel_params(n, p, q, r);

  /**
   * Checking infeasibility condition that states that if
   *
   *  (p) = ( beta + sum of (U_i * alpha_i) ) > 0
   *
   *                 and
   *
   *  (q) = ( 1 - sum of U_i ) == 0
   *
   *  then the IP-KERN is infeasible.
   *
   * */

  printf(
      "\n--- Checking infeasibility condition && initializing i_start ---\n");
  if (p > 0 && q == 0) {
    res.feasible = INFEASIBLE_SYS;
    return INFEASIBLE_SYS;
  }

  if (q == 0) {
    i_start = 1;
  } else {
    i_start = 0;
  }
  printf("i_start = %d\n", i_start);

  printf("\n--- Sorting pi in a nonincreasing order of y[pi[i]] ---\n");

  printf("pi={");
  for (h = 0; h < n; h++) {
    printf("%d ", pi[h]);
  }
  printf("} -> {");

  /* Sorting pi in a nonincreasing order using the key y[pi[h]]*/
  for (h = 1; h < n; h++) {
    int tmp = pi[h];
    int prev_h = h - 1;

    while (prev_h >= 0 && y[pi[prev_h]] < y[tmp]) {
      pi[prev_h + 1] = pi[prev_h];
      --prev_h;
    }
    pi[prev_h + 1] = tmp;
  }

  for (h = 0; h < n; h++) {
    printf("%d ", pi[h]);
  }
  printf("}\n");

  printf(" -> y[pi[i]]={");
  for (h = 0; h < n; h++) {
    printf("%f ", y[pi[h]]);
  }
  printf("}\n");

  while (iter < 4) {
    p = res.opt_sol = r;
    q = 1.0;
    i = n - 1;

    /* Solve LP*/
    printf("\n--- Computing optimal value for LP (iteration n.%d) ---",
           iter + 1);

    printf("\n1. Setting p=%f, q=%f, i=%d.\n", p, q, i);
    while (i >= i_start) {
      k = pi[i];
      printf("2. k=%d:\n", k);
      printf("\tTask %d with C=%.3f, T=%.3f, alpha=%.3f\n", k, cur_ts.wcet[k],
             cur_ts.per[k], alpha[k]);
      printf("\t%.3f <= %.3f*%.3f=%.3f ?", p, q, y[k], q * y[k]);
      if (p <= q * y[k]) {
        printf(" YES! Possible solution =%.3f\n", p / q);
        res.opt_sol = p / q;
        break;
      } else {
        printf("\tNO\n");
      }

      p -= cur_ts.wcet[k] * x[k] + cur_ts.wcet[k] * alpha[k] / cur_ts.per[k];
      q -= cur_ts.wcet[k] / cur_ts.per[k];
      i--;

      printf("\tUpdated p=%f, q=%f, i=%d\n", p, q, i);
    }

    printf("3. i=%d, i_start=%d: (i==i_start)? ", i, i_start);
    if (i == i_start - 1) {
      res.opt_sol = p / q;
      printf("YES, possibile solution: %f.\n", res.opt_sol);
    } else {
      printf("NO\n");
    }

    if (res.opt_sol <= a) {
      res.opt_sol = a;
      res.feasible = FEASIBLE_SYS;
      printf("4a. Solution found (res.opt_sol <= a): %f\n", a);
      return FEASIBLE_SYS;
    }

    printf("4b. Checking if value is greater than b: %f > %f ? ", res.opt_sol,
           b);
    if (res.opt_sol > b) {
      printf("YES!\n");
      return INFEASIBLE_SYS;
    } else {
      printf("NO!\n");
    }

    printf("4c. Checking if i == (n-1): %d == %d ? ", i, n - 1);
    if (i == n - 1) {
      if (res.opt_sol == NO_SOL) {
        printf("Error: feasible but no solution found!\n");
      }
      printf("YES\n");
      res.feasible = FEASIBLE_SYS;
      return FEASIBLE_SYS;
    } else {
      printf("NO\n");
    }

    for (j = i + 1; j < n; j++) {
      k = pi[j];
      d = ceil((res.opt_sol + alpha[k]) / cur_ts.per[k]) - x[k];
      printf("5. Adding d = ceil ((%.3f + %.3f) / %.3f) - %.3f = %.3f\n",
             res.opt_sol, alpha[k], cur_ts.per[k], x[k], d);
      x[k] += d;
      y[k] += cur_ts.per[k] * d;
      r += cur_ts.wcet[k] * d;
    }

    i = i + 1;
    printf("5b. New values:\n\ti=%d\n\tx={", i);

    for (h = 0; h < n; h++) {
      printf("%.3f ", x[h]);
    }
    printf("}\n\ty={");
    for (h = 0; h < n; h++) {
      printf("%.3f ", y[h]);
    }
    printf("}\n");

    h1 = 0, h2 = 0;
    /*Copy pi in pi1, pi2*/
    for (h = 0; h < n; h++) {
      if (h < i) {
        pi1[h1] = pi[h];
        h1++;
      } else {
        pi2[h2] = pi[h];
        h2++;
      }
    }

    /*
    printf("6. Creating pi1={");
    for (h = 0; h < h1; h++) {
      printf("%d ", pi1[h]);
    }
    printf("}, pi2={");
    for (h = 0; h < h2; h++) {
      printf("%d ", pi2[h]);
    }
    printf("}\n");
    */

    /* Print array pi*/
    printf("6. Sorting pi in the range %d to %d : [", i, n - 1);
    for (h = 0; h < n; h++) {
      printf("%d ", pi[h]);
    }
    printf("] -> [");

    /* Sort pi2 in a nonincreasing order using the key y[pi2[h]]*/
    for (h = 0; h < n - i; h++) {
      int tmp = pi2[h];
      int prev_h = h - 1;

      while (prev_h >= 0 && y[pi2[prev_h]] < y[tmp]) {
        pi2[prev_h + 1] = pi2[prev_h];
        --prev_h;
      }
      pi2[prev_h + 1] = tmp;
    }

    /* Print array pi2*/
    for (h = 0; h < n - i; h++) {
      printf("%d ", pi2[h]);
    }
    printf("]\n");

    printf("7. Merging with i=%d and n=%d, ", i, n);

    merge_pi(n, i, pi1, pi2, pi);

    printf("pi1 and pi2, result in y: [");
    for (h = 0; h < n; h++) {
      printf("%.3f ", y[pi[h]]);
    }
    printf("]\n");

    iter++;
  }

  printf("Exited while loop with iter = %d\n", iter);
}

void cp_kern_res_set_zero(cp_kern_res_t *res) {
  res->feasible = UNDET_SYS;
  res->opt_sol = NO_SOL; /* Valore non sicuro?*/
}

/* Kernel's parameters for EDF schedulability test:
 * - n = k
 * - alpha =  D̂ - T
 * - beta = 1
 * - a = - b_k + 1
 * - b = - a_k
 * */
int cp_kern_edf(int k) {
  double a, b;
  int i = 0, h = 0;
  int h1, h2;
  int j, iter = 0;
  int i_start; /*smallest point of the domain of f*/
  int n = k + 1;
  double p, q, r, d;

  a = -B[k] + 1;
  b = -A[k];

  printf("\n --- Initializing parameters in CP-KERN over [a,b]=[%f,%f]----\n",
         a, b);

  if (a > b) {
    res.feasible = INFEASIBLE_SYS;
    return INFEASIBLE_SYS;
  }

  /*None tasks*/
  if (n == 1) {
    res.feasible = FEASIBLE_SYS;
    res.opt_sol = a;
    return FEASIBLE_SYS;
  }

  /*Initializing p,q and r*/
  p = r = beta;
  q = 1;

  /* Initializing x, y and pi*/
  for (h = 0; h < n; h++) {
    x[h] = ceil((a + alpha[h]) / cur_ts.per[h]);
    if (x[h] == -0.0) {
      x[h] = 0.0;
    }
    printf("x_%d = ceil((%.2f + %.2f) / %.2f) = %.2f\n", i, a, alpha[h],
           cur_ts.per[h], x[h]);

    y[h] = cur_ts.per[h] * x[h] - alpha[h];
    printf("y_% = %.2f * %.2f - %.2f = %.2f\n", i, cur_ts.per[h], x[h],
           alpha[h], y[h]);
    pi[h] = h; /* pi stores the indexes of y so that for any i <= j y[pi[i]]
                      <= y[pi[j]], right now is the identity permutation*/
    p += (cur_ts.wcet[h] * alpha[h]) / cur_ts.per[h];
    q -= (cur_ts.wcet[h] / cur_ts.per[h]);
    r += cur_ts.wcet[h] * x[h];
  }

  print_kernel_params(n, p, q, r);

  /**
   * Checking infeasibility condition that states that if
   *
   *  (p) = ( beta + sum of (U_i * alpha_i) ) > 0
   *
   *                 and
   *
   *  (q) = ( 1 - sum of U_i ) == 0
   *
   *  then the IP-KERN is infeasible.
   *
   * */

  printf(
      "\n--- Checking infeasibility condition && initializing i_start ---\n");
  if (p > 0 && q == 0) {
    res.feasible = INFEASIBLE_SYS;
    return INFEASIBLE_SYS;
  }

  if (q == 0) {
    i_start = 1;
  } else {
    i_start = 0;
  }
  printf("i_start = %d\n", i_start);

  printf("\n--- Sorting pi in a nonincreasing order of y[pi[i]] ---\n");

  printf("pi={");
  for (h = 0; h < n; h++) {
    printf("%d ", pi[h]);
  }
  printf("} -> {");

  /* Sorting pi in a nonincreasing order using the key y[pi[h]]*/
  for (h = 1; h < n; h++) {
    int tmp = pi[h];
    int prev_h = h - 1;

    while (prev_h >= 0 && y[pi[prev_h]] < y[tmp]) {
      pi[prev_h + 1] = pi[prev_h];
      --prev_h;
    }
    pi[prev_h + 1] = tmp;
  }

  for (h = 0; h < n; h++) {
    printf("%d ", pi[h]);
  }
  printf("}\n");

  printf(" -> y[pi[i]]={");
  for (h = 0; h < n; h++) {
    printf("%f ", y[pi[h]]);
  }
  printf("}\n");

  while (iter < 2) {
    p = res.opt_sol = r;
    q = 1.0;
    i = n - 1;

    /* Solve LP*/
    printf("\n--- Computing optimal value for LP ---");

    printf("\n1. Setting p=%f, q=%f, i=%d.\n", p, q, i);
    while (i >= i_start) {
      k = pi[i];
      printf("2. k=%d:\n", k);
      printf("\tTask %d with C=%.3f, T=%.3f, alpha=%.3f\n", k, cur_ts.wcet[k],
             cur_ts.per[k], alpha[k]);
      printf("\t%.3f <= %.3f*%.3f=%.3f ?", p, q, y[k], q * y[k]);
      if (p <= q * y[k]) {
        printf(" YES! Possible solution =%.3f\n", p / q);
        res.opt_sol = p / q;
        break;
      } else {
        printf("\tNO\n");
      }

      p -= cur_ts.wcet[k] * x[k] + cur_ts.wcet[k] * alpha[k] / cur_ts.per[k];
      q -= cur_ts.wcet[k] / cur_ts.per[k];
      i--;

      printf("\tUpdated p=%f, q=%f, i=%d\n", p, q, i);
    }

    printf("3. i=%d, i_start=%d: (i==i_start)? ", i, i_start);
    if (i == i_start - 1) {
      res.opt_sol = p / q;
      printf("YES, possibile solution: %f.\n", res.opt_sol);
    } else {
      printf("NO\n");
    }

    if (res.opt_sol <= a) {
      res.opt_sol = a;
      res.feasible = FEASIBLE_SYS;
      printf("4a. Solution found (res.opt_sol <= a): %f\n", a);
      return FEASIBLE_SYS;
    }

    printf("4b. Checking if value is greater than b: %f > %f ? ", res.opt_sol,
           b);
    if (res.opt_sol > b) {
      printf("YES!\n");
      return INFEASIBLE_SYS;
    } else {
      printf("NO!\n");
    }

    printf("4c. Checking if i == (n-1): %d == %d ? ", i, n - 1);
    if (i == n - 1) {
      if (res.opt_sol == NO_SOL) {
        printf("Error: feasible but no solution found!\n");
      }
      printf("YES\n");
      res.feasible = FEASIBLE_SYS;
      return FEASIBLE_SYS;
    } else {
      printf("NO\n");
    }

    for (j = i + 1; j < n; j++) {
      k = pi[j];
      d = ceil((res.opt_sol + alpha[k]) / cur_ts.per[k]) - x[k];
      printf("5. Adding d = ceil ((%.3f + %.3f) / %.3f) - %.3f = %.3f\n",
             res.opt_sol, alpha[k], cur_ts.per[k], x[k], d);
      x[k] += d;
      y[k] += cur_ts.per[k] * d;
      r += cur_ts.wcet[k] * d;
    }

    i = i + 1;
    printf("5b. New values:\n\ti=%d\n\tx={", i);

    for (h = 0; h < n; h++) {
      printf("%.3f ", x[h]);
    }
    printf("}\n\ty={");
    for (h = 0; h < n; h++) {
      printf("%.3f ", y[h]);
    }
    printf("}\n");

    h1 = 0, h2 = 0;
    /*Copy pi in pi1, pi2*/
    for (h = 0; h < n; h++) {
      if (h < i) {
        pi1[h1] = pi[h];
        h1++;
      } else {
        pi2[h2] = pi[h];
        h2++;
      }
    }

    /*
    printf("6. Creating pi1={");
    for (h = 0; h < h1; h++) {
      printf("%d ", pi1[h]);
    }
    printf("}, pi2={");
    for (h = 0; h < h2; h++) {
      printf("%d ", pi2[h]);
    }
    printf("}\n");
    */

    /* Print array pi*/
    printf("6. Sorting pi in the range %d to %d : [", i, n - 1);
    for (h = 0; h < n; h++) {
      printf("%d ", pi[h]);
    }
    printf("] -> [");

    /* Sort pi2 in a nonincreasing order using the key y[pi2[h]]*/
    for (h = 0; h < n - i; h++) {
      int tmp = pi2[h];
      int prev_h = h - 1;

      while (prev_h >= 0 && y[pi2[prev_h]] < y[tmp]) {
        pi2[prev_h + 1] = pi2[prev_h];
        --prev_h;
      }
      pi2[prev_h + 1] = tmp;
    }

    /* Print array pi2*/
    for (h = 0; h < n - i; h++) {
      printf("%d ", pi2[h]);
    }
    printf("]\n");

    printf("7. Merging with i=%d and n=%d, ", i, n);

    merge_pi(n, i, pi1, pi2, pi);

    printf("pi1 and pi2, result in y: [");
    for (h = 0; h < n; h++) {
      printf("%.3f ", y[pi[h]]);
    }
    printf("]\n");

    iter++;
  }

  printf("Exited while loop with iter = %d\n", iter);
}

void cp_kern_alloc(int n) {

  x = malloc(n * sizeof(double));
  y = malloc(n * sizeof(double));
  pi = malloc(n * sizeof(int));
  pi1 = malloc(n * sizeof(int));
  pi2 = malloc(n * sizeof(int));
}

void cp_kern_free() {
  ts_free(&cur_ts);
  free(alpha);
  free(A);
  free(B);
  free(x);
  free(y);
  free(pi);
  free(pi1);
  free(pi2);
}

void print_kernel_params(int n, double p, double q, double r) {
  int h;
  printf("alpha={");
  for (h = 0; h < cur_ts.num; h++) {
    printf("%f ", alpha[h]);
  }
  printf("}\nx={");
  for (h = 0; h < n; h++) {
    printf("%f ", x[h]);
  }
  printf("}\ny={");
  for (h = 0; h < n; h++) {
    printf("%f ", y[h]);
  }
  printf("}\npi={");
  for (h = 0; h < n; h++) {
    printf("%d ", pi[h]);
  }
  printf("}\np=%f, q=%f, r=%f, n=%d\n", p, q, r, n);
}

void merge_pi(int n, int i, int *pi1, int *pi2, int *pi) {

  int h1 = 0, h2 = 0, h = 0;

  while (h1 < i && h2 < n - i && h < n) {
    if (y[pi1[h1]] > y[pi2[h2]]) {
      pi[h] = pi1[h1];
      h1++;
    } else {
      pi[h] = pi2[h2];
      h2++;
    }
    h++;
  }

  while (h1 < i) {
    pi[h] = pi1[h1];
    h1++;
    h++;
  }

  while (h2 < n - i) {
    pi[h] = pi2[h2];
    h2++;
    h++;
  }
}
