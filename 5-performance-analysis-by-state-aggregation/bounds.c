#include<stdio.h>
#include<stdlib.h>
#include<string.h>

#define K 2     // number of job types



/*** VARIABLES ***/

int S;          // number of servers
int m;          // number of customer classes of each type
int r[K];       // degree of each type
long double varrho[K];  // traffic intensity of each type
long double *p[K];      // probability that a given server cannot serve
                        // a given number of jobs of each type

long double **pi;       // unnormalized stationary measure, with pi[0][0] = 1.
long double **piL0;     // unnormalized product pi L_0
long double **piL1;     // unnormalized product pi L_1




/*** COMPUTATIONS ***/

void compute_mean (long double *L, long double *psi) {
  int i, j, k;
  long double rate;

  // initialization
  for (i = 0 ; i <= m ; ++i) {
    for (j = 0 ; j <= m ; ++j) {
      pi[i][j] = 0.;
      piL0[i][j] = 0.;
      piL1[i][j] = 0.;
    }
  }
  pi[0][0] = 1.;

  (*psi) = 0.;
  for (k = 0 ; k < K ; ++k) L[k] = 0.;

  // recursion
  for (i = 0 ; i <= m ; ++i) {
    for (j = 0 ; j <= m ; ++j) {
      rate = S * (1. - p[0][i] * p[1][j]) - i * varrho[0] - j * varrho[1];

      if (i > 0) pi[i][j] += (m - i + 1) * varrho[0] * pi[i-1][j];
      if (j > 0) pi[i][j] += (m - j + 1) * varrho[1] * pi[i][j-1];
      if (i + j > 0) pi[i][j] /= rate;
      (*psi) += pi[i][j];

      if (i > 0) {
        piL0[i][j] = i * varrho[0] * pi[i][j];
        piL0[i][j] += (m - i + 1) * varrho[0] * (pi[i-1][j] + piL0[i-1][j]);
        if (j > 0) piL0[i][j] += (m - j + 1) * varrho[1] * piL0[i][j-1];
        piL0[i][j] /= rate;
        L[0] += piL0[i][j];
      }

      if (j > 0) {
        piL1[i][j] = j * varrho[1] * pi[i][j];
        if (i > 0) piL1[i][j] += (m - i + 1) * varrho[0] * piL1[i-1][j];
        piL1[i][j] += (m - j + 1) * varrho[1] * (pi[i][j-1] + piL1[i][j-1]);
        piL1[i][j] /= rate;
        L[1] += piL1[i][j];
      }
    }
  }

  // normalization
  (*psi) = 1. / (*psi);
  L[0] /= m;
  L[1] /= m;
}



int main(int argc, char **argv) {
  int i, k;
  long double alpha, h, epsilon,
         *L_plus, *L_minus,
         psi_plus, psi_minus;
  char file_name[500], type[10];
  FILE *file = NULL;

  // read parameters
  if (argc < 6) {
    fprintf(stderr, "Error: you should have at least 5 arguments.\n");
    exit(0);
  }
  S = atoi(argv[1]);
  m = atoi(argv[2]);

  for (k = 0 ; k < K ; ++k) {
    r[k] = atoi(argv[3+k]);
    p[k] = malloc((m + 1) * sizeof(long double));
    p[k][0] = 1.;
    for (i = 1 ; i <= m ; ++i) {
      p[k][i] = p[k][i-1] * (1. - 1. * r[k] / S);
    }
  }
  h = S * (1. - p[0][m] * p[1][m]);

  epsilon = atof(argv[5]);

  // initialize internal variables
  pi = malloc((m + 1) * sizeof(long double*));
  piL0 = malloc((m + 1) * sizeof(long double*));
  piL1 = malloc((m + 1) * sizeof(long double*));
  for (i = 0 ; i <= m ; ++i) {
    pi[i] = malloc((m + 1) * sizeof(long double));
    piL0[i] = malloc((m + 1) * sizeof(long double));
    piL1[i] = malloc((m + 1) * sizeof(long double));
  }

  // initialize outputs
  L_plus = malloc(K * sizeof(long double)); 
  L_minus = malloc(K * sizeof(long double)); 

  // print stdout
  printf("Bounds - K = %d customer types - epsilon = %Le\n", K, epsilon);
  printf("%d servers - %d classes per type\n", S, m);
  for (k = 0 ; k < K ; ++k) {
    printf("Degree of type-%d customers: %d\n", k+1, r[k]);
  }

  // create the file
  sprintf(file_name, "data/bounds-%s-%s-%s-%s-%s.csv",
      argv[1], argv[2], argv[3], argv[4], argv[5]);
  file = fopen(file_name, "w");

  // write headers in the file
  fprintf(file, "alpha");
  for (k = 0 ; k < K ; ++k) {
    fprintf(file, ",upper%d,lower%d", k+1, k+1);
  }
  fprintf(file, "\n");

  // computations
  for (alpha = .005 ; alpha < .98 ; alpha += .005) {

    // rho / (1 + epsilon)
    varrho[0] = alpha * r[0] * h / ((1. + epsilon) * (r[0] + r[1]) * m);
    varrho[1] = alpha * r[1] * h / ((1. + epsilon) * (r[0] + r[1]) * m);
    compute_mean(L_plus, &psi_plus);

    // rho / (1 - epsilon)
    varrho[0] = alpha * r[0] * h / ((1. - epsilon) * (r[0] + r[1]) * m);
    varrho[1] = alpha * r[1] * h / ((1. - epsilon) * (r[0] + r[1]) * m);
    compute_mean(L_minus, &psi_minus);

    // rho
    varrho[0] = alpha * r[0] * h / ((r[0] + r[1]) * m);
    varrho[1] = alpha * r[1] * h / ((r[0] + r[1]) * m);

    // write in the file
    fprintf(file, "%Le", alpha);
    for (k = 0 ; k < K ; ++k) {
      fprintf(file, ",%Le,%Le",
          varrho[k] / (psi_minus * L_plus[k]),
          varrho[k] / (psi_plus * L_minus[k]));
    }
    fprintf(file, "\n");

  }

  // narrow data points near the critical load rho = 1
  for (alpha = .98 ; alpha < 1. - epsilon ; alpha += .0001) {

    // rho / (1 + epsilon)
    varrho[0] = alpha * r[0] * h / ((1. + epsilon) * (r[0] + r[1]) * m);
    varrho[1] = alpha * r[1] * h / ((1. + epsilon) * (r[0] + r[1]) * m);
    compute_mean(L_plus, &psi_plus);

    // rho / (1 - epsilon)
    varrho[0] = alpha * r[0] * h / ((1. - epsilon) * (r[0] + r[1]) * m);
    varrho[1] = alpha * r[1] * h / ((1. - epsilon) * (r[0] + r[1]) * m);
    compute_mean(L_minus, &psi_minus);

    // rho
    varrho[0] = alpha * r[0] * h / ((r[0] + r[1]) * m);
    varrho[1] = alpha * r[1] * h / ((r[0] + r[1]) * m);

    // write in the file
    fprintf(file, "%Le", alpha);
    for (k = 0 ; k < K ; ++k) {
      fprintf(file, ",%Le,%Le",
          varrho[k] / (psi_minus * L_plus[k]),
          varrho[k] / (psi_plus * L_minus[k]));
    }
    fprintf(file, "\n");

  }

  printf("\n");

  fclose(file);

  free(L_plus); free(L_minus);
  free(pi); free(piL0); free(piL1);

  return 0;
}
