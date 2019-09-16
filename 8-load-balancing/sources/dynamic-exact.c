#include "utils.h"
// defines the constants MAX_I and MAX_K

#define ONE 1UL
#define INT_LENGTH 64   /* length of the unsigned integers used to encode the state */






/*** VARIABLES ***/



/* Inputs */

uint64_t I;       /* number of computers */
uint64_t K;       /* number of job types */

uint64_t elli[MAX_I];   /* number of tokens of each computer */
uint64_t compatible[MAX_I][MAX_K]; /* compatiblity graph between computers and types,
                                      encoded as an adjacency matrix */
uint64_t notcompatible[MAX_K][MAX_I + 1]; /* list of computers each type
                                             is **not** compatible with */

double mui[MAX_I];      /* normalized service capacity of each computer */
double mu;              /* overall capacity of the cluster */
double nuk[MAX_K];      /* normalized arrival rate of each job type */
double *nuA;            /* normalized arrival rate function */



/* Variables for the computations */

double Rho;       /* largest load to consider */
double delta;     /* step by which we increase the load */

uint64_t length;  /* number of bits to encode the number of available tokens of each computer */
uint64_t mask;    /* mask to detect the number of available tokens of a given computer */
uint64_t ymax;    /* maximum aggregate state to consider */
uint64_t Amax;    /* maximum active set to consider */
uint64_t *ei;     /* I-dimensional unit vectors */
double *piy;      /* unnormalized stationary measure, with piy[0] = 1 */
double pi;        /* normalizating constant of the stationary distribution */



/* Outputs */

char name[40];    /* name of the output file */





/*** UTILITARY ***/

/* print an unsigned integer in binary */
void printb (uint64_t U) {
  uint64_t i;

  for (i = 0 ; i < INT_LENGTH ; ++i) {
    if ( !((INT_LENGTH - i) % length) ) printf(" ");
    printf("%lu", U >> (INT_LENGTH - 1 - i) & ONE);
  }
}

/* verifies whether the bitmap y is valid */
static inline int is_valid (uint64_t y) {
  uint64_t i;

  for (i = 0 ; i < I ; ++i) {
    if ( ((y >> (length * i)) & mask) > elli[i] )
      return 0;
  }

  return 1;
}



/*** COMPUTATIONS ***/

/* basic operations
 * access the number of available tokens from computer i:
 *   y >> (length * i) & mask
 * check whether there is at least one available token from computer i:
 *   ((y >> (length * i)) & mask) > 0
 * remove one available token of computer i:
 *   y -= (ONE << (length * i))
 */

void stationary (double rho) {
  uint64_t i, y, A;

  /* initialization */
  piy[0] = 1.;
  pi = 1.;

  /* recursion step */
  for (y = ONE ; y < ymax ; ++y) {
    piy[y] = 0.;

    if (is_valid(y)) {
      A = 0;
      for (i = 0 ; i < I ; ++i) {
        if ( (y >> (length * i)) & mask ) {
          /* there is at least one available token from computer i */
          piy[y] += mui[i] * piy[y - ei[i]];
          A += ONE << i;
        }
      }
      piy[y] /= (nuA[A] * rho);

      pi += piy[y];
    }
  }
}

double loss (int k) {
  uint64_t i, t, y, carry;
  double betak;

  betak = 0.;
  y = 0;
  carry = 0;

  while (!carry) {
    /* update the loss probability */
    betak += piy[y];

    /* update y */
    carry = 1;
    t = 0;
    while ( carry && ((i = notcompatible[k][t]) != -1) ) {
      if ( ((y >> (i * length)) & mask) < elli[i] ) {
        y += ONE << (i * length);
        carry = 0;
      } else y -= elli[i] << (i * length);
      ++t;
    }
  }

  return betak / pi;
}

double idling (int i) {
  uint64_t j, y, carry;
  double psii;

  psii = 0.;
  y = elli[i] << (i * length);
  carry = 0;

  while (!carry) {
    /* update the idling probability */
    psii += piy[y];

    /* update y */
    carry = 1;
    j = 0;
    while ( carry && (j < I) ) {
      if (j != i) {
        if ( ((y >> (j * length)) & mask) < elli[j] ) {
          y += ONE << (j * length);
          carry = 0;
        } else y -= elli[j] << (j * length);
      }
      ++j;
    }
  }

   return psii / pi;
}

double mean (int i) {
  uint64_t y;
  double Li;

  Li = 0.;
  for (y = ONE ; y < ymax ; ++y) {
    if (is_valid(y)) {
      Li += ((y >> (length * i)) & mask) * piy[y];
    }
  }

  return elli[i] - (Li / pi);
}





/*** INPUTS AND OUTPUTS ***/

void read_inputs (int argc, char **argv) {
  int c;
  uint64_t i, k, t, current, total;
  char *arg;
  double norm;

  /* default initialization */
  sprintf(name, "toto");
  K = 0; I = 0;
  for (i = 0 ; i < MAX_I ; ++i) {
    elli[i] = 0;
    mui[i] = 0.;
    for (k = 0 ; k < MAX_K ; ++k)
      compatible[i][k] = 0;
  }
  for (k = 0 ; k < MAX_K ; ++k) {
    nuk[k] = 0.;
    notcompatible[k][0] = -1;
  }
  Rho = 5.;
  delta = .01;

  /* read arguments */
  opterr = 0;
  while ( (c = getopt(argc, argv, "f:k:i:l:c:r:m:d:h")) != -1 ) {
    switch(c) {
      case 'f':
        /* output file name */
        strncpy(name, optarg, 40);
        break;
      case 'k':
        /* number of job types */
        K = atoi(optarg);
        if (K > MAX_K) {
          fprintf(stderr, "Option -k: the number of job types"
              "cannot be larger than %d.\n", MAX_K);
          exit(EXIT_FAILURE);
        }
        break;
      case 'i':
        /* number of computers */
        I = atoi(optarg);
        if (I > MAX_I) {
          fprintf(stderr, "Option -i: the number of computers"
              "cannot be larger than %d.\n", MAX_I);
          exit(EXIT_FAILURE);
        }
        break;
      case 'l':
        /* number of tokens of each computer */
        total = 0;
        arg = strtok(optarg, " ");
        i = 0;
        while (i < I && arg != NULL) {
          elli[i] = atoi(arg);
          total += elli[i];
          ++i;
          arg = strtok(NULL, " ");
        }

        if (i == 0) {
          fprintf(stderr, "Option -l: the numbers of tokens of each computer.\n");
          exit(EXIT_FAILURE);
        }

        while (i < I) {
          elli[i] = elli[i-1];
          total += elli[i];
          ++i;
        }

        /* number of bits to encode the number of tokens of each computer */
        length = 0;
        for (i = 0 ; i < I ; ++i) {
          current = INT_LENGTH;
          while ( (current > length) && !(elli[i] >> (current - 1) & ONE) ) --current;
          if (current > length) length = current;
        }

        if (I * length >= INT_LENGTH) {
          fprintf(stderr, "Option -l: the unsigned integer type does not have enough bits to encode all states.\n");
          exit(EXIT_FAILURE);
        }

        /* other lengths */
        mask = (ONE << length) - ONE;
        ymax = ONE << (I * length);
        Amax = ONE << I;

        break;
      case 'c':
        /* computer-to-type compatibilities,
         * encoded as an adjacency matrix between computers and types */
        arg = strtok(optarg, " ");

        for (i = 0 ; i < I ; ++i) {
          for (k = 0 ; k < K ; ++k) {
            if (arg == NULL) {
              fprintf(stderr, "Option -c: you didn't specify all the compatibilities.\n");
              exit(EXIT_FAILURE);
            }

            compatible[i][k] = atoi(arg);
            arg = strtok(NULL, " ");
          }
        }
        break;
      case 'r':
        arg = strtok(optarg, " ");

        /* arrival rates */
        norm = 0.;
        for (k = 0 ; k < K ; ++k) {
          if (arg == NULL) {
            fprintf(stderr, "Option -r: you didn't specify all the arrival rates.\n");
            exit(EXIT_FAILURE);
          }

          nuk[k] = atof(arg);
          norm += nuk[k];
          arg = strtok(NULL, " ");
        }

        /* normalize the arrival rates */
        for (k = 0 ; k < K ; ++k) nuk[k] /= norm;

        /* service rates */
        mu = 0.;
        for (i = 0 ; i < I ; ++i) {
          if (arg == NULL) {
            fprintf(stderr, "Option -r: you didn't specify all the service rates.\n");
            exit(EXIT_FAILURE);
          }

          mui[i] = atof(arg);
          mu += mui[i];
          arg = strtok(NULL, " ");
        }

        /* normalize the service rates */
        for (i = 0 ; i < I ; ++i) mui[i] /= mu;

        break;
      case 'm':
        /* maximum load that is considered */
        Rho = atof(optarg);
        break;
      case 'd':
        /* step by which we increase the load */
        delta = atof(optarg);
        break;
      case 'h':
        /* help */
        printf("\nOptions:\n"
            "-f: output file name\n"
            "-k: number of job types\n"
            "-i: number of computers\n"
            "-l: number of tokens of each class\n"
            "-c: computer-to-type compatibilities, encoded as an adjacency matrix\n"
            "-r: per-type arrival rates and per-computer capacities\n"
            "-m: maximum load to consider\n"
            "    (default value: 5)\n"
            "-d: step by which the load is increased\n"
            "    (default value: 0.01)\n"
            "\n");
        exit(0);
        break;
      case '?':
        /* error */
        if (optopt == 'c')
          fprintf(stderr, "Option -%c requires an argument.\n", optopt);
        else if (isprint(optopt))
          fprintf(stderr, "Unknown option `-%c'.\n", optopt);
        else
          fprintf(stderr, "Unknown option character `\\x%x'.\n", optopt);
        exit(EXIT_FAILURE);
      default:
        abort();
    }
}

  for (i = optind; i < argc; ++i)
    printf ("Non-option argument %s\n", argv[i]);
}

void print_inputs () {
  uint64_t i, k;

  printf("\n--- DYNAMIC - EXACT ---\n\n");

  /* output file prefix name */
  printf("Output file prefix name: %s.csv\n\n", name);

  /* computation parameters */
  printf("Maximum load to be considered:      %.3e\n", Rho);
  printf("Step by which we increase the load: %.3e\n\n", delta);

  /* numbers of types and computers */
  printf("%lu job type(s) and %lu computer(s)\n\n", K, I);

  /* numbers of tokens of each computer */
  for (i = 0 ; i < I ; ++i) {
    printf("Computer %lu: %lu token(s)\n", i, elli[i]);
  }
  printf("\n");

  /* computer-to-type compatibilities */
  printf("Adjacency matrix of the compatiblity graph between computers and types:\n");
  for (i = 0 ; i < I ; ++i) {
    for (k = 0 ; k < K ; ++k) printf("%lu ", compatible[i][k]);
    printf("\n");
  }
  printf("\n");

  /* arrival rates */
  printf("Normalized arrival rates:\n");
  for (k = 0 ; k < K ; ++k) printf("%.3e  ", nuk[k]);
  printf("\n\n");

  /* service rates */
  printf("Normalized service rates:\n");
  for (i = 0 ; i < I ; ++i) printf("%.3e  ", mui[i]);
  printf("\nTotal service rate: %.3e\n\n", mu);
}

FILE *open_csv_file () {
  char file_name[44];
  FILE *file = NULL;

  sprintf(file_name, "%s.csv", name);
  file = fopen(file_name, "w");

  if (file == NULL) {
    fprintf(stderr, "Couldn't open the output file %s.\n\n", file_name);
    exit(EXIT_FAILURE);
  } else printf("Save data in the file '%s'.\n\n", file_name);

  return file;
}

void initialize_arrays () {
  uint64_t i, k, t, A;

  /* stationary measure */
  piy = malloc(ymax * sizeof(double));
  if (piy == NULL) {
    fprintf(stderr, "Couldn't allocate the variable piy.\n\n");
    exit(EXIT_FAILURE);
  }

  /* unit vectors */
  ei = malloc(I * sizeof(uint64_t));
  if (ei == NULL) {
    fprintf(stderr, "Couldn't allocate the variable ei.\n\n");
    exit(EXIT_FAILURE);
  }

  for (i = 0 ; i < I ; ++i) ei[i] = ONE << (length * i);

  /* rate function */
  nuA = malloc(Amax * sizeof(double));
  if (nuA == NULL) {
    fprintf(stderr, "Couldn't allocate the variable nuA.\n\n");
    exit(EXIT_FAILURE);
  }

  nuA[0] = 0.;
  for (A = ONE ; A < Amax ; ++A) {
    nuA[A] = 0.;
    for (k = 0 ; k < K ; ++k) {
      for (i = 0 ; i < I ; ++i) {
        if ( ((A >> i) & ONE) && compatible[i][k] ) {
          nuA[A] += nuk[k];
          break;
        }
      }
    }
  }

  /* list of adjacency lists */
  for (k = 0 ; k < K ; ++k) {
    t = 0;
    for (i = 0 ; i < I ; ++i) {
      if (!compatible[i][k]) {
        notcompatible[k][t] = i;
        ++t;
      }
    }
    notcompatible[k][t] = -1;
  }
}





/*** MAIN ***/

int main (int argc, char **argv) {
  uint64_t i, k;
  double rho, betak, beta, etai, eta, Li, L;
  FILE *file;

  /* input arguments */
  read_inputs(argc, argv);
  print_inputs();
  initialize_arrays();

  /* open the output file */
  file = open_csv_file();

  /* print the column names */
  fprintf(file, "rho");
  for (k = 0 ; k < K ; ++k) {
    fprintf(file, ",betak%lu", k+1);
  }
  for (i = 0 ; i < I ; ++i) {
    fprintf(file, ",etai%lu,Li%lu,gammai%lu", i+1, i+1, i+1);
  }
  fprintf(file, ",beta,eta,L,gamma\n");

  for (rho = delta ; rho < Rho + .5 * delta ; rho += delta) {
    /* computations */
    printf("rho = %.3e\n", rho);
    stationary(rho);

    /* print data in the output file */
    fprintf(file, "%e", rho);

    /* per-type loss probability */
    beta = 0.;
    for (k = 0 ; k < K ; ++k) {
      betak = loss(k);
      fprintf(file, ",%e", betak);
      beta += nuk[k] * betak;
    }

    /* per-computer idling probability, mean number of jobs, and mean service rate */
    eta = 0.;
    L = 0.;
    for (i = 0 ; i < I ; ++i) {
      etai = 1. - idling(i);
      eta += mui[i] * etai;

      Li = mean(i);
      L += Li;

      if (Li > 0.) fprintf(file, ",%e,%e,%e", etai, Li, mu * mui[i] * etai / Li);
      else fprintf(file, ",%e,%e,%e", etai, 0., 0.);
    }

    /* system-wide metrics */
    fprintf(file, ",%e,%e,%e,%e\n",
        beta, eta, L, rho * mu * (1. - beta) / L);
  }

  /* close the file and release memory */
  fclose(file);
  free(piy); free(ei); free(nuA);
}
