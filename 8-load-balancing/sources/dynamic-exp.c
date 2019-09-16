#include "utils.h"
// defines the constants MAX_I, MAX_K, MAX_R, and MAX_L





/*** GLOBAL VARIABLES ***/



/* Inputs */

int K;                  /* number of job types */
int I;                  /* number of computers */

int elli[MAX_I];        /* number of tokens at each computer */
int Ki[MAX_I][MAX_K + 1];   /* adjacency list of the compatibility
                               graph between computers and types */

double mui[MAX_I];      /* normalized service capacity of each computer */
double mu;              /* overall capacity of the cluster */
double nuk[MAX_K];      /* normalized arrival rate of each job type */
double event_rate;      /* event rate */
double threshold;       /* probability that the next event is an arrival
                         rather than a departure */



/* Variables of the simulation */

double Rho;             /* largest load to consider */
double delta;           /* step by which we increase the load */
int R;                  /* number of simulation runs for each load */
double alpha = 1.96;    /* P(-alpha < X < alpha) = 95% with X ~ N(0,1) */

int warmup;             /* number of warm-up steps */
int steady;             /* number of steps in steady state */

int xi[MAX_I];          /* number of jobs at each computer */
int xk[MAX_K];          /* number of jobs of each type */
int xik[MAX_I][MAX_K];  /* number of jobs of each type at each computer */



/* Circular buffer encoding the sequence of available tokens */

int buffer[MAX_L];      /* circular buffer that contains the available tokens */
int head;               /* position of the buffer head */
int tail;               /* position of the buffer tail */

int assignment[MAX_K];
/* current assignment of each job type,
 * i.e., the computer an incoming job of each type would be assigned to */



/* Outputs */

char name[40];        /* name of the output file */

int timeout[MAX_I];   /* cumulative number of times the timer of each computer expires */
int empty[MAX_I];     /* cumulative number of times the timer of each computer expires
                         while the computer is empty */
double etai[MAX_I];   /* empirical activity probability of each computer */

int arrival[MAX_K];   /* cumulative number of arrivals within each type */
int loss[MAX_K];      /* cumulative number of lost jobs within each type */
double betak[MAX_K];  /* empirical loss probability within each type */

double T;             /* total time in steady state */
double Lk[MAX_K];     /* empirical expected number of jobs of each type */
double Li[MAX_I];     /* empirical expected number of jobs at each computer */





/*** PROTOTYPES ***/

/* Queue of available tokens */

void initialize_queue ();   /* initialize the queue */
void print_queue ();        /* print the queue and the future assignment decisions */

/* Simulations */

void simulation ();

/* Inputs and outputs */

void read_inputs (int argc, char **argv); /* read the input arguments, given as options */
void print_inputs ();       /* print all input arguments */
FILE *open_csv_file ();     /* open the output data file */





/*** QUEUE OF AVAILABLE TOKENS ***/

void initialize_queue () {
  int i, k, t;

  /* start with an empty queue */
  head = 0; tail = 0;
  for (k = 0 ; k < K ; ++k) assignment[k] = -1;

  /* progressively append tokens */
  for (i = 0 ; i < I ; ++i) {
    if (elli[i] > 0) {
      /* append all computer-i tokens to the queue */
      for (t = 0 ; t < elli[i] ; ++t) {
        buffer[tail] = i;
        tail = (tail + 1) % MAX_L;
      }

      /* assign the job types */
      t = 0;
      while ( (k = Ki[i][t]) != -1 ) {
        if (assignment[k] == -1) assignment[k] = i;
        ++t;
      }
    }
  }
}

void print_queue () {
  int k, p;

  printf("Buffer: [ ");
  for (p = head ; p != tail ; p = (p + 1) % MAX_L)
    printf("%d ", buffer[p]);
  printf("]\n\n");

  printf("Assignments:\n");
  for (k = 0 ; k < K ; ++k)
    printf("Incoming jobs of type %d would be assigned to computer %d\n", k, assignment[k]);
  printf("\n");
}

static inline void append_to_queue (int i) {
  int k, p, t;

  /* append the computer token */
  buffer[tail] = i;
  tail = (tail + 1) % MAX_L;

  /* update the assignments */
  t = 0;
  while ( (k = Ki[i][t]) != -1 ) {
    if (assignment[k] == -1) assignment[k] = i;
    ++t;
  }
}

static inline void remove_from_queue (int i) {
  int k, p, q, t, current, next, nb_released;

  /* move the tokens until we find a token of computer i */
  p = head;
  q = (head = (head + 1) % MAX_L);
  current = buffer[p];
  next = buffer[q];
  while (current != i) {
    p = q;
    q = (q + 1) % MAX_L;
    buffer[p] = current;
    current = next;
    next = buffer[q];
  }

  /* release the job types that would have been
   * assigned to computer i upon arrival */
  nb_released = 0;
  t = 0;
  while ( (k = Ki[i][t]) != -1 ) {
    if (assignment[k] == i) {
      assignment[k] = -1;
      ++nb_released;
    }
    ++t;
  }

  /* reassign the released job types */
  if (head != tail) {
    p = (p + 1) % MAX_L;
    while ( nb_released && (p != tail) ) {
      current = buffer[p];
      t = 0;
      while ( (k = Ki[current][t]) != -1 ) {
        if (assignment[k] == -1) {
          assignment[k] = current;
          --nb_released;
        }
        ++t;
      }
      p = (p + 1) % MAX_L;
    }
  }
}





/*** SIMULATIONS ***/

static inline int draw_arrival_type () {
  int k = 0;
  double u = uniform();
  double v = nuk[0];

  while (u > v) {
    ++k;
    v += nuk[k];
  }

  return k;
}

static inline int draw_departure_computer () {
  int i = 0;
  double u = uniform();
  double v = mui[0];

  while (u > v) {
    ++i;
    v += mui[i];
  }

  return i;
}

static inline int draw_departure_type (int i) {
  int k = 0;
  int u = random() % xi[i];
  int v = xik[i][0];

  while (u >= v) {
    ++k;
    v += xik[i][k];
  }

  return k;
}

static inline int jump () {
  int i, k;

  if (uniform() < threshold) {
    /* choose the incoming job type */
    k = draw_arrival_type();

    /* update the system state */
    if ( (i = assignment[k]) != -1 ) {
      remove_from_queue(i);
      ++xi[i]; ++xik[i][k]; ++xk[k];
      return 1;
    }
  } else {
    /* choose the computer */
    i = draw_departure_computer();

    /* update the system state */
    if (xi[i] > 0) {
      k = draw_departure_type(i);
      --xi[i]; --xik[i][k]; --xk[k];
      append_to_queue(i);
      return 1;
    }
  }

  return 0;
}

static inline int jump_and_update_counters () {
  int i, k;

  if (uniform() < threshold) {
    /* choose the incoming job type */
    k = draw_arrival_type();

    /* update the system state and the counters */
    ++arrival[k];
    if ( (i = assignment[k]) != -1 ) {
      remove_from_queue(i);
      ++xi[i]; ++xik[i][k]; ++xk[k];
      return 1;
    } else ++loss[k];
  } else {
    /* choose the computer */
    i = draw_departure_computer();

    /* update the system state and the counters */
    ++timeout[i];
    if (xi[i] > 0) {
      k = draw_departure_type(i);
      --xi[i]; --xik[i][k]; --xk[k];
      append_to_queue(i);
      return 1;
    } else ++empty[i];
  }

  return 0;
}

static inline void update_metrics () {
  int i, k;
  double t;
 
  T += ( t = exponential(event_rate) );
  for (k = 0 ; k < K ; ++k) Lk[k] += xk[k] * t;
  for (i = 0 ; i < I ; ++i) Li[i] += xi[i] * t;
}

void simulation () {
  int i, k, step;

  /* initialization */
  initialize_queue();
  T = 0.;
  for (i = 0 ; i < I ; ++i) {
    xi[i] = 0;
    for (k = 0 ; k < K ; ++k) xik[i][k] = 0;
    Li[i] = 0.;
    timeout[i] = 0;
    empty[i] = 0;
  }

  for (k = 0 ; k < K ; ++k) {
    xk[k] = 0;
    Lk[k] = 0.;
    arrival[k] = 0;
    loss[k] = 0;
  }

  /* warmup */
  step = 0;
  while (step < warmup) step += jump();

  /* steady state */
  step = 0;
  while (step < steady) {
    step += jump_and_update_counters();
    update_metrics();
  }

  /* compute metrics */
  for (k = 0 ; k < K ; ++k) {
    betak[k] = 1. * loss[k] / arrival[k];
    Lk[k] /= T;
  }
  for (i = 0 ; i < I ; ++i) {
    Li[i] /= T;
    etai[i] = 1. - 1. * empty[i] / timeout[i];
  }
}





/*** INPUTS AND OUTPUTS ***/

void read_inputs (int argc, char **argv) {
  int c, i, k, t, L;
  char *arg;
  double norm;

  /* default initialization */
  sprintf(name, "toto");
  K = 0; I = 0;
  for (i = 0 ; i < MAX_I ; ++i) {
    elli[i] = 0;
    Ki[i][0] = -1;
    mui[i] = 0.;
  }
  for (k = 0 ; k < MAX_K ; ++k) nuk[k] = 0.;
  Rho = 5.;
  delta = .2;
  R = 100;
  warmup = 1000000;
  steady = 1000000;

  /* read arguments */
  opterr = 0;
  while ( (c = getopt(argc, argv, "f:k:i:l:c:r:m:d:R:w:t:h")) != -1 ) {
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
        L = 0;
        arg = strtok(optarg, " ");
        i = 0;
        while (i < I && arg != NULL) {
          elli[i] = atoi(arg);
          L += elli[i];
          ++i;
          arg = strtok(NULL, " ");
        }

        if (i == 0) {
          fprintf(stderr, "Option -l: you should specify the number of tokens.\n");
          exit(EXIT_FAILURE);
        }

        while (i < I) {
          elli[i] = elli[i-1];
          L += elli[i];
          ++i;
        }

        if (L >= MAX_L) {
          fprintf(stderr, "Option -l: the number of tokens is too large"
              " for the buffer length.\n");
          exit(EXIT_FAILURE);
        }
        break;
      case 'c':
        /* computer-to-type compatibilities,
         * encoded as an adjacency matrix between computers and types */
        arg = strtok(optarg, " ");

        for (i = 0 ; i < I ; ++i) {
          t = 0;

          for (k = 0 ; k < K ; ++k) {
            if (arg == NULL) {
              fprintf(stderr, "Option -c: you didn't specify all the compatibilities.\n");
              exit(EXIT_FAILURE);
            }

            if (atoi(arg)) {
              Ki[i][t] = k;
              ++t;
            }

            arg = strtok(NULL, " ");
          }

          Ki[i][t] = -1;
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
      case 'R':
        R = atoi(optarg);

        if (R > MAX_R) {
          fprintf(stderr, "Option -R: you cannot run more than %d independent runs.\n", MAX_R);
          exit(EXIT_FAILURE);
        } else if (R <= 1) {
          fprintf(stderr, "Option -R: you cannot run less than 2 independent runs.\n");
          exit(EXIT_FAILURE);
        }

        break;
      case 'w':
        /* number of warmup steps */
        warmup = pow(10, atoi(optarg));
        break;
      case 't':
        /* number of steps in steady state */
        steady = pow(10, atoi(optarg));
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
            "    (default value: 0.2)\n"
            "-R: number of simulation runs per load\n"
            "    (default value: 10; shouldn't be less than 2)\n"
            "-w: number of warmup steps = 10^argument\n"
            "    (default value: 10^6)\n"
            "-t: number of steps in steady state = 10^argument\n"
            "    (default value: 10^6)\n"
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
  int i, k, t;

  printf("\n--- DYNAMIC - SIMULATION - EXPONENTIAL ---\n\n");

  /* output file prefix name */
  printf("Output file prefix name:         %s.csv\n\n", name);

  /* simulation parameters */
  printf("Maximum load to be considered:      %.3e\n", Rho);
  printf("Step by which we increase the load: %.3e\n", delta);
  printf("Number of simulation runs per load: %d\n", R);
  printf("Number of warm-up steps:            %d\n", warmup);
  printf("Number of steps in steady-state:    %d\n\n", steady);

  /* numbers of types and computers */
  printf("%d job type(s) and %d computer(s)\n\n", K, I);

  /* numbers of tokens of each computer */
  for (i = 0 ; i < I ; ++i) printf("Computer %d: %d token(s)\n", i, elli[i]);
  printf("\n");

  /* computer-to-type compatibilities */
  printf("Compatibilities of each computer:\n");
  for (i = 0 ; i < I ; ++i) {
    printf("Computer %d <- ", i);
    t = 0;
    while ( (k = Ki[i][t]) != -1) {
      printf("%d ", k);
      ++t;
    }
    printf("\n");
  }
  printf("\n");

  /* arrival rates */
  printf("Distribution of the arrivals:\n");
  for (k = 0 ; k < K ; ++k) {
    printf("%.3e  ", nuk[k]);
    if (!(k % 6)) printf("\n");
  }
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





/*** MAIN  ***/

int main (int argc, char **argv) {
  int i, k, r;
  FILE *file;
  double rho;

  double mean_betak[MAX_K], dev_betak[MAX_K],
         mean_Lk[MAX_K], dev_Lk[MAX_K],
         mean_gammak[MAX_K], dev_gammak[MAX_K];

  double mean_etai[MAX_I], dev_etai[MAX_I],
         mean_Li[MAX_I], dev_Li[MAX_I],
         mean_gammai[MAX_I], dev_gammai[MAX_I];

  double mean_beta, dev_beta,
         mean_eta, dev_eta,
         mean_L, dev_L,
         mean_gamma, dev_gamma;

  double betark[MAX_R][MAX_K],
         Lrk[MAX_R][MAX_K],
         gammark[MAX_R][MAX_K];

  double etari[MAX_R][MAX_I],
         Lri[MAX_R][MAX_I],
         gammari[MAX_R][MAX_I];

  double betar[MAX_R], etar[MAX_R], Lr[MAX_R], gammar[MAX_R];

  /* input arguments */
  read_inputs(argc, argv);
  print_inputs();

  /* open the output file */
  file = open_csv_file();
  fprintf(file, "rho");
  for (k = 0 ; k < K ; ++k)
    fprintf(file, ",betak%d,wbetak%d,Lk%d,wLk%d,gammak%d,wgammak%d",
        k+1, k+1, k+1, k+1, k+1, k+1);
  for (i = 0 ; i < I ; ++i)
    fprintf(file, ",etai%d,wetai%d,Li%d,wLi%d,gammai%d,wgammai%d",
        i+1, i+1, i+1, i+1, i+1, i+1);
  fprintf(file, ",beta,wbeta,eta,weta,L,wL,gamma,wgamma\n");

  for (rho = delta ; rho < Rho + .5 * delta ; rho += delta) {
    printf("rho = %.3e\n", rho);
    event_rate = rho + 1.;
    threshold = rho / event_rate;

    /* initialiaze the means and the deviations */
    for (k = 0 ; k < K ; ++k) {
      mean_betak[k] = 0.; dev_betak[k] = 0.;
      mean_Lk[k] = 0.; dev_Lk[k] = 0.;
      mean_gammak[k] = 0.; dev_gammak[k] = 0.;
    }

    for (i = 0 ; i < I ; ++i) {
      mean_etai[i] = 0.; dev_etai[i] = 0.;
      mean_Li[i] = 0.; dev_Li[i] = 0.;
      mean_gammai[i] = 0.; dev_gammai[i] = 0.;
    }

    mean_beta = 0.; dev_beta = 0.;
    mean_eta = 0.; dev_eta = 0.;
    mean_L = 0.;  dev_L = 0.;
    mean_gamma = 0.; dev_gamma = 0.;

    for (r = 0 ; r < R ; ++r) {
      /* run the simulation */
      srandom(time(NULL));
      simulation();

      /* save the results */
      betar[r] = 0.;
      Lr[r] = 0.;
      etar[r] = 0.;

      for (k = 0 ; k < K ; ++k) {
        betark[r][k] = betak[k];
        betar[r] += nuk[k] * betak[k];
        Lrk[r][k] = Lk[k];
        Lr[r] += Lk[k];
        gammark[r][k] = rho * mu * nuk[k] * (1. - betak[k]) / Lk[k];
      }

      for (i = 0 ; i < I ; ++i) {
        etari[r][i] = etai[i];
        etar[r] += mui[i] * etai[i];
        Lri[r][i] = Li[i];
        gammari[r][i] = mu * mui[i] * etai[i] / Li[i];
      }

      gammar[r] = rho * mu * (1. - betar[r]) / Lr[r];

      /* update the means */
      for (k = 0 ; k < K ; ++k) {
        mean_betak[k] += betark[r][k];
        mean_Lk[k] += Lrk[r][k];
        mean_gammak[k] += gammark[r][k];
      }

      for (i = 0 ; i < I ; ++i) {
        mean_etai[i] += etari[r][i];
        mean_Li[i] += Lri[r][i];
        mean_gammai[i] += gammari[r][i];
      }

      mean_beta += betar[r];
      mean_eta += etar[r];
      mean_L += Lr[r];
      mean_gamma += gammar[r];
    }

    /* normalize the means */

    for (k = 0 ; k < K ; ++k) {
      mean_betak[k] /= R;
      mean_Lk[k] /= R;
      mean_gammak[k] /= R;
    }

    for (i = 0 ; i < I ; ++i) {
      mean_etai[i] /= R;
      mean_Li[i] /= R;
      mean_gammai[i] /= R;
    }

    mean_beta /= R;
    mean_eta /= R;
    mean_L /= R;
    mean_gamma /= R;

    /* compute the interval widths */

    for (r = 0 ; r < R ; ++r) {
      for (k = 0 ; k < K ; ++k) {
        dev_betak[k] += pow(betark[r][k] - mean_betak[k], 2);
        dev_Lk[k] += pow(Lrk[r][k] - mean_Lk[k], 2);
        dev_gammak[k] += pow(gammark[r][k] - mean_gammak[k], 2);
      }

      for (i = 0 ; i < I ; ++i) {
        dev_etai[i] += pow(etari[r][i] - mean_etai[i], 2);
        dev_Li[i] += pow(Lri[r][i] - mean_Li[i], 2);
        dev_gammai[i] += pow(gammari[r][i] - mean_gammai[i], 2);
      }

      dev_beta += pow(betar[r] - mean_beta, 2);
      dev_eta += pow(etar[r] - mean_eta, 2);
      dev_L += pow(Lr[r] - mean_L, 2);
      dev_gamma += pow(gammar[r] - mean_gamma, 2);
    }

    /* print data in the output file */
    fprintf(file, "%e", rho);

    for (k = 0 ; k < K ; ++k) {
      fprintf(file, ",%e,%e", mean_betak[k], alpha * sqrt(dev_betak[k] / (R-1)) / R);
      fprintf(file, ",%e,%e", mean_Lk[k], alpha * sqrt(dev_Lk[k] / (R-1)) / R);
      fprintf(file, ",%e,%e", mean_gammak[k], alpha * sqrt(dev_gammak[k] / (R-1)) / R);
    }

    for (i = 0 ; i < I ; ++i) {
      fprintf(file, ",%e,%e", mean_etai[i], alpha * sqrt(dev_etai[i] / (R-1)) / R);
      fprintf(file, ",%e,%e", mean_Li[i], alpha * sqrt(dev_Li[i] / (R-1)) / R);
      fprintf(file, ",%e,%e", mean_gammai[i], alpha * sqrt(dev_gammai[i] / (R-1)) / R);
    }

    fprintf(file, ",%e,%e", mean_beta, alpha * sqrt(dev_beta / (R-1)) / R);
    fprintf(file, ",%e,%e", mean_eta, alpha * sqrt(dev_eta / (R-1)) / R);
    fprintf(file, ",%e,%e", mean_L, alpha * sqrt(dev_L / (R-1)) / R);
    fprintf(file, ",%e,%e\n", mean_gamma, alpha * sqrt(dev_gamma / (R-1)) / R);
  }

  /* close the output file */
  fclose(file);
}
