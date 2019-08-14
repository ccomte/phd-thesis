#include "utils.h"

#define N 2           // number of job classes
#define S 3           // number of servers
#define MAX 100000    // queue length
#define R 5           // number of estimates per point





/*** INPUT ***/

char name[500];     // prefix of the output file

long double proportion[N];   // distribution of the arrivals between classes

int degree[N];      // number of servers that can process each class
int* neighbors[N];  // array of servers which can serve the jobs of each class

// Hyperexponential job sizes
// A proportion p of jobs have an exponentially distributed size with mean sigma1
// and a proportion 1-p of jobs have an exponentially distributed size with mean sigma2
long double p, sigma1, sigma2;





/*** VARIABLES OF THE SIMULATION ***/

int warmup;         // number of warm-up steps
int steady;         // number of steps in steady-state
long double delta;  // scaling of the interval

int x[N];         // number of jobs of each class in the queue

long double lambda[N];        // per-class arrival rates
long double total_lambda;     // total arrival rate
long double total_completion; // total completion rates

int is_served[N]; // says if a job of the class was already served
int is_active[S]; // tracks the activity of the servers





/*** OUTPUTS ***/

long double T;                  // total simulation time
long double xT[N];              // cumulative number of jobs of each class
long double delay[N][R];        // delays experienced at each run
long double average_delay[N];   // average delay experienced over all runs
long double service_rate[N][R];       // service rates experienced at each run
long double average_service_rate[N];  // average service rate experienced over all runs





/*** PROTOTYPES ***/

// Queue

void arrival (int i, long double sigma);
// add a job of class i with mean size sigma at the end of the queue
void departure (int k);     // remove the job in position k of the queue
void update_rates ();       // update the service rates of the jobs

// Simulation

void simulation ();
void jump ();

// Inputs and outputs

void read_inputs_from_string (char **argv);
void print_inputs ();
FILE *open_input_file (char *type);
void initialize_lambda (long double rho);





/*** QUEUE ***/

int n;            // number of jobs in the queue
int class[MAX];   // sequence of job classes
long double size[MAX];        // sequence of job sizes
long double completion[MAX];  // sequence of jobs phase completion rates


void arrival (int i, long double sigma) {
  // update the macrostate
  ++x[i];

  // update the microstate
  class[n] = i;
  size[n] = sigma;
  ++n;

  // check the length of the queue
  if (n == MAX) {
    fprintf(stderr, "THE QUEUE IS TOO SHORT.\n");
    exit(EXIT_FAILURE);
  }
}

void departure (int k) {
  int l;

  // update the macrostate
  --x[class[k]];

  // update the microstate
  --n;
  for (l = k ; l < n ; ++l) {
    class[l] = class[l+1];
    size[l] = size[l+1];
  }
}

void update_rates () {
  int i, j, k, s, total_active_servers;

  // initialization
  total_completion = 0.;
  for (i = 0 ; i < N ; ++i) is_served[i] = 0;
  for (s = 0 ; s < S ; ++s) is_active[s] = 0;
  total_active_servers = 0;

  // greedy server allocation
  k = 0;
  while (k < n && total_active_servers < S) {
    completion[k] = 0.;
    i = class[k];

    if (!is_served[i]) {
      is_served[i] = 1;

      // look for free servers
      for (j = 0 ; j < degree[i] ; ++j) {
        s = neighbors[i][j];

        if (!is_active[s]) {
          is_active[s] = 1;
          ++total_active_servers;
          completion[k] += 1. / size[k];
        }
      }

      total_completion += completion[k];
    }

    ++k;
  }
}





/*** SIMULATIONS ***/

void simulation () {
  int i, step;
  long double t, T;

  // initialization
  T = 0.;
  n = 0;
  total_completion = 0.;
  for (i = 0 ; i < N ; ++i) {
    x[i] = 0;
    xT[i] = 0.;
  }

  // warmup
  for (step = 0 ; step < warmup ; ++step) {
    jump();
    update_rates();
  }

  // steady state
  for (step = 0 ; step < steady ; ++step) {
    jump();
    update_rates();

    // metric update
    t = exponential(total_lambda + total_completion);
    T += t;
    for (i = 0 ; i < N ; ++i) xT[i] += x[i] * t;
  }

  // metrics
  for (i = 0 ; i < N ; ++i) xT[i] /= T;
}

void jump () {
  int i, j, k;
  long double u, v;

  u = (total_lambda + total_completion) * uniform();

  if (u <= total_lambda) {
    // external arrival

    v = lambda[0];
    i = 0;
    while (u > v) {
      ++i;
      v += lambda[i];
    }

    // size
    if (uniform() < p) arrival(i, sigma1);
    else arrival(i, sigma2);

  } else {
    // completion of a job

    v = total_lambda + completion[0];
    k = 0;
    while (u > v) {
      ++k;
      v += completion[k];
    }

    departure(k);
  }
}





/*** INPUTS AND OUTPUTS ***/

void read_inputs_from_string (char **argv) {
  int i, j, m, r, s;
  long double total;

  r = 0;

  // output file prefix name
  sprintf(name, "%s", argv[++r]);

  // arrival rates
  total = 0.;
  for (i = 0 ; i < N ; ++i) {
    proportion[i] = atof(argv[++r]);
    total += proportion[i];
  }

  for (i = 0 ; i < N ; ++i) proportion[i] /= total;

  // compatibility constraints
  for (i = 0 ; i < N ; ++i) {
    degree[i] = atoi(argv[++r]);
    neighbors[i] = malloc(degree[i] * sizeof(int));
    for (j = 0 ; j < degree[i] ; ++j) {
      neighbors[i][j] = atoi(argv[++r]);
    }
  }

  // parameter of the hyperexponential distribution
  p = atof(argv[++r]);
  p /= p + atof(argv[++r]);
  sigma1 = atof(argv[++r]);
  sigma2 = atof(argv[++r]);
}

void print_inputs () {
  int i, j, s;

  printf("\n\n|| ----- SIMULATION - GRAPH - HYPEREXPONENTIAL ----- ||\n\n");

  printf("Output file prefix name: %s\n\n", name);

  printf("N = %d\tS = %d\n\n", N, S);

  printf("Normalized arrival rates:\n");
  for (i = 0 ; i < N ; ++i) printf("Class %d:\t%Le\n", i, proportion[i]);
  printf("\n");

  printf("Compatibility constraints:\n");
  for (i = 0 ; i < N ; ++i) {
    printf("Class %d:\t", i);
    for (j = 0 ; j < degree[i] ; ++j)
      printf("%d\t", neighbors[i][j]);
    printf("\n");
  }
  printf("\n");

  printf("Parameter of the hyperexponential distribution:\n");
  printf("p = %Le\tsigma1 = %Le\tsigma2 = %Le\n", p, sigma1, sigma2);
  printf("Mean job size = %Le\n\n", p * sigma1 + (1.-p) * sigma2);
}

FILE *open_output_file (char *type) {
  int i;
  char file_name[500];
  FILE *file = NULL;

  // build the file name
  sprintf(file_name, "%s-%s", name, type);
  file = fopen(file_name, "w");

  // exit with an error if the file couldn't be opened
  if (file == NULL) {
    fprintf(stderr, "Couldn't open the output file %s.\n\n\n", file_name);
    exit(EXIT_FAILURE);
  } else printf("Save %s in the file %s.\n", type, file_name);

  // print the column headers
  fprintf(file, "rho");
  for (i = 0 ; i < N ; ++i) {
    fprintf(file, ",performance%d,interval%d", i, i);
  }
  fprintf(file, "\n");

  return file;
}

void initialize_lambda (long double rho) {
  int i;

  total_lambda = 0.;

  for (i = 0 ; i < N ; ++i) {
    lambda[i] = rho * S * proportion[i] / (p * sigma1 + (1-p) * sigma2);
    total_lambda += lambda[i];
  }
}





int main(int argc, char **argv) {
  int i, r;
  long double rho, deviation;
  FILE *delay_file = NULL, *service_rate_file = NULL;

  srand(time(NULL));

  read_inputs_from_string(argv);
  print_inputs();

  delay_file = open_output_file("delay");
  service_rate_file = open_output_file("service-rate");
  printf("\n");

  warmup = 10000;
  steady = 10000;
  delta = 1.96;     // P(-delta < X < delta) = 95% (with X ~ N(0,1))

  for (rho = .05 ; rho < .995 ; rho += .05) {
    printf("rho = %Le\n", rho);

    // initialization
    initialize_lambda(rho);
    for (i = 0 ; i < N ; ++i) {
      average_delay[i] = 0.;
      average_service_rate[i] = 0.;
    }

    // R runs of simulation
    for (r = 0 ; r < R ; ++r) {
      simulation();

      for (i = 0 ; i < N ; ++i) {
        delay[i][r] = xT[i] / lambda[i];
        average_delay[i] += delay[i][r];

        service_rate[i][r] = lambda[i] * (p * sigma1 + (1-p) * sigma2) / xT[i];
        average_service_rate[i] += service_rate[i][r];
      }
    }

    // output delay
    fprintf(delay_file, "%Le", rho);
    for (i = 0 ; i < N ; ++i) {
      average_delay[i] /= R;
      deviation = 0.;
      for (r = 0 ; r < R ; ++r)
        deviation += pow(delay[i][r] - average_delay[i], 2);
      fprintf(delay_file, ",%Le,%Le",
          average_delay[i], delta * sqrt(deviation / (R-1)) / sqrt(R));
    }
    fprintf(delay_file, "\n");

    // output service rate
    fprintf(service_rate_file, "%Le", rho);
    for (i = 0 ; i < N ; ++i) {
      average_service_rate[i] /= R;
      deviation = 0.;
      for (r = 0 ; r < R ; ++r)
        deviation += pow(service_rate[i][r] - average_service_rate[i], 2);
      fprintf(service_rate_file, ",%Le,%Le",
          average_service_rate[i], delta * sqrt(deviation / (R-1)) / sqrt(R));
    }
    fprintf(service_rate_file, "\n");
  }

  for (i = 0 ; i < N ; ++i) free(neighbors[i]);
  fclose(delay_file);
  fclose(service_rate_file);

  printf("\n\n");
  exit(0);
}
