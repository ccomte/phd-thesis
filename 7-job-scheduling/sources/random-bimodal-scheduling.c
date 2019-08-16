#include "utils.h"
// defines the constants S, MAX, R, WARMUP and STEADY





/*** INPUT ***/

char name[500];     // prefix of the output file

int d;              // degree of parallelism
long double tau;    // timeout rate of the servers

// job sizes with a bimodal number of exponentially distributed phases
// a proportion p of jobs have sigma1 phases
// and a proportion 1-p of jobs have sigma2 phases
long double p;
int sigma1, sigma2;





/*** VARIABLES OF THE SIMULATION ***/

long double delta;  // scaling of the confidence interval

long double lambda;           // external arrival rate
long double total_timeout;    // total timeout rates
long double total_completion; // total completion rates

int is_active[S]; // tracks the activity of the servers





/*** OUTPUTS ***/

long double T;                // total simulation time
long double xT;               // cumulative number of jobs
long double delay[R];         // delays experienced at each run
long double average_delay;    // average delay experienced over all runs
long double service_rate[R];        // service rates experienced at each run
long double average_service_rate;   // average service rate experienced over all runs





/*** PROTOTYPES ***/

// Queue

void arrival (int sigma);
// add a job with size sigma at the end of the queue and draws its servers
void departure (int k);     // remove the job in position k of the queue
void interruption (int k);  // move the job in position k to the end of the queue
void update_rates ();       // update the service rates of the jobs

// Simulation

void simulation ();
void jump ();

// Inputs and outputs

void read_inputs_from_string (char **argv);
void print_inputs ();
FILE *open_input_file (char *type);





/*** QUEUE ***/

int n;                // number of jobs in the queue
int* class[MAX];      // sequence of job server sets
int size[MAX];        // sequence of job sizes (number of remaining phases)
long double timeout[MAX];     // sequence of job timeout rates
long double completion[MAX];  // sequence of jobs phase completion rates


void arrival (int sigma) {
  int j, s, count, *servers;

  // draw the servers of the new job
  servers = malloc(d * sizeof(int));
  count = 0;
  while (count < d) {
    s = floor(S * uniform());
    j = 0;
    while (j < count && servers[j] != s) ++j;
    if (j == count) {
      servers[count] = s;
      ++count;
    }
  }

  // update the microstate
  class[n] = servers;
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

  // free the server array
  free(class[k]);

  // update the microstate
  --n;
  for (l = k ; l < n ; ++l) {
    class[l] = class[l+1];
    size[l] = size[l+1];
  }
}

void interruption (int k) {
  int l, sigma, *servers;

  servers = class[k];
  sigma = size[k];

  // update the microstate
  for (l = k ; l < n-1 ; ++l) {
    class[l] = class[l+1];
    size[l] = size[l+1];
  }

  class[n-1] = servers;
  size[n-1] = sigma;
}

void update_rates () {
  int j, k, s, total_active_servers;
  int *servers;

  // initialization
  total_timeout = 0.;
  total_completion = 0.;
  for (s = 0 ; s < S ; ++s) is_active[s] = 0;
  total_active_servers = 0;

  // greedy server allocation
  k = 0;
  while (k < n && total_active_servers < S) {
    timeout[k] = 0.;
    completion[k] = 0.;
    servers = class[k];

    // look for free servers
    for (j = 0 ; j < d ; ++j) {
      s = servers[j];

      if (!is_active[s]) {
        is_active[s] = 1;
        ++total_active_servers;
        timeout[k] += tau;
        completion[k] += 1.;
      }
    }

    total_timeout += timeout[k];
    total_completion += completion[k];

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
  total_timeout = 0.;
  total_completion = 0.;
  xT = 0.;

  // warmup
  for (step = 0 ; step < WARMUP ; ++step) {
    jump();
    update_rates();
  }

  // steady state
  for (step = 0 ; step < STEADY ; ++step) {
    jump();
    update_rates();

    // metric update
    t = exponential(lambda + total_timeout + total_completion);
    T += t;
    xT += n * t;
  }

  // metrics
  xT /= T;
}

void jump () {
  int k;
  long double u, v;

  u = (lambda + total_timeout + total_completion) * uniform();

  if (u <= lambda) {
    // external arrival

    if (uniform() < p) arrival(sigma1);
    else arrival(sigma2);

  } else if (u <= lambda + total_timeout) {
    // the timeout of a server expired

    v = lambda + timeout[0];
    k = 0;
    while (u > v) {
      ++k;
      v += timeout[k];
    }

    // move the job to the end of the queue
    interruption(k);

  } else {
    // completion of the phase of a job

    v = lambda + total_timeout + completion[0];
    k = 0;
    while (u > v) {
      ++k;
      v += completion[k];
    }

    // update the size
    --size[k];
    if (size[k] == 0) departure(k);
  }
}





/*** INPUTS AND OUTPUTS ***/

void read_inputs_from_string (char **argv) {
  int i, j, m, r, s;
  long double total;

  r = 0;

  // output file prefix name
  sprintf(name, "%s", argv[++r]);

  // degree of parallelism
  d = atoi(argv[++r]);

  // parameter of the "bimodal" distribution
  p = atof(argv[++r]);
  p /= p + atof(argv[++r]);
  sigma1 = atoi(argv[++r]);
  sigma2 = atoi(argv[++r]);

  // timer rate
  m = atoi(argv[++r]);
  tau = m / (p * sigma1 + (1.-p) * sigma2);
}

void print_inputs () {
  int i, j, s;

  printf("\n\n|| ----- SIMULATION - RANDOM - BIMODAL ----- ||\n\n");

  printf("Output file prefix name: %s\n\n", name);

  printf("S = %d\t\td = %d\n\n", S, d);

  printf("Parameter of the \"bimodal\" distribution:\n");
  printf("p = %Le\tsigma1 = %d\tsigma2 = %d\n\n", p, sigma1, sigma2);

  printf("Mean number of interruptions per job = %Le\n",
      tau * (p * sigma1 + (1.-p) * sigma2));
  printf("Timer rate = %Le\n\n", tau);
}

FILE *open_output_file (char *type) {
  char file_name[500];
  FILE *file = NULL;

  // build the file name
  sprintf(file_name, "%s-%s.csv", name, type);
  file = fopen(file_name, "w");

  // exit with an error if the file couldn't be opened
  if (file == NULL) {
    fprintf(stderr, "Couldn't open the output file %s.\n\n\n", file_name);
    exit(EXIT_FAILURE);
  } else printf("Save %s in the file %s.\n", type, file_name);

  // print the column headers
  fprintf(file, "rho,performance,interval\n");

  return file;
}





int main(int argc, char **argv) {
  int r;
  long double rho, deviation;
  FILE *delay_file = NULL, *service_rate_file = NULL;

  srand(time(NULL));

  read_inputs_from_string(argv);
  print_inputs();

  delay_file = open_output_file("delay");
  service_rate_file = open_output_file("service-rate");
  printf("\n");

  delta = 1.96;     // P(-delta < X < delta) = 95% (with X ~ N(0,1))

  for (rho = .05 ; rho < .995 ; rho += .05) {
    // initialization
    lambda = rho * S / (p * sigma1 + (1.-p) * sigma2);
    average_delay = 0.;
    average_service_rate = 0.;

    // R runs of simulation
    for (r = 0 ; r < R ; ++r) {
      simulation();

      delay[r] = xT / lambda;
      average_delay += delay[r];

      service_rate[r] = lambda * (p * sigma1 + (1-p) * sigma2) / xT;
      average_service_rate += service_rate[r];
    }

    // output delay
    average_delay /= R;
    deviation = 0.;
    for (r = 0 ; r < R ; ++r) {
      deviation += pow(delay[r] - average_delay, 2);
    }
    fprintf(delay_file, "%Le,%Le,%Le\n", rho,
        average_delay, delta * sqrt(deviation / (R-1)) / sqrt(R));

    // output service rate
    average_service_rate /= R;
    deviation = 0.;
    for (r = 0 ; r < R ; ++r) {
      deviation += pow(service_rate[r] - average_service_rate, 2);
    }
    fprintf(service_rate_file, "%Le,%Le,%Le\n", rho,
        average_service_rate, delta * sqrt(deviation / (R-1)) / sqrt(R));
  }

  fclose(delay_file);
  fclose(service_rate_file);

  printf("\n");
  exit(0);
}
