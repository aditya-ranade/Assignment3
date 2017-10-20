/**
 * Parallel VLSI Wire Routing via OpenMP
 * Name 1(andrew_id 1), Name 2(andrew_id 2)
 */

#include "wireroute.h"
#include <time.h>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <assert.h>
#include <omp.h>
#include "mic.h"
#include <limits.h>
#include <random>
#define BUFSIZE 1024

#define MIN(X, Y) (((X) < (Y)) ? (X) : (Y))
#define MAX(X, Y) (((X) < (Y)) ? (X) : (Y))

static int _argc;
static const char **_argv;

/* Starter code function, don't touch */
const char *get_option_string(const char *option_name,
			      const char *default_value)
{
  for (int i = _argc - 2; i >= 0; i -= 2)
    if (strcmp(_argv[i], option_name) == 0)
      return _argv[i + 1];
  return default_value;
}

/* Starter code function, do not touch */
int get_option_int(const char *option_name, int default_value)
{
  for (int i = _argc - 2; i >= 0; i -= 2)
    if (strcmp(_argv[i], option_name) == 0)
      return atoi(_argv[i + 1]);
  return default_value;
}

/* Starter code function, do not touch */
float get_option_float(const char *option_name, float default_value)
{
  for (int i = _argc - 2; i >= 0; i -= 2)
    if (strcmp(_argv[i], option_name) == 0)
      return (float)atof(_argv[i + 1]);
  return default_value;
}

/* Starter code function, do not touch */
static void show_help(const char *program_path)
{
    printf("Usage: %s OPTIONS\n", program_path);
    printf("\n");
    printf("OPTIONS:\n");
    printf("\t-f <input_filename> (required)\n");
    printf("\t-n <num_of_threads> (required)\n");
    printf("\t-p <SA_prob>\n");
    printf("\t-i <SA_iters>\n");
}





int increment_cost_x(cost_t *costs, cost_t *costs_t, int x1, int x2, int y, int dim_x, int dim_y, bool write, int sum) {
	int max_cost = 0;
	int multiplier = 1;
	int val = 0;
	if (x1 > x2) multiplier = -1;
	while (x1 != x2) {	
		val = costs[y*dim_x + x1] + sum;
		if (write) {
			costs[y*dim_x + x1] += sum;
			costs_t[x1*dim_y + y] += sum;
		}
		if (val > max_cost) max_cost = val;
		x1 += multiplier;
	}
	return max_cost;
}

int increment_cost_y(cost_t *costs, cost_t *costs_t, int y1, int y2, int x, int dim_x, int dim_y, bool write, int sum) {
	int max_cost = 0;
	int multiplier = 1;
	int val = 0;
	if (y1 > y2) multiplier = -1;
	while (y1 != y2) {
		val = costs_t[x*dim_y + y1] + sum; 
		if (write) {
			costs[y1*dim_x + x] += sum;
			costs_t[x*dim_y + y1] += sum;
		}
		if (val > max_cost) max_cost = val;
		y1 += multiplier;
	}
 return max_cost;
}

int max(int x, int y) {
	return (x < y) ? y : x;	
}


int traverse_path(int x1, int y1, int x2, int y2, int path, int bend, cost_t *costs, cost_t* costs_t, int dim_x,int dim_y, bool write, int sum) {
	int cost1; int cost2; int cost3; int cost4;
	if (path == 0) {
		cost1 = increment_cost_x(costs, costs_t, x1, bend, y1, dim_x, dim_y,  write, sum);
		cost2 = increment_cost_y(costs, costs_t, y1, y2, bend, dim_x, dim_y, write, sum);
		cost3 = increment_cost_x(costs, costs_t, bend, x2, y2, dim_x, dim_y, write, sum);
	}

	else {
		cost1 = increment_cost_y(costs, costs_t, y1, bend, x1, dim_x, dim_y, write, sum);
		cost2 = increment_cost_x(costs, costs_t, x1, x2, bend, dim_x, dim_y, write, sum);
		cost3 = increment_cost_y(costs, costs_t, bend, y2, x2, dim_x, dim_y, write, sum);
	}
	cost4 = costs[y2*dim_x + x2] + sum;	
	if (write) {
		costs[y2*dim_x + x2] = cost4;
		costs_t[x2*dim_y + y2] = cost4;
	}
	
	return max(max(max(cost1, cost2), cost3), cost4);
}


void traverse_path1(int x1, int y1, int x2, int y2, int path, int bend, cost_t *costs, cost_t *costs_t, int dimx,int dimy, int sum, int *min_costs, int* agg_costs, int i) {
	int val = 0;
	int x_multiplier = (x1 > x2) ? -1 : 1;
	int y_multiplier = (y1 > y2) ? -1 : 1;
	int maxCost = 0;
	int aggCost =0;
	int cost4;

	int val1, val2, val3, val4;
	int aggCost1, aggCost2, aggCost3, aggCost4;
	int maxCost1, maxCost2, maxCost3;
	

	if (path == 0) {	
		while (x1 != bend) {
			val1 = costs[y1*dimx + x1] + 1;
			aggCost1 += val1;
			if (val1 > maxCost1) maxCost1 = val1;		
			x1 += x_multiplier;	
		}
		while (y1 != y2) {
			val2 = costs_t[x1*dimy + y1] + 1;
			aggCost2 += val2;
			if (val2 > maxCost2) maxCost2 = val2;
			y1 += y_multiplier;		
		}
		while (x1 != x2) {
			val3 = costs[y2*dimx + x1] + 1;
			aggCost3 += val3;
			if (val3 > maxCost3) maxCost3 = val3;
			x1 += x_multiplier;
		}
	}
	else {
		while (y1 != bend) {
			val1 = costs_t[x1*dimy + y1] + 1;
			aggCost1 += val1;
			if (val1 > maxCost1) maxCost1 = val1;		
			y1 += y_multiplier;	
		}
		while (x1 != x2) {
			val2 = costs[y1*dimx + x1] + 1;
			aggCost2 += val2;
			if (val2 > maxCost2) maxCost2 = val2;
			x1 += x_multiplier;		
		}
		while (y1 != y2) {
			val3 = costs_t[x1*dimy + y1] + 1;
			aggCost3 += val3;
			if (val3 > maxCost3) maxCost3 = val3;
			y1 += y_multiplier;
		}
	}
	val4 = costs[y2*dimx + x2] + 1;
	aggCost4 = val4;
	min_costs[i] = max(max(max(maxCost1, maxCost2), maxCost3), val4);
	agg_costs[i] = aggCost1 + aggCost2 + aggCost3 + aggCost4;
}




void all_paths(int x1, int y1, int x2, int y2, int *res) {
	int i = 0;
	int multiplier = (x1 > x2) ? -1 : 1;
	
	int counter = x1 + multiplier;
	if (x1 != x2) {
		while (counter != x2) {
			res[i] = counter;
			counter += multiplier;
			i++;
		}
		res[i] = counter;
		i++;
	}
	if (y1 != y2) {
		multiplier = (y1 > y2) ? -1 : 1;
		counter = y1 + multiplier;
		while (counter != y2) {
			res[i] = counter;
			counter += multiplier;
			i++;
		}
		res[i] = counter;
		i++;
	}
}

void check_paths(int x1, int y1, int x2, int y2, cost_t *costs, cost_t *costs_t, int dimx, int dimy,  int *min_costs, int *agg_costs)  {
	
	int path, bend, i, curr_cost;
	int num_paths_x = abs(x1 - x2);
	int num_paths_y = abs(y1 - y2);
	int num_paths = num_paths_x + num_paths_y;
	int x_multiplier = (x1 > x2) ? -1 : 1;
	int y_multiplier = (y1 > y2) ? -1 : 1;
	
	
#pragma omp parallel for num_threads(64) schedule(dynamic) private(i, path, bend)
	for (i = 0; i < num_paths; i++) {
		if (i < num_paths_x) {
			path = 0;
			bend = x1 + (x_multiplier * (i + 1)); 
			traverse_path1(x1, y1, x2, y2, path, bend, costs, costs_t, dimx, dimy, 1, min_costs, agg_costs, i);
		}
		else {
			path = 1;
			bend = y1 + (y_multiplier * (i - num_paths_x + 1));
			traverse_path1(x1, y1, x2, y2, path, bend, costs, costs_t, dimx,dimy, 1, min_costs, agg_costs, i);
		}
	}
}
	
	
	
	
/*	
	 	
	int thread_index = omp_get_thread_num();
	int start = thread_index * chunk_size;
	int i = 0;
	int min_cost = INT_MAX;:w

	int best_path;
	int bend;
	int curr_cost = min_cost;
	while (i < chunk_size && (start + i) < w.num_paths) {
		if ((start + i) < num_paths_x) {
			w.bend = w.c1[0] + (x_multiplier * (start + i + 1)); 
			w.path = 0;
			curr_cost = traverse_path(w, costs, dimx, false, 1);
		}
		else {
		  int offset = start + i - num_paths_x;
			w.path = 1;
			w.bend = w.c1[1] + (y_multiplier * (offset + 1));
			curr_cost = traverse_path(w, costs, dimx, false, 1);
		}
		if (curr_cost < min_cost) {
			min_cost = curr_cost;
			best_path = w.path;
			bend = w.bend;
		}
		i++;
	}

	min_costs[thread_index] = min_cost;
	paths[thread_index] = best_path;
	bends[thread_index] = bend;
}
*/

void set_random_path(wire_t *w) {
		int x_multiplier = (w->c1[0] > w->c2[0]) ? -1 : 1;
		int y_multiplier = (w->c1[1] > w->c2[1]) ? -1 : 1;
		int num_paths_x = abs(w->c1[0] - w->c2[0]);
		int num_paths_y = abs(w->c1[1] - w->c2[1]);
	
		int number = rand() % (num_paths_x + num_paths_y);
		if (number < num_paths_x) {
			w->path = 0;
			w->bend = w->c1[0] + x_multiplier*(number + 1);
		}
		else {
			w->path = 1;
			int offset = number - num_paths_x;
			w->bend = w->c1[1] + y_multiplier*(offset + 1);
		}
}

void fill_costs(cost_t *costs, cost_t *costs_t, int dimx, int dimy, wire_t *wires, int num_of_wires, double prob) {
	int minCost = std::numeric_limits<int>::max();
	float P = (float)prob*100;
	srand (static_cast <unsigned> (time(0)));
	int *min_costs = (int *) malloc(sizeof(int)*dimx*dimy);
	int *agg_costs = (int *)malloc(sizeof(int)*dimx*dimy);
			
	
	for (int i = 0; i < num_of_wires; i++) {
		int best_path_counter;
		int best_path;
		wire_t *w = &wires[i];
		traverse_path(w->c1[0], w->c1[1], w->c2[0], w->c2[1], w->path, w->bend, costs, costs_t, dimx,dimy, true, -1);
		
		float choice = static_cast <float> (rand()) / (static_cast <float> (RAND_MAX/100));
	  int num_paths_x = abs(w->c1[0] - w->c2[0]);	
		int x_multiplier = (w->c1[0] > w->c2[0]) ? -1 : 1;
		int y_multiplier = (w->c1[1] > w->c2[1]) ? -1 : 1;
		if (choice > 10) {
			
			check_paths(w->c1[0], w->c1[1], w->c2[0], w->c2[1], costs, costs_t, dimx,dimy, min_costs, agg_costs);
			int index;
			int agg_cost;
			int try_bend, try_path, new_agg_cost;
			for (int i = 0; i < w->num_paths; i++) {
				if (min_costs[i] < minCost) {
					minCost = min_costs[i];
					index = i;
				}
				if (min_costs[i] == minCost) {
					if (agg_costs[i] < agg_costs[index]) index = i;
				}
			}
			
			if (index < num_paths_x) {
				 w->path = 0;
				 w->bend = w->c1[0] + (x_multiplier *(index+ 1));
			}
			else {
				w->path = 1;
				w->bend = w->c1[1] + (y_multiplier * (index - num_paths_x + 1));

			}
		}
		else {
			set_random_path(w);	
		}
		traverse_path(w->c1[0], w->c1[1], w->c2[0], w->c2[1], w->path, w->bend, costs,costs_t, dimx, dimy, true, 1);
	}

	free(min_costs);
	free(agg_costs);
	/*int agg = 0;
	int maxi = 0;
	int curr_max = 0;
	for (int i = 0; i < dimx*dimy; i++) {
		if (costs[i] > 1) agg += costs[i];
		if (costs[i] > maxi) maxi = costs[i];
	}
	fprintf(stdout, "aggregate = %d, max = %d \n", agg, maxi);
	*/
}

int main(int argc, const char *argv[])
{
  using namespace std::chrono;
  typedef std::chrono::high_resolution_clock Clock;
  typedef std::chrono::duration<double> dsec;

  auto init_start = Clock::now();
  double init_time = 0;

  _argc = argc - 1;
  _argv = argv + 1;

  /* You'll want to use these parameters in your algorithm */
  const char *input_filename = get_option_string("-f", NULL);
  int num_of_threads = get_option_int("-n", 1);
  double SA_prob = get_option_float("-p", 0.1f);
  int SA_iters = get_option_int("-i", 5);

  int error = 0;

  if (input_filename == NULL) {
    printf("Error: You need to specify -f.\n");
    error = 1;
  }

  if (error) {
    show_help(argv[0]);
    return 1;
  }

  printf("Number of threads: %d\n", num_of_threads);
  printf("Probability parameter for simulated annealing: %lf.\n", SA_prob);
  printf("Number of simulated anneling iterations: %d\n", SA_iters);
  printf("Input file: %s\n", input_filename);

  FILE *input = fopen(input_filename, "r");

  if (!input) {
    printf("Unable to open file: %s.\n", input_filename);
    return -1;
  }

  int dim_x, dim_y;
  int num_of_wires;

  fscanf(input, "%d %d\n", &dim_x, &dim_y);
  fscanf(input, "%d\n", &num_of_wires);

  wire_t *wires = (wire_t *)calloc(num_of_wires, sizeof(wire_t));
  (void)wires;

  int x1, y1, x2, y2;
  int i = 0;
  while (fscanf(input, "%d %d %d %d\n", &x1, &y1, &x2, &y2) != EOF) {
    /* PARSE THE INPUT FILE HERE.
     * Define wire_t in wireroute.h and store
     * x1, x2, y1, and y2 into the wires array allocated above
     * based on your wire_t definition. */

    wires[i].c1[0] = x1;
    wires[i].c2[0] = x2;
    wires[i].c1[1] = y1;
    wires[i].c2[1] = y2;
		wires[i].bend = 0;
    wires[i].path = 0;
		wires[i].num_paths = abs(x1 - x2) + abs(y1 - y2);
    set_random_path(&wires[i]);
		i++;
  }

  if (i != num_of_wires) {
    printf("Error: wire count mismatch");
    return -1;
  }
  cost_t *costs = (cost_t *)calloc(dim_x * dim_y, sizeof(cost_t));
  cost_t *costs_t = (cost_t *)calloc(dim_x * dim_y, sizeof(cost_t));
	(void)costs;
  /* INITIALIZE YOUR COST MATRIX HERE */

  /* Initialize additional data structures needed in the algorithm
   * here if you feel it's needed. */
	
	for (int i = 0; i < num_of_wires; i++) {
		traverse_path(wires[i].c1[0], wires[i].c1[1], wires[i].c2[0], wires[i].c2[1], wires[i].path, wires[i].bend, costs, costs_t, dim_x, dim_y, true, 1);	
  }
  error = 0;

  init_time += duration_cast<dsec>(Clock::now() - init_start).count();
  printf("Initialization Time: %lf.\n", init_time);

  auto compute_start = Clock::now();
  double compute_time = 0;


#ifdef RUN_MIC /* Use RUN_MIC to distinguish between the target of compilation */

  /* This pragma means we want the code in the following block be executed in
   * Xeon Phi.
   */
#pragma offload target(mic)  \
	inout(wires: length(num_of_wires) INOUT)    \
  inout(costs: length(dim_x*dim_y) INOUT)    \
	inout (costs_t: length(dim_x*dim_y) INOUT) 
#endif
  {
    /* Implement the wire routing algorithm here
     * Feel free to structure the algorithm into different functions
     * Don't use global variables.
     * Use OpenMP to parallelize the algorithm.
     * You should really implement as much of this (if not all of it) in
     * helper functions. */
	
	for (int i = 0; i < SA_iters; i++) {

		fill_costs(costs, costs_t, dim_x, dim_y, wires, num_of_wires, SA_prob);	
	}  
	}

  compute_time += duration_cast<dsec>(Clock::now() - compute_start).count();
  printf("Computation Time: %lf.\n", compute_time);

  /* OUTPUT YOUR RESULTS TO FILES HERE
   * When you're ready to output your data to files, uncommment this chunk of
   * code and fill in the specified blanks indicated by comments. More about
   * this in the README. */
  char input_filename_cpy[BUFSIZE];
  strcpy(input_filename_cpy, input_filename);
  char *filename = basename(input_filename_cpy);
  char output_filename[BUFSIZE];


  sprintf(output_filename, "costs_%s_%d.txt", filename, num_of_threads);
  FILE *output_costs_file = fopen(output_filename, "w");
  if (!output_costs_file) {
    printf("Error: couldn't output costs file");
    return -1;
  }

  fprintf(output_costs_file, "%d %d\n", dim_x, dim_y);
	int max=0;
  for (int i = 0; i < dim_y; i++) {
	  for (int j = 0; j < dim_x; j++) {
		  fprintf(output_costs_file, "%d ", costs[i*dim_x + j]);
		  if ((i*dim_x + j + 1) % dim_x == 0) fprintf(output_costs_file, "\n");
			if (costs[i*dim_x+j]>max) max=costs[i*dim_x+j];
	  }
  }
	fprintf(stdout, " Max: %d", max);
  fclose(output_costs_file);


  sprintf(output_filename, "output_%s_%d.txt", filename, num_of_threads);
  FILE *output_routes_file = fopen(output_filename, "w");
  if (!output_routes_file) {
    printf("Error: couldn't output routes file");
    return -1;
  }

  fprintf(output_routes_file, "%d %d\n", dim_x, dim_y);
  fprintf(output_routes_file, "%d\n", num_of_wires);

  for (int i = 0; i < num_of_wires; i++) {
	wire_t w = wires[i];
	int x1 = w.c1[0];
	int x2 = w.c2[0];
	int y1 = w.c1[1];
	int y2 = w.c2[1];
	//fprintf(stdout, "path = %d, bend = %d \n", w.path, w.bend);
	if (x1 == x2 || y1 == y2) fprintf(output_routes_file, "%d %d %d %d\n", x1, y1, x2, y2);
	else {
		if (w.path == 0) {
			if (w.bend == x1) 
				fprintf(output_routes_file, "%d %d %d %d %d %d\n", x1, y1, x1, y2, x2, y2);
			else if (w.bend == x2) 
				fprintf(output_routes_file, "%d %d %d %d %d %d\n", x1, y1, x2, y1, x2, y2);
			else fprintf(output_routes_file, "%d %d %d %d %d %d %d %d\n", x1, y1, 
				w.bend, y1, w.bend, y2, x2, y2);
		}
		else {
			if (w.bend == y1) 
				fprintf(output_routes_file, "%d %d %d %d %d %d\n", x1, y1, x2, y1, x2, y2);
			else if (w.bend == y2) 
				fprintf(output_routes_file, "%d %d %d %d %d %d\n", x1, y1, x1, y2, x2, y2);
			else fprintf(output_routes_file, "%d %d %d %d %d %d %d %d\n", x1, y1, 
				x1, w.bend, x2, w.bend, x2, y2);
		}
	}
  }

  fclose(output_routes_file);

  return 0;
}
