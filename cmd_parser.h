#include "ode_solvers.h" // Need Methods enum
#include "vects.h"  // NDIM is set here

// cmd_parser.h
#ifndef CMD_PARSER_H
#define CMD_PARSER_H

typedef struct {
    char *output;          // Path to output file
    int num_steps;         // Number of integration steps
    double time_step;      // time step
    double init_v[NDIM];   // initial velocity
    double init_r[NDIM];   // initial position

    double kill_range[2];  // We can use this to look for periodic orbits
    int kill_condition;
    int compute_closest_only;

    enum Method method;    // Solver method enum
    Integrator integrator; // Pointer to choice function from ode_solvers.c

    char forceStr[256];        // Choice of force function

    int save_skip;         // save on each save_skip time step
    int binary_output;     // Option to control format of output file

} CommandLineOptions;

enum Method get_method(char *option);

double  parse_double(char *option);
int     parse_int(char *option);
int parse_vec(char *option, double *vec);

void print_help();
void print_methods();

void process_arguments(int argc, char *argv[], CommandLineOptions *options);

#endif // CMD_PARSER_H
