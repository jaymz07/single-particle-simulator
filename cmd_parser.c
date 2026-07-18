/* cmd_parser.c
 * This is a platform-dependent file.
 * It makes use of the getopt() function, which assumes
 * a POSIX compliant system.
 * Tested to be working with Ubuntu 22.04 and GCC 11.4.0
 */
#include "cmd_parser.h"
#include "ode_solvers.h" // Need Methods enum, NDIM
#include "data_utils.h"
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>

#define MAX_VEC 10

enum Method get_method(char *option) {
    if(strcmp(option, "newton") == 0)
        return NEWTON;
    if(strcmp(option, "euler") == 0)
        return EULER;
    if(strcmp(option, "cromer") == 0)
        return CROMER;
    if(strcmp(option, "verlet") == 0)
        return VERLET;
    if(strcmp(option, "verlet_v") == 0)
        return VERLET_V;
    if(strcmp(option, "rk2") == 0)
        return RK2;
    if(strcmp(option, "rkn4") == 0)
        return RKN4;
    if(strcmp(option, "t4_vv") == 0)
        return T4_VV;
    if(strcmp(option, "t4_pv") == 0)
        return T4_PV;
    if(strcmp(option, "fr_pv") == 0)
        return FR_PV;
    printf("Error! No such method: %s\n", option);
    print_methods();
    printf("Defaulting to Verlet.\n");
    return VERLET;
}

Integrator get_integrator(enum Method method) {
    Integrator out;
    switch (method) {
        case EULER:
            printf("Using Euler method.\n");
            out = euler;
            break;
        case NEWTON:
            printf("Using Newton method.\n");
            out = newton;
            break;
        case CROMER:
            printf("Using Cromer method.\n");
            out = cromer;
            break;
        case VERLET:
            printf("Using Position Verlet method.\n");
            out = verlet_p;
            break;
        case VERLET_V:
            printf("Using Velocity Verlet method.\n");
            out = verlet_v;
            break;
        case RK2:
            printf("Using RK2 method.\n");
            out = rk2;
            break;
        case RKN4:
            printf("Using RKN4 method.\n");
            out = rkn4;
            break;
        case T4_VV:
            printf("Using T4_vv method.\n");
            out = t4_vv;
            break;
        case T4_PV:
            printf("Using T4_pv method.\n");
            out = t4_pv;
            break;
        case FR_PV:
            printf("Using Forrest-Ruth method.\n");
            out = fr_pv;
            break;
        default:
            printf("Unrecognized Integration Method!\n");
            return NULL;
    }
    return out;
}

double parse_double(char *option) {
    char *endptr;
    double val;
    val = strtod(option, &endptr);
    // Check for errors in conversion
    if (*endptr != '\0') {
        fprintf(stderr, "Error: Invalid floating point number '%s'\n", option);
        exit(EXIT_FAILURE);
    }
    return val;
}

int parse_int(char *option) {
    char *endptr;
    long val;
    int int_val;
    val = strtol(option, &endptr, 10);
    if (endptr == option) {
        fprintf(stderr, "Error parsing Integer: %s\n", option);
        exit(EXIT_FAILURE);
    }
    int_val = (int) val;
    return int_val;
}

int parse_vec(char *option, double *vec) {
    int num_vals = split_doubles(option, vec, MAX_VEC);
    if(num_vals != NDIM) {
        printf("Error parsing vector: %s\n", option);
        return -1;
    }
    return num_vals;
}

void print_help() {
    //printf("HELP!\n");
    printf("Usage:\n");
    printf("Commamnd line options:\n");
    printf("\t-m [method]\n");
    printf("\t-o [output file]\n");
    printf("\t-n [num steps]\n");
    printf("\t-t [time step]\n");
    printf("\t-v [initial v_x],[initial v_y]\n");
    printf("\t-r [initial r_x],[initial r_y]\n");
    printf("\n");
}

void print_methods() {
    printf("Available methods:\nnewton\nrk2\nverlet\neuler\ncromer\n");
}

void process_arguments(int argc, char *argv[], CommandLineOptions *options) {
    int opt, val_int;

    // Initialize default options
    options->method = VERLET;
    options->num_steps = 5000;
    options->time_step = 0.01;
    options->output = NULL;
    options->binary_output=0;   // Default to text file output
    options->save_skip=1;       // Default to outputting every time-step
    options->kill_condition=0;  // Default to no kill condition
    options->compute_closest_only=0; // Default to output of entire trajectory
    strcpy(options->forceStr, "kepler");

    options->init_r[0] = 10.0;
    options->init_r[1] = 0.0;
    options->init_v[0] = 0.0;
    options->init_v[1] = 0.1;

    // Define short options
    const char *short_options = "hcm:n:t:o:v:r:s:b:k:f:";

    if(argc < 2) {
        print_help();
        exit(EXIT_SUCCESS);
    }
    while ((opt = getopt(argc, argv, short_options)) != -1) {
        switch (opt) {

            case 'h':  // Help
                print_help();
                exit(EXIT_SUCCESS);
                break;
            case 'c': // Compute closest approach
                options->compute_closest_only=1;
                break;

            // Solver options -------
            case 'm':  // Method Selection
                options->method = get_method(optarg);
                break;
            case 'n':  // Number of steps
                options->num_steps = parse_int(optarg);
                printf("Num steps set to %d\n", options->num_steps);
                break;
            case 't':  // Time num_step
                options->time_step = parse_double(optarg);
                printf("Time step set to %f\n", options->time_step);
                break;
            case 'f': // Choice of force function
                strcpy(options->forceStr, optarg);
                break;

            // Initial Conditions --------
            case 'r':  //Initial Position
                parse_vec(optarg, options->init_r);
//                 printf("r0:");
//                 for(int i=0;i<NDIM;i++)
//                     printf("%f\t", options->init_r[i]);
//                 printf("\n");
                break;
            case 'v':  //Initial Velocity
                parse_vec(optarg, options->init_v);
//                 printf("v0:");
//                 for(int i=0;i<NDIM;i++)
//                     printf("%f\t", options->init_v[i]);
//                 printf("\n");
                break;

            // Program Control ------
            case 'k': // Kill Condition
                parse_vec(optarg, options->kill_range);
                options->kill_condition = 1;
                break;

            // Output Control ------
            case 'o':  // File Output
                options->output = optarg;
                break;
            case 's':  // Save skip
                val_int = parse_int(optarg);
                if(val_int < 1)
                    printf("Invalid save_skip value!\n");
                options->save_skip = val_int;
                break;
            case 'b':  // Binary file output option
                val_int = parse_int(optarg);
                if(val_int == 0)
                    printf("Using text output mode.\n");
                else if(val_int == 1)
                    printf("Using binary output mode.\n");
                else {
                    printf("Invalid specifier for output format. Should be 0 or 1.\n");
                    printf("Defaulting to text output\n");
                    val_int = 0;

                }
                options->binary_output=val_int;
                break;
            case '?':  // Unknown option
                printf("Unknown option: %s\n", optarg);
                break;
            default:
                printf("Error in processing command line arguments.\n");
                break;
        }
    }
    options->integrator = get_integrator(options->method);
}
