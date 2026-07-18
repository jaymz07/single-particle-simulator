#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "data_utils.h"
#include "cmd_parser.h"
#include "ode_solvers.h"
#include "vects.h"

//-------------Force Definition-------------------------
// This is the Kepler problem
struct Vec kepler(struct Point p, double t) {
    double r3 = pow(vec_mag(p.r), 3);
    //double r3 = 1;     // For spring problem
    struct Vec out;
    out = vec_scale(p.r, -1.0/r3);
    return out;
}

// Fixed two center gravitational field
struct Vec grav_two_center(struct Point p, double t) {
    struct Vec r1, r2;
    r1.r[0] = -0.5;
    r1.r[1] =  0.0;
    r2.r[0] =  0.5;
    r2.r[1] =  0.0;

    struct Vec delta1 = vec_diff(p.r, r1);
    struct Vec delta2 = vec_diff(p.r, r2);
    double d1=vec_mag(delta1);
    double d2=vec_mag(delta2);

    struct Vec f1, f2;
    f1 = vec_scale(delta1, -0.5/pow(d1, 3));
    f2 = vec_scale(delta2, -0.5/pow(d2, 3));

    return vec_sum(f1, f2);

}

// Rotating two center gravitational field
struct Vec grav_two_center_rot(struct Point p, double t) {
    struct Vec r1, r2;
    r1.r[0] =  -0.5*cos(t);
    r1.r[1] =  -0.5*sin(t);
    r2.r[0] =   0.5*cos(t);
    r2.r[1] =   0.5*sin(t);

    struct Vec delta1 = vec_diff(p.r, r1);
    struct Vec delta2 = vec_diff(p.r, r2);
    double d1=vec_mag_sqr(delta1);
    double d2=vec_mag_sqr(delta2);

    struct Vec f1, f2;
    f1 = vec_scale(delta1, -0.5/pow(d1, 3.0/2));
    f2 = vec_scale(delta2, -0.5/pow(d2, 3.0/2));

    return vec_sum(f1, f2);

}

struct Vec grav_two_center_rotframe(struct Point p, double t) {
    struct Vec r1, r2;
    r1.r[0] =  -0.5;
    r1.r[1] =   0;
    r2.r[0] =   0.5;
    r2.r[1] =   0;

    struct Vec delta1 = vec_diff(p.r, r1);
    struct Vec delta2 = vec_diff(p.r, r2);
    double d1=vec_mag_sqr(delta1);
    double d2=vec_mag_sqr(delta2);

    struct Vec f1, f2, frot;
    f1 = vec_scale(delta1, -0.5/pow(d1, 3.0/2));
    f2 = vec_scale(delta2, -0.5/pow(d2, 3.0/2));

    frot = vec_sum(f1, f2);
    frot = vec_sum(frot, vec_scale(p.r, 1));
    struct Vec vp;
    vp.r[0] =  -2*p.v.r[1];
    vp.r[1] =   2*p.v.r[0];
    return vec_sum(frot, vp);
}

struct Vec forceFunc(struct Point p, double t) {
    //return kepler(p, t);
    //return grav_two_center(p, t);
    return grav_two_center_rot(p, t);
}

// ------------------Main------------------------------
int main(int argc, char *argv[]) {
    char* filename_out = "out/output.txt";
    Integrator integrator;
    CommandLineOptions options;
    struct Vec (*forcePtr)(struct Point, double);
    double saveOutput[2];

    printf("-------------------------------------\n");
    process_arguments(argc, argv, &options);

    // Choice of method is taken from command line arguments
    integrator = options.integrator;

    if(strcmp(options.forceStr, "circular_restricted") == 0)
        forcePtr = grav_two_center_rot;
    else if(strcmp(options.forceStr, "rot_frame") == 0)
        forcePtr = grav_two_center_rotframe;
    else if(strcmp(options.forceStr, "fixed_two_center") == 0)
        forcePtr = grav_two_center;
    else if(strcmp(options.forceStr, "kepler") == 0)
        forcePtr = kepler;
    else {
        printf("Unrecognized force option: %s\n", options.forceStr);
        strcpy(options.forceStr,"circular_restricted");
        forcePtr = grav_two_center_rot;
    }
    printf("Force mode set to: %s\n", options.forceStr);



    if(options.output == NULL) {
        printf("No output file specified!\n");
        //print_help();
        // return 1;
    }
    else
        filename_out = options.output;
    printf("Output will be written to: %s\n", filename_out);

    // -----------Parameters---------------

    int numSteps = options.num_steps;     // Number of points to calculate
    double tstep = options.time_step;     // Time step

    double t_start = 0.0;                 // Starting time
    double t_end = numSteps*tstep;        // Ending time

    int save_skip = options.save_skip;

    // -----------Program Start---------------
    printf("Integrating from ti=%.1f to tf=%.1f\n", t_start, t_end);

    // Use linspace function to initialize array
    //double* t = linspace(0, t_end, numSteps);

    // Dynamically initialize the output arrays
    double* t;
    double* x;
    double* y;
    double* vx;
    double* vy;

    double t_global = t_start;    // Time variable

    struct Point traj; // Current trajectory in phase space
    struct Point traj_0;
    struct Vec r0;

    r0 = vec_from_array(options.init_r);

    traj.r = vec_from_array(options.init_r);
    traj.v = vec_from_array(options.init_v);
    traj_0 = point_copy(traj);

/*  The proceeding double dot structure is perhaps the result of
    too much abstraction. I choose to leave it for now, as otherwise the
    2-D nature of the problem is not apparent here and would be
    tucked away in the NDIM constant contained in the "vects.h" file. */

    // Perform the selected method
    int i, j;
    struct Vec drvec;
    double dr;
    printf("Integrating.");
    if(! options.compute_closest_only) {

        t = initializeArray(numSteps);
        x = initializeArray(numSteps);
        y = initializeArray(numSteps);
        vx = initializeArray(numSteps);
        vy = initializeArray(numSteps);

        t[0] = t_global;
        x[0] = traj.r.r[0];  // Store the initial conditions in the output array
        y[0] = traj.r.r[1];
        vx[0] = traj.v.r[0];
        vy[0] = traj.v.r[1];

        for (i = 0; i < numSteps-1; i++) {
            for(j = 0; j < save_skip; j++) {
                traj = integrator(traj, tstep/save_skip, t_global, forcePtr);
                t_global += tstep/save_skip;
            }
            t[i+1] = t_global;
            x[i+1] = traj.r.r[0];
            y[i+1] = traj.r.r[1];
            vx[i+1] = traj.v.r[0];
            vy[i+1] = traj.v.r[1];
            if(i > 500 && options.kill_condition) {
                drvec = vec_diff(traj.r, r0);
                dr = vec_mag_sqr(drvec);
                if(dr > options.kill_range[0] &&
                        dr < options.kill_range[1]) {
                    break;
                }
            }
        }
    }
    else {
        double r, r_closest, t_closest=0;
        struct Point delta, p_closest;

        // Integrate some time to get away from the initial position
        for (i = 0; i < 1000; i++) { // Skip first N points
            traj = integrator(traj, tstep, t_global, forceFunc);
            t_global += tstep;
        }
        delta = point_diff(traj, traj_0);
        r = point_mag_sqr(delta);
        r_closest = r;
        // Start looking for the closest approach to where we started.
        for (i = i; i < numSteps-1; i++) {
            traj = integrator(traj, tstep, t_global, forceFunc);
            t_global += tstep;

            delta = point_diff(traj, traj_0);
            r = point_mag_sqr(delta);
            if(r < r_closest) {
                r_closest = r;
                t_closest = t_global;
                p_closest = point_copy(traj);
            }
        }
        //printf("CPA: %.18e\n", r_closest);
        //printf("CPT: %.18e\n", t_global);
        saveOutput[0] = r_closest;
        saveOutput[1] = t_closest;
        printf("rc: %f\t\n", sqrt(r_closest));
        //printf("dt: %f\t\n", tstep);
        writeArray_bin(saveOutput, 2, filename_out);

        return 0;
    }

    // Define output file structure
    int numColumns = 4;
    double* dataArray[] = {x, y, vx, vy};

    // Write the time and output arrays to a text file
    printf("Writing output file...\n");

    // Here i=numsteps-1, unless loop breaks.
    if(options.binary_output)
        writeFile_bin(t, dataArray, i-1, numColumns, filename_out);
    else
        writeFile(t, dataArray, i-1, numColumns, filename_out);
    printf("Complete.\n");

    // Free dynamically allocated memory
    free(t);
    free(x);
    free(y);
    free(vx);
    free(vy);

    return 0;
}
