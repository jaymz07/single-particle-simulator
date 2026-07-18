#include "vects.h"

#ifndef ODE_SOLVERS_H
#define ODE_SOLVERS_H

// Allow us to keep track of the method we use
enum Method {
    EULER,
    NEWTON,
    CROMER,
    VERLET,
    VERLET_V,
    RK2,
    RKN4,
    T4_VV,
    T4_PV,
    FR_PV
};

// Function pointer type we call Integrator
// Allows us to conveniently change methods
// Takes as argument a pointer to accel(Point, double) function
typedef struct Point (*Integrator)(struct Point, double, double, struct Vec (*)(struct Point, double t));

// Solver methods
struct Point euler(struct Point p, double tstep, double time, struct Vec (*accel)(struct Point, double t));
struct Point cromer(struct Point p, double tstep, double time, struct Vec (*accel)(struct Point, double t));
struct Point newton(struct Point p, double tstep, double time, struct Vec (*accel)(struct Point, double t));
struct Point verlet_p(struct Point p, double tstep, double time, struct Vec(*accel)(struct Point, double t));
struct Point verlet_v(struct Point p, double tstep, double time, struct Vec(*accel)(struct Point, double t));
struct Point rk2(struct Point p, double tstep, double time, struct Vec(*accel)(struct Point, double t));
struct Point rkn4(struct Point p, double tstep, double time, struct Vec(*accel)(struct Point, double t));
struct Point t4_vv(struct Point p, double tstep, double time, struct Vec(*accel)(struct Point, double t));
struct Point t4_pv(struct Point p, double tstep, double time, struct Vec(*accel)(struct Point, double t));
struct Point fr_pv(struct Point p, double tstep, double time, struct Vec(*accel)(struct Point, double t));

#endif // ODE_SOLVERS_H
