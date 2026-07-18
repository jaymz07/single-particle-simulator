#include "ode_solvers.h" // Defines NDIM
#include "vects.h"

// ------------Integrators------------------
// maps input point p to new point using tstep=timestep, and accel(Point)=forceFunc(Point)
// written so that an arbitary number of dimensions can be used.

struct Point euler(struct Point p, double tstep, double time, struct Vec (*accel)(struct Point, double t)) {
    struct Point new;
    struct Vec a = accel(p,time);
    new.v = vec_sum(p.v, vec_scale(a,   tstep));
    new.r = vec_sum(p.r, vec_scale(p.v, tstep));
    return new;
}

struct Point cromer(struct Point p, double tstep, double time, struct Vec (*accel)(struct Point, double t)) {
    struct Point new;
    struct Vec a = accel(p,time);
    new.v = vec_sum(p.v, vec_scale(a,   tstep));
    new.r = vec_sum(p.r, vec_scale(new.v, tstep));
    return new;
}

struct Point newton(struct Point p, double tstep, double time, struct Vec (*accel)(struct Point, double t)) {
    struct Point new;
    new.r = vec_sum(p.r, vec_scale(p.v, tstep));
    struct Vec a = accel(new,time);
    new.v = vec_sum(p.v, vec_scale(a,   tstep));
    return new;
}

struct Point verlet_p(struct Point p, double tstep, double time, struct Vec(*accel)(struct Point, double t)) {
    struct Point new;
    double t_local=time;
    new.r = vec_sum(p.r, vec_scale(p.v, 0.5*tstep));
    t_local += tstep/2;
    struct Vec a = accel(new, t_local);    // a(r1)
    new.v = vec_sum(p.v, vec_scale(a, tstep));
    new.r = vec_sum(new.r, vec_scale(new.v, 0.5*tstep));
    // t_local += tstep/2;
    return new;
}

struct Point verlet_v(struct Point p, double tstep, double time, struct Vec(*accel)(struct Point, double t)) {
    struct Point new;
    double t_local=time;
    struct Vec a = accel(p, t_local);    // a(r)
    new.v = vec_sum(p.v, vec_scale(a,    0.5*tstep));
    new.r = vec_sum(p.r, vec_scale(new.v,  tstep));
    t_local += tstep;
    a = accel(new, t_local);    // a(r1)
    new.v = vec_sum(new.v, vec_scale(a, 0.5*tstep));
    return new;
}

struct Point rk2(struct Point p, double tstep, double time, struct Vec(*accel)(struct Point, double t)) {
    struct Point new;
    struct Point ra;
    struct Vec a_n = accel(p,time);
    new.r = vec_sum(vec_sum(p.r, vec_scale(p.v, tstep)), vec_scale(a_n, 0.5*tstep*tstep));
    ra.r = vec_sum(p.r, vec_scale(p.v, tstep));

    struct Vec a_n1 = accel(ra,time);  // a(rn + vn*dt)
    new.v = vec_sum(p.v, vec_scale(vec_sum(a_n, a_n1), 0.5*tstep));
    return new;
}

struct Point rkn4(struct Point p, double tstep, double time, struct Vec(*accel)(struct Point, double t)) {
    struct Point k1, k2, k3, k4;
    k1.v = accel(p,time);
    k1.r = vec_copy(p.v);
    k2.v = accel(point_sum(p, point_scale(k1, 0.5*tstep)), time);
    k2.r = vec_sum(p.v,       vec_scale(k1.v, 0.5*tstep));
    k3.v = accel(point_sum(p, point_scale(k2, 0.5*tstep)), time);
    k3.r = vec_sum(p.v,       vec_scale(k2.v, 0.5*tstep));
    k4.v = accel(point_sum(p, point_scale(k3, tstep)), time);
    k4.r = vec_sum(p.v,       vec_scale(k3.v, tstep));
    struct Point k;
    k = point_sum(k1, point_scale(k2, 2));
    k = point_sum(k,  point_scale(k3, 2));
    k = point_sum(k,  k4);
    return point_sum(p, point_scale(k, tstep/6));
}

struct Point t4_vv(struct Point p, double tstep, double time, struct Vec(*accel)(struct Point, double t)) {
    struct Point pn1, pn2;
    pn1 = verlet_v(p, tstep/2, time, accel);
    pn1 = verlet_v(pn1, tstep/2, time, accel);
    pn1 = point_scale(pn1, 4.0/3);

    pn2 = verlet_v(p, tstep, time, accel);
    pn2 = point_scale(pn2, -1.0/3);

    return point_sum(pn1, pn2);
}

struct Point t4_pv(struct Point p, double tstep, double time, struct Vec(*accel)(struct Point, double t)) {
    struct Point pn1, pn2;
    pn1 = verlet_p(p, tstep/2, time, accel);
    pn1 = verlet_p(pn1, tstep/2, time, accel);
    pn1 = point_scale(pn1, 4.0/3);

    pn2 = verlet_p(p, tstep, time, accel);
    pn2 = point_scale(pn2, -1.0/3);

    return point_sum(pn1, pn2);
}

struct Point fr_pv(struct Point p, double tstep, double time, struct Vec(*accel)(struct Point, double t)) {
    const double s = 1.2599210498948732; // cbrt(2)
    double s1=1.0/(2-s), s2=s/(2-s);
    double t_local=time;

    struct Point out;
    out = verlet_p(p,    s1*tstep, t_local, accel);
    t_local += tstep*s1;
    out = verlet_p(out, -s2*tstep, t_local, accel);
    t_local -= tstep*s2;
    out = verlet_p(out,  s1*tstep, t_local, accel);
    //t_local += tstep*s1;

    return out;
}
