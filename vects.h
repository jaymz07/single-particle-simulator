#ifndef VECTS_H
#define VECTS_H

// Number of dimensions, 2 for x,y solution
#define NDIM 2

// ----------------Data Structures---------------
struct Vec {        // General Vector
    double r[NDIM];
};

struct Point {       // Phase space point
    struct Vec r;
    struct Vec v;
};

// ---------------Vector Functions----------------
double vec_mag_sqr(struct Vec v);
double vec_mag(struct Vec v);
struct Vec vec_sum(struct Vec a, struct Vec b);
struct Vec vec_diff(struct Vec a, struct Vec b);
struct Vec vec_scale(struct Vec v, double scale);
struct Vec vec_rotate(struct Vec v, double angle);
struct Vec vec_hat(struct Vec v);

struct Vec vec_zero();
struct Vec vec_copy(struct Vec v);
struct Vec vec_from_point(struct Point p);
struct Vec vec_from_array(double* arr);

// ---------------Point Functions----------------
double point_mag(struct Point p);
double point_mag_sqr(struct Point p);
struct Point point_sum(struct Point p1, struct Point p2);
struct Point point_diff(struct Point p1, struct Point p2);
struct Point point_rotate(struct Point p, double angle); // Only 2-D!!
struct Point point_scale(struct Point p, double scale);
struct Point point_copy(struct Point p);

#endif // VECTS_H
