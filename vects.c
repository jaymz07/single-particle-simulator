#include "vects.h"
#include <math.h>

struct Vec vec_zero() {
    struct Vec out;
    for (int i=0; i<NDIM; i++)
        out.r[i] = 0;
    return out;
}

double vec_mag_sqr(struct Vec v) {
    double out = 0;
    for (int i=0; i<NDIM; i++)
        out += v.r[i]*v.r[i];
    return out;
}

double vec_mag(struct Vec v) {
    return sqrt(vec_mag_sqr(v));
}

struct Vec vec_sum(struct Vec a, struct Vec b) {
    struct Vec out;
    for (int i=0; i<NDIM; i++)
        out.r[i] = a.r[i] + b.r[i];
    return out;
}

struct Vec vec_diff(struct Vec a, struct Vec b) {
    struct Vec out;
    for (int i=0; i<NDIM; i++)
        out.r[i] = a.r[i] - b.r[i];
    return out;
}

struct Vec vec_scale(struct Vec v, double scale) {
    struct Vec out;
    for (int i=0; i<NDIM; i++)
        out.r[i] = v.r[i]*scale;
    return out;
}

struct Vec vec_hat(struct Vec v) {
    double mag = vec_mag(v);
    return vec_scale(v, 1.0/mag);
}

struct Vec vec_rotate(struct Vec v, double angle) { // Only 2-D!!
    struct Vec out;
    double s=sin(angle), c=cos(angle);
    out.r[0] = v.r[0]*c - v.r[1]*s;
    out.r[1] = v.r[1]*s + v.r[1]*c;
    return out;
}

struct Vec vec_copy(struct Vec v) {
    struct Vec out;
    for (int i=0; i<NDIM; i++)
        out.r[i] = v.r[i];
    return out;
}

struct Vec vec_from_point(struct Point p) {
    struct Vec out;
    for (int i=0; i<NDIM; i++)
        out.r[i] = p.r.r[i];
    return out;
}

struct Vec vec_from_array(double* arr) {
    struct Vec out;
    for (int i=0; i<NDIM; i++)
        out.r[i] = arr[i];
    return out;
}

double point_mag(struct Point p) {
    return sqrt(point_mag_sqr(p));
}

double point_mag_sqr(struct Point p) {
    return vec_mag_sqr(p.r) + vec_mag_sqr(p.v);
}

struct Point point_scale(struct Point p, double scale) {
    struct Point out;
    out.r = vec_scale(p.r, scale);
    out.v = vec_scale(p.v, scale);
    return out;
}

struct Point point_sum(struct Point p1, struct Point p2) {
    struct Point out;
    out.r = vec_sum(p1.r, p2.r);
    out.v = vec_sum(p1.v, p2.v);
    return out;
}

struct Point point_diff(struct Point p1, struct Point p2) {
    struct Point out;
    out.r = vec_diff(p1.r, p2.r);
    out.v = vec_diff(p1.v, p2.v);
    return out;
}

struct Point point_rotate(struct Point p, double angle) { // Only 2-D!!
    struct Point out;
    double s=sin(angle), c=cos(angle);
    out.r.r[0] = p.r.r[0]*c - p.r.r[1]*s;
    out.r.r[1] = p.r.r[1]*s + p.r.r[1]*c;
    out.v.r[0] = p.v.r[0]*c - p.v.r[1]*s;
    out.v.r[1] = p.v.r[1]*s + p.v.r[1]*c;
    return out;
}

struct Point point_copy(struct Point p) {
    struct Point out;
    out.r = vec_copy(p.r);
    out.v = vec_copy(p.v);
    return out;
}
