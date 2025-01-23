#ifndef _TYPES_H_
#define _TYPES_H_

#include <stdio.h>
#include <stdint.h>
#include <float.h>

// enum DeviceType { cpu, gpu };

typedef struct _Trace {
    int64_t   n_samples, ref_i, pwin_i, pwin_l, swin_i, swin_l;
    float* obs;
} Trace;

void read_trace(Trace* tr, FILE* fp);
void write_trace(const Trace* tr, FILE* fp);
void free_trace(Trace* tr);

typedef struct _MomentTensor {
    float m11, m22, m33, m12, m13, m23;
} MomentTensor;

void read_mt(MomentTensor* mt, FILE* fp);
void write_mt(const MomentTensor* mt, FILE* fp);

typedef struct _GFtrace {
    int64_t   n_samples, p_i,  s_i;
    float *g11, *     g22, *g33, *g12, *g13, *g23;
} GFtrace;

void read_gf_trace(GFtrace* gf, FILE* fp);
void write_gf_trace(const GFtrace* gf, FILE* fp);
void free_gf_trace(GFtrace* gf);

#endif // _TYPES_H_
