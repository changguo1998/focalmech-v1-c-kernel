#ifndef _KERNEL_H_
#define _KERNEL_H_

#include "types.h"

#ifdef GPU
__global__
#endif
void kernel(
#ifndef GPU
    int imisfit,
#endif
    int           n_trace,
    Trace*        tr,
    int           n_mt,
    MomentTensor* mt,
    int           n_loc,
    GFtrace*      gf,
    float*        pl2,
    int*          pshift,
    float*        sl2,
    int*          sshift,
    float*        pol,
    float*        psr);

#endif // _KERNEL_H_
