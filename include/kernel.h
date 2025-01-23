#ifndef _KERNEL_H_
#define _KERNEL_H_

#include <stdint.h>
#include "types.h"

void kernel(Trace         tr,
            MomentTensor  mt,
            GFtrace       gf,
            float*        pl2,
            int64_t*      pshift,
            float*        sl2,
            int64_t*      sshift,
            float*        pol,
            float*        psr);

#endif // _KERNEL_H_
