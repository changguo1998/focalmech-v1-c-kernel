#ifndef _IO_H_
#define _IO_H_

#include "types.h"

void read_input(const char*    fname,
                int*           n_obs_traces,
                int*           n_mt,
                int*           n_srcloc,
                Trace**        traces,
                MomentTensor** mts,
                GFtrace**      gf_database);

void write_result(const char*  fname,
                  const int    n_obs_traces,
                  const int    n_mt,
                  const int    n_srcloc,
                  const float* misfit_pl2,
                  const int*   misfit_pshift,
                  const float* misfit_sl2,
                  const int*   misfit_sshift,
                  const float* misfit_polarity,
                  const float* misfit_psr);

#endif // _IO_H_
