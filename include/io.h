#ifndef _IO_H_
#define _IO_H_

#include "types.h"

void read_input(char*          fname,
                int64_t*       n_obs_traces,
                int64_t*       n_mt,
                int64_t*       n_srcloc,
                Trace**        traces,
                MomentTensor** mts,
                GFtrace**      gf_database);

void write_result(char*          fname,
                  const int64_t  n_obs_traces,
                  const int64_t  n_mt,
                  const int64_t  n_srcloc,
                  const float*   misfit_pl2,
                  const int64_t* misfit_pshift,
                  const float*   misfit_sl2,
                  const int64_t* misfit_sshift,
                  const float*   misfit_polarity,
                  const float*   misfit_psr);

#endif // _IO_H_
