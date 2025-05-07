#include <math.h>

#include "kernel.h"

#define LARGE_FLOAT 9999.0

// #define DEBUG 1

__device__
void cal_normalized_l2_gpu(float* nm,
                       int    l,
                       int    i1,
                       float* w1,
                       int    i2,
                       float* g11,
                       float* g22,
                       float* g33,
                       float* g12,
                       float* g13,
                       float* g23,
                       float  m11,
                       float  m22,
                       float  m33,
                       float  m12,
                       float  m13,
                       float  m23) {
    float sum = 0, s1 = 0, s2 = 0, t = 0, g = 0;
    int   i;
    for(i = 0; i < l; i++) {
        s1 += w1[i1 + i] * w1[i1 + i];
        g = g11[i2 + i] * m11 +
            g22[i2 + i] * m22 +
            g33[i2 + i] * m33 +
            g12[i2 + i] * m12 +
            g13[i2 + i] * m13 +
            g23[i2 + i] * m23;
        s2 += g * g;
    }
    s1 = sqrtf(s1);
    s2 = sqrtf(s2);
    t = s1 * s2;
    for(i = 0; i < l; i++) {
        g = g11[i2 + i] * m11 +
            g22[i2 + i] * m22 +
            g33[i2 + i] * m33 +
            g12[i2 + i] * m12 +
            g13[i2 + i] * m13 +
            g23[i2 + i] * m23;
        sum += w1[i1 + i] * g / t;
    }
    *nm = sqrtf(0.5 - 0.5 * sum);
}

__device__
void min_l2_gpu(Trace*        tr,
            MomentTensor* mt,
            GFtrace*      gf,
            int           iobs,
            int           isyn,
            int           wlen,
            int           maxshift,
            float*        misfit,
            int*          shift) {
    int   tshift;
    float tmisfit;
    *misfit = LARGE_FLOAT;
    *shift = 0;
#ifdef DEBUG
    printf("win: %d, obs: %d, syn: %d, maxshift: %d\n", wlen, iobs, isyn, maxshift);
#endif
    for(tshift = -maxshift; tshift < maxshift; tshift++) {
        cal_normalized_l2_gpu(&tmisfit, wlen, iobs, tr->obs, isyn + tshift,
            gf->g11, gf->g22, gf->g33, gf->g12, gf->g13, gf->g23,
            mt->m11, mt->m22, mt->m33, mt->m12, mt->m13, mt->m23);
        if(tmisfit < *misfit) {
            *misfit = tmisfit;
            *shift = tshift;
        }
    }
}

#define calc_syn(pgf, igf, pmt) \
    (pgf)->g11[(igf)] * (pmt)->m11 + \
    (pgf)->g22[(igf)] * (pmt)->m22 + \
    (pgf)->g33[(igf)] * (pmt)->m33 + \
    (pgf)->g12[(igf)] * (pmt)->m12 + \
    (pgf)->g13[(igf)] * (pmt)->m13 + \
    (pgf)->g23[(igf)] * (pmt)->m23

__global__
void kernel_gpu(
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
    float*        psr) {
    int itr, imt, iloc, igf, it, pmaxshift, smaxshift;
    int imisfit, bid, tid;
    bid = blockIdx.x;
    tid = threadIdx.x;
    imisfit = tid + bid * blockDim.x;

    if(imisfit >= (n_trace * n_mt * n_loc)) return;

    itr = imisfit % n_trace;
    it = (imisfit - itr) / n_trace;
    imt = it % n_mt;
    iloc = (it - imt) / n_mt;
    igf = iloc + n_loc * itr;

    Trace        *ptr = &tr[itr];
    MomentTensor *pmt = &mt[imt];
    GFtrace      *pgf = &gf[igf];

    pmaxshift = ptr->pwin_l;
    smaxshift = ptr->swin_l;

    // l2 norm
#ifdef DEBUG
    printf("L2 P\n");
#endif
    min_l2_gpu(ptr, pmt, pgf, ptr->pwin_i, pgf->p_i, ptr->pwin_l, pmaxshift, &pl2[imisfit], &pshift[imisfit]);

#ifdef DEBUG
    printf("L2 S\n");
#endif
    min_l2_gpu(ptr, pmt, pgf, ptr->swin_i, pgf->s_i, ptr->swin_l, smaxshift, &sl2[imisfit], &sshift[imisfit]);

    // polarity
#ifdef DEBUG
    printf("Polarity\n");
#endif
    pmaxshift = ptr->pwin_l / 8;
    pmaxshift = (pmaxshift > 0) ? pmaxshift : 5;
    float pol_obs = 0, pol_syn = 0, syn;
    for(it = 0; it < pmaxshift; it++) {
        pol_obs += ptr->obs[it + ptr->pwin_i];
        syn = calc_syn(pgf, it+pgf->p_i, pmt);
        pol_syn += syn;
    }
    pol[imisfit] = pol_obs * pol_syn >= 0.0 ? abs(pol_syn) : -abs(pol_syn);


    // ps ratio
#ifdef DEBUG
    printf("PSR\n");
#endif
    double obs_p_amp = 0, obs_s_amp = 0, syn_p_amp = 0, syn_s_amp = 0;
    for(it = 0; it < ptr->pwin_l; it++) {
        obs_p_amp += ptr->obs[it + ptr->pwin_i] * ptr->obs[it + ptr->pwin_i];
        syn = calc_syn(pgf, it+pgf->p_i+*pshift, pmt);
        syn_p_amp += syn * syn;
    }
    for(it = 0; it < ptr->swin_l; it++) {
        obs_s_amp += ptr->obs[it + ptr->swin_i] * ptr->obs[it + ptr->swin_i];
        syn = calc_syn(pgf, it+pgf->s_i+*sshift, pmt);
        syn_s_amp += syn * syn;
    }
    psr[imisfit] = log10((obs_s_amp * syn_p_amp) / (obs_p_amp * syn_s_amp));
}

#undef DEBUG
