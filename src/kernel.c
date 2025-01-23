#include <math.h>

#include "kernel.h"

// #define DEBUG 1

float cal_normalized_l2(int64_t l,
                        int64_t i1,
                        float*  w1,
                        int64_t i2,
                        float*  g11,
                        float*  g22,
                        float*  g33,
                        float*  g12,
                        float*  g13,
                        float*  g23,
                        float   m11,
                        float   m22,
                        float   m33,
                        float   m12,
                        float   m13,
                        float   m23) {
    double  sum = 0.0, s1 = 0.0, s2 = 0.0, t, g = 0;
    int64_t i;
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
    s1 = sqrt(s1);
    s2 = sqrt(s2);
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
    return sqrt(1.0 - sum);
}

void minl2(Trace        tr,
           MomentTensor mt,
           GFtrace      gf,
           int64_t      iobs,
           int64_t      isyn,
           int64_t      wlen,
           int64_t      maxshift,
           float*       misfit,
           int64_t*     shift) {
    int64_t tshift;
    float   tmisfit;
    *misfit = 9999.0;
    *shift = 0;
#ifdef DEBUG
    printf("win: %lld, obs: %lld, syn: %lld, maxshift: %lld\n", wlen, iobs, isyn, maxshift);
#endif
    for(tshift = -maxshift; tshift < maxshift; tshift++) {
        tmisfit = cal_normalized_l2(wlen, iobs, tr.obs, isyn + tshift,
            gf.g11, gf.g22, gf.g33, gf.g12, gf.g13, gf.g23, mt.m11, mt.m22, mt.m33, mt.m12, mt.m13, mt.m23);
        if(tmisfit < *misfit) {
            *misfit = tmisfit;
            *shift = tshift;
        }
    }
}

void kernel(Trace        tr,
            MomentTensor mt,
            GFtrace      gf,
            float*       pl2,
            int64_t*     pshift,
            float*       sl2,
            int64_t*     sshift,
            float*       pol,
            float*       psr) {
    int64_t it, pmaxshift, smaxshift;
    pmaxshift = tr.pwin_l;
    smaxshift = tr.swin_l;

    // l2 norm
#ifdef DEBUG
    printf("L2 P\n");
#endif
    minl2(tr, mt, gf, tr.pwin_i, gf.p_i, tr.pwin_l, pmaxshift, pl2, pshift);

#ifdef DEBUG
    printf("L2 S\n");
#endif
    minl2(tr, mt, gf, tr.swin_i, gf.s_i, tr.swin_l, smaxshift, sl2, sshift);

    // polarity
#ifdef DEBUG
    printf("Polarity\n");
#endif
    pmaxshift = tr.pwin_l / 8;
    pmaxshift = (pmaxshift > 0) ? pmaxshift : 5;
    double pol_obs = 0, pol_syn = 0, syn;
    for(it = 0; it < pmaxshift; it++) {
        pol_obs += tr.obs[it + tr.pwin_i];
        syn = gf.g11[it + gf.p_i] * mt.m11 +
            gf.g22[it + gf.p_i] * mt.m22 +
            gf.g33[it + gf.p_i] * mt.m33 +
            gf.g12[it + gf.p_i] * mt.m12 +
            gf.g13[it + gf.p_i] * mt.m13 +
            gf.g23[it + gf.p_i] * mt.m23;
        pol_syn += syn;
    }
    *pol = pol_obs * pol_syn >= 0.0 ? abs(pol_syn) : -abs(pol_syn);


    // ps ratio
#ifdef DEBUG
    printf("PSR\n");
#endif
    double obs_p_amp = 0, obs_s_amp = 0, syn_p_amp = 0, syn_s_amp = 0;
    for(it = 0; it < tr.pwin_l; it++) {
        obs_p_amp += tr.obs[it + tr.pwin_i] * tr.obs[it + tr.pwin_i];
        syn = gf.g11[it + gf.p_i + *pshift] * mt.m11 +
            gf.g22[it + gf.p_i + *pshift] * mt.m22 +
            gf.g33[it + gf.p_i + *pshift] * mt.m33 +
            gf.g12[it + gf.p_i + *pshift] * mt.m12 +
            gf.g13[it + gf.p_i + *pshift] * mt.m13 +
            gf.g23[it + gf.p_i + *pshift] * mt.m23;
        syn_p_amp += syn * syn;
    }
    for(it = 0; it < tr.swin_l; it++) {
        obs_s_amp += tr.obs[it + tr.swin_i] * tr.obs[it + tr.swin_i];
        syn = gf.g11[it + gf.s_i + *sshift] * mt.m11 +
            gf.g22[it + gf.s_i + *sshift] * mt.m22 +
            gf.g33[it + gf.s_i + *sshift] * mt.m33 +
            gf.g12[it + gf.s_i + *sshift] * mt.m12 +
            gf.g13[it + gf.s_i + *sshift] * mt.m13 +
            gf.g23[it + gf.s_i + *sshift] * mt.m23;
        syn_s_amp += syn * syn;
    }
    *psr = log10((obs_s_amp * syn_p_amp) / (obs_p_amp * syn_s_amp));
}

#undef DEBUG
