/*
    Kalign - a multiple sequence alignment program

    Copyright 2006, 2019 Timo Lassmann

    This file is part of kalign.

    Kalign is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.

*/

#include "euclidean_dist.h"
#include "rng.h"
#include <xmmintrin.h>
#include <immintrin.h>
#include "float.h"
/* These functions were taken from:  */
/* https://stackoverflow.com/questions/6996764/fastest-way-to-do-horizontal-float-vector-sum-on-x86 */
float hsum256_ps_avx(__m256 v);
float hsum_ps_sse3(__m128 v);

int edist_serial(const float* a,const float* b,const int len, float* ret)
{
        int i;
        float d = 0.0f;
        float t;

        for(i = 0; i < len;i++){
                t = (a[i] - b[i]);
                d += t *t;
        }

        *ret = sqrtf(d);
        return OK;
}


int edist_serial_d(const double* a,const double* b,const int len, double* ret)
{
        int i;
        double d = 0.0f;
        double t;

        for(i = 0; i < len;i++){
                t = (a[i] - b[i]);
                d += t *t;
        }

        *ret = sqrt(d);
        return OK;
}

#ifdef HAVE_AVX2

int edist_256(const float* a,const float* b, const int len, float* ret)
{

        float d = 0.0f;
        register int i;
        __m256 xmm1;// = _mm256_load_ps(a);
        __m256 xmm2;// = _mm256_load_ps(b);
        __m256 r = _mm256_set1_ps(0.0f);
        for(i = 0;i < len;i+=8){
                xmm1 = _mm256_load_ps(a);
                xmm2 = _mm256_load_ps(b);

                xmm1 =  _mm256_sub_ps(xmm1, xmm2);

                xmm1 = _mm256_mul_ps(xmm1, xmm1);

                r = _mm256_add_ps(r, xmm1);
                a+=8;
                b+=8;
        }
        d = hsum256_ps_avx(r);

        *ret = sqrtf(d);
        return OK;
}




float hsum256_ps_avx(__m256 v)
{
        __m128 vlow  = _mm256_castps256_ps128(v);
        __m128 vhigh = _mm256_extractf128_ps(v, 1); // high 128
        vlow  = _mm_add_ps(vlow, vhigh);     // add the low 128
        return hsum_ps_sse3(vlow);         // and inline the sse3 version, which is optimal for AVX
        // (no wasted instructions, and all of them are the 4B minimum)
}

float hsum_ps_sse3(__m128 v)
{
        __m128 shuf = _mm_movehdup_ps(v);        // broadcast elements 3,1 to 2,0
        __m128 sums = _mm_add_ps(v, shuf);
        shuf        = _mm_movehl_ps(shuf, sums); // high half -> low half
        sums        = _mm_add_ss(sums, shuf);
        return        _mm_cvtss_f32(sums);
}

#endif
