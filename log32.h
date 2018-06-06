#ifndef __LOG32__
#define __LOG32__

#include <stdint.h>

//------------------------------------------------------------------------------
// number of leading zeros (signed)

#if ( __ADSPBLACKFIN__ )

inline int nlzs(uint32_t x)
{
    int8_t n;
    asm( "%0.l = signbits %1;" : "=d" (n) : "d" (x) );
    return n;
}

#elif( __xcore__ )

inline int nlzs(uint32_t x)
{
    int n;
    n = __builtin_clz(x) - 1;
    return n;
}

#else

inline int nlzs(uint32_t x)
{
    // Henry S. Warren, "Hacker's Delight (2nd Edition)", Addison-Wesley, 2012.
    // number of leading zeros, binary search, p. 99.

    int n;
    if (x == 0) return 32;

    n = 0;
    if (x <= 0x0000FFFF) {n = n +16; x = x <<16;}
    if (x <= 0x00FFFFFF) {n = n + 8; x = x << 8;}
    if (x <= 0x0FFFFFFF) {n = n + 4; x = x << 4;}
    if (x <= 0x3FFFFFFF) {n = n + 2; x = x << 2;}
    if (x <= 0x7FFFFFFF) {n = n + 1;}
    return n - 1;
}

#endif

//------------------------------------------------------------------------------
// 32-bit logarithm
#define L_20LOG10   0x06054609  /* y = 20*log10(x) */
#define L_10LOG10   0x0302A305  /* y = 10*log10(x) */
#define L_LOG10     0x004D104D  /* y = log10(x)    */
#define L_LOG       0x00B17218  /* y = log(x)      */
#define L_LOG2      0x01000000  /* y = log2(x)     */

inline int32_t log32(int32_t in_lin, int32_t conv_coeff = L_20LOG10)
{
    // log-function: input in_lin (Q1.31 format), output out_log (Q16.16 format).
    // coeffs[] contains the coefficients of the Taylor series (Q2.14 format).
    // conv_coeff (Q8.24 format) for conversion to log, log10...
    // (c) Hagen Jaeger, Uwe Simmer, April 2014.

    static int16_t coeffs[] = {23637, -11819, 7879, -5909, 4727,
                               -3940, 3377, -2955, 2626, -2364};

    int taylor_cnt, shift_cnt, taylor_deg = 10;
    int16_t product, x1;
    int32_t out_log;

    shift_cnt = nlzs(in_lin);       // number of leading zeros

    in_lin = in_lin << shift_cnt;   // normalization

    out_log = -shift_cnt << 16;     // conversion to Q16.16

    x1 = in_lin >> 16;              // conversion from Q1.31 to Q1.15
    x1 = x1 + 0x8000;               // (in_lin - 1)

    product = x1;                   // x^1,	Q1.15 format


    // Taylor series of log(x)
    for (taylor_cnt = 0; taylor_cnt < taylor_deg; taylor_cnt++)
    {
        out_log += ((int32_t) product * (int32_t) coeffs[taylor_cnt]) >> 13;
        // Q1.15 format * Q2.14 = Q3.29 format, Q3.29 >> 13 = Q16.16 format

        product = ((int32_t) product * (int32_t) x1) >> 15;
        // power of x for Taylor series: x, x^2, x^3, etc.
        // Q1.15 format * Q1.15 = Q2.30, Q2.30 >> 15 = Q1.15 format
    }

    // product of out_log (Q16.16 format) and conv_coeff (Q8.24 format)
    out_log = ((int64_t) out_log * (int64_t) conv_coeff) >> 24;

    return out_log;
}

//------------------------------------------------------------------------------
// 32-bit exponential function
#define E_20LOG10   0x002A854B  /* y = 10^(x/20) */
#define E_10LOG10   0x00550A97  /* y = 10^(x/10) */
#define E_LOG10     0x035269E1  /* y = 10^x      */
#define E_LOG       0x01715476  /* y = e^x       */
#define E_LOG2      0x01000000  /* y = 2^x       */

inline int32_t exp32(int32_t in_log, int32_t conv_coeff = E_20LOG10)
{
    // e^x-Funktion: input in_log (Q16.16 format), output out_lin (Q1.31 format).
    // coeffs[] contains the coefficients of the Taylor series (Q0.16 format).
    // conv_coeff (Q8.24 format) for conversion to anti-log, anti-log10...
    // (c) Hagen Jaeger, Uwe Simmer, April 2014.

    static int16_t coeffs[] = {10923, 2731, 546};
    int shift_cnt, taylor_cnt, taylor_deg = 3;
    int16_t product, x, hw;
    uint16_t lw;
    int16_t ln2 = 0x58B9;
    int32_t temp32;
    uint32_t out_lin;

    // product of in_log (Q16.16 format) and conv_coeff (Q8.24 format)
    in_log = ((int64_t) in_log * conv_coeff) >> 24;

    hw = in_log >> 16;          // extraction of high word, Q1.15 format
    lw = in_log & 0xFFFF;       // extraction of low  word, Q0.16 format

    shift_cnt = -hw;            // number of right shifts

    if (shift_cnt > 31)
        return 0;

    x = lw >> 1;                // conversion to Q1.15 format

    out_lin = 0x7FFFFFFF;                   // out = 1
    temp32 = ((int32_t) x * (int32_t) ln2); // x' = x*ln(2)
    out_lin += temp32 << 1;                 // out = 1 + x'
    x = (int16_t) (temp32 >> 15);           // x': 1.15
    temp32 = ((int32_t) x * (int32_t) x);   // x' * x'
    out_lin += temp32;                      // out = 1 + x' + x'*x'/2
    product = (int16_t) (temp32 >> 15);

    // Taylor series of exp(x), starting at order 3
    for (taylor_cnt = 0; taylor_cnt < taylor_deg; taylor_cnt++)
    {
        product = (int16_t) (((int32_t) product * (int32_t) x) >> 15);
        // power of x for Taylor series: x^3, x^4, x^5
        // Q1.15 format * Q1.15 = Q2.30, Q2.30 >> 15 = Q1.15 format

        out_lin += (int32_t) product * (int32_t) coeffs[taylor_cnt];
        // Q1.15 format * Q0.16 format = Q1.31 format
    }

    // out = taylor(frac(in)) * 2^int(in)
    out_lin = out_lin >> shift_cnt;

    return out_lin;
}

#endif  // __LOG32__

//--------------------- License -----------------------------------------------

// Copyright (c) 2014-2016 Hagen Jaeger, Uwe Simmer,
// Institute for Hearing Technology and Audiology,
// Jade University of Applied Sciences Oldenburg.

// Permission is hereby granted, free of charge, to any person obtaining
// a copy of this software and associated documentation files
// (the "Software"), to deal in the Software without restriction, including
// without limitation the rights to use, copy, modify, merge, publish,
// distribute, sublicense, and/or sell copies of the Software, and to
// permit persons to whom the Software is furnished to do so, subject
// to the following conditions:
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.

// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
// EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
// OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
// IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
// CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
// TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
// SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
