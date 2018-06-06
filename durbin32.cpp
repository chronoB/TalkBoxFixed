#include <stdlib.h>
#include "durbin32.h"

/*---------------------------------------------------------------------------*\
|   Fixed-Point Version of the Durbin Algorithm                               |
|                                                                             |
|   Authors:                                                                  |
|       Finn Bayer, Christoph Eike, Uwe Simmer, 24. Jun. 2016.                |
|                                                                             |
|   Reference:                                                                |
|   [1] Robert Bristow-Johnson,                                               |
|       Fixed-point implementation of levinson durbin algorithm,              |
|       comp.dsp, 04.01.2011                                                  |
\*---------------------------------------------------------------------------*/

#define N 128

int32_t durbin32(int32_t *r, int32_t *a, int n, int fractional_digits,
                 int32_t k_max)
{
                            // r, k_max: 1.31 format
                            // a: 8.24 format
    int32_t a_temp[N],      // 8.24 format
            ki,             // 8.24 format
            alpha;          // 1.31 format
    int64_t epsilon;        // 9.55 format
    int32_t temp32;
    int i, j;

    /* n <= N = constant */
    if (n > N)
    {
        return 0;
    }

    // temp32 = 1.0
    temp32 = 1L << fractional_digits;
    k_max = (int32_t) (((int64_t) k_max * temp32) >> 31);

    for (i = 0; i < n; i++)
    {
        a[i] = 0;
    }

    alpha = r[0];

    for (i = 0; i < n; i++)
    {
        /* epsilon = a[0] * r[i]; */
        epsilon = ((int64_t) r[i+1]) << fractional_digits;
        for (j = 0; j<i; j++)
        {
            epsilon += (int64_t) a[j] * r[i - j];
        }

        ki = (int32_t) (-epsilon / alpha);

        if (labs(ki) > k_max)
        {
            return alpha;
        }

        a[i] = ki;  // 8.24 format

        temp32 = (0x7FFFFFFF - (int32_t) (((int64_t) ki * ki) >> (2 * fractional_digits - 31)));

        alpha = ((int64_t) alpha * temp32) >> 31;

        for (j = 0; j<i; j++)
        {
            /* update a[] array into temporary array */
            a_temp[j] = a[j] + (int32_t) (((int64_t) ki * a[i - j - 1]) >> fractional_digits);
        }

        for (j = 0; j<i; j++)
        {
            /* update a[] array */
            a[j] = a_temp[j];
        }
    }

    return alpha;
}

//--------------------- License ------------------------------------------------

// Copyright (c) 2016 Finn Bayer, Christoph Eike, Uwe Simmer

// Permission is hereby granted, free of charge, to any person obtaining 
// a copy of this software and associated documentation files 
// (the "Software"), to deal in the Software without restriction, 
// including without limitation the rights to use, copy, modify, merge, 
// publish, distribute, sublicense, and/or sell copies of the Software, 
// and to permit persons to whom the Software is furnished to do so, 
// subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included 
// in all copies or substantial portions of the Software.

// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, 
// EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF 
// MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. 
// IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY 
// CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, 
// TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE 
// SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

//------------------------------------------------------------------------------
