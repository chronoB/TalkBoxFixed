#include <stdlib.h>
#include "calcAutoCoeff32.h"

void calcAutoCoeff32(int32_t *acf, int num_acf, int32_t *signal, int num_signal)
{
    int i, k, n_shift;
    int32_t max_value;
    int32_t abs_value;
    int64_t temp64;
    int32_t temp32;

    // integer base 2 logarithm
    n_shift = 0;
    for (i = 1; i < num_signal; i *= 2)
        n_shift++;
    n_shift = (n_shift + 1) / 2;

    // max(abs(signal))
    max_value = 0;
    for (i = 0; i < num_signal; i++)
    {
        abs_value = labs(signal[i]);

        if (max_value < abs_value)
            max_value = abs_value;
    }

    // number of leading signals of signal
    for (i = 0; i < 32; i++)
    {
        if (max_value >= 0x40000000)
            break;

        max_value = max_value << 1;
        n_shift--;
    }

    // normalization of signal
    if (n_shift > 0)
    {
        for (i = 0; i < num_signal; i++)
            signal[i] = signal[i] >> n_shift;
    }
    else
    {
        n_shift = -n_shift;
        for (i = 0; i < num_signal; i++)
            signal[i] = signal[i] << n_shift;
    }

    // acf[0]
    temp64 = 0;
    for (i = 0; i < num_signal; i++)
        temp64 += ((int64_t) signal[i] * signal[i]);
    temp32 = (int32_t) (temp64 >> 32);

    if (temp32 == 0)
    {
        acf[0] = 0x7FFFFFFF;
        for (int k = 1; k <num_acf; k++)
            acf[k] = 0;
        return;
    }

    // number of leading zeros of acf[0]
    for (i = 0; i < 32; i++)
    {
        if (temp32 >= 0x20000000)
            break;

        temp32 = temp32 << 1;
    }

    // 32 - nlz(acf[0])
    n_shift = 32 - i;

    // autocorrelation function
    for (k = 0; k <num_acf; k++)
    {
        temp64 = 0;
        for (i = 0; i < num_signal - k; i++)
        {
            temp64 += ((int64_t) signal[i + k] * signal[i]);
        }
        acf[k] = (int32_t) (temp64 >> n_shift);
    }

    // 1/acf[0] in 4.28 format, 5.59 / 1.31 = 4.28
    int32_t inv_acf0 = (int32_t) ((1LL << 59) / acf[0]);

    const int64_t max_acf = (1ll << 59)-1;

    // acf[i] = acf[i] / acf[0];
    for (k = 0; k < num_acf; k++)
    {
        temp64 = ((int64_t) acf[k] * inv_acf0);

        if (temp64 > max_acf)
            temp64 = max_acf;

        acf[k] = (int32_t) (temp64 >> 28);
    }
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