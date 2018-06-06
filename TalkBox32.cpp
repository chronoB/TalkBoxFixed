#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "TalkBox32.h"
#include "calcAutoCoeff32.h"
#include "durbin32.h"
#include "lpcFilter32.h"
#include "log32.h"

#define M_PI    3.14159265358979323846

const int32_t k_max = (int32_t) (0.99 * 0x7FFFFFFF);

/* a tunable high-pass filter based on a first order allpass */

inline int32_t highpass32(int32_t in, int32_t coeff, int32_t *mem)
{
    // coeff in 1.31 format
    int64_t temp64;
    int32_t out;

    in = in >> 1;

    temp64 = (int64_t) coeff * (in - mem[1]);

    out = (int32_t) (temp64 >> 31); // quantization
    out += mem[0];

    mem[0] = in;                    // non-recursive state
    mem[1] = out;                   // recursive state

    return (in - out);
}

TalkBox32::TalkBox32(double fs)
{
    this->fs = fs;

    // parameter for smoothing
    setSmoothingTime(0.03f);

    // gate off
    gate_level = 0;

    // integer base 2 logarithm of block_length
    n_shift_block = 0;
    for (int i=1; i<block_length; i*=2)
        n_shift_block++;

    // integer base 2 logarithm of memory_rms_size
    n_shift_memory = 0;
    for (int i=1; i<memory_rms_size; i*=2)
        n_shift_memory++;

    // high pass design
    double ftan = tan(M_PI * 20000. / fs);
    high_pass_coeff = (int32_t) ((ftan-1) / (ftan+1) * 0x7FFFFFFF);

    // set states to null
    resetStates();

    sample_buffer = input_buffer0;
    block_buffer  = input_buffer1;

    acf_index = 0;
}

TalkBox32::~TalkBox32(void)
{
}

void TalkBox32::process(int32_t samples[])
{
    int32_t temp32;

    // synthesizer signal
    temp32 = samples[0];

    // input * gain
    temp32 = ((int64_t) error_gain * temp32) >> 31;

    // input * voice_rms
    temp32 = ((int64_t) voice_rms * temp32) >> 31;

    // all-pole filter
    std::unique_lock<std::mutex> locker(a_coeff_mutex, std::defer_lock);
    locker.lock();

    samples[0] = lpcFilter32(temp32, a32, memory_lpc, num_coeffs, fractional_digits);

    locker.unlock();

    // voice signal
    sample_buffer[buffer_position++] = samples[1];

    if (buffer_position >= block_length)
    {
        buffer_position = 0;

        if (block_ready == true)
            printf("timing error\n");

        // swap buffer
        int32_t *tmp_ptr = block_buffer;
        block_buffer = sample_buffer;
        sample_buffer = tmp_ptr;

        block_ready = true;
    }
}

void TalkBox32::calculateLPCcoefficients(void)
{
    int32_t temp32;
    int32_t abs_voice;
    int32_t error_power32;

    // new input block available?
    if (block_ready == false)
        return;

    abs_voice = 0;
    for (int i=0; i<block_length; i++)
    {
        temp32 = block_buffer[i];

        // voise rms
        abs_voice += (labs(temp32) >> n_shift_block);

        // high pass
        temp32 = highpass32(temp32, high_pass_coeff, memory_hp);

        block_buffer[i] = temp32;
    }

    // RMS (FIR)
    for (int i = memory_rms_size - 1; i > 0; i--)
    {
        memory_rms32[i] =  memory_rms32[i - 1];
    }
    memory_rms32[0] = abs_voice;

    voice_rms = 0;
    for (int i = 0; i < memory_rms_size; i++)
    {
        voice_rms += (memory_rms32[i] >> n_shift_memory);
    }

    if (voice_rms < (1L << 29))
        voice_rms <<= 2;
    else
        voice_rms = 0x7FFFFFFF;

    if (voice_rms < gate_level)     // gate
    {
        voice_rms = 0;
    }

    calcAutoCoeff32(acf32[acf_index], num_coeffs+1, block_buffer, block_length);

    // averaging of acfs
    for (int i = 0; i < num_coeffs + 1; i++)
        acf32[acf_index][i] = (acf32[0][i] >> 2) + (acf32[1][i] >> 2) + (acf32[2][i] >> 2) + (acf32[3][i] >> 2);

    // smoothing of acf
    for (int i = 0; i < num_coeffs + 1; i++)
        acf32_smooth[i] = (((int64_t) acf32_smooth[i] * acf_alpha0) + ((int64_t) acf32[acf_index][i] * acf_alpha1)) >> 31;

    if (voice_rms)
    {
        error_power32 = durbin32(acf32_smooth, a32_temp, num_coeffs, fractional_digits, k_max);

        // sqrt(error_power32)
        int32_t log_gain = log32(error_power32);
        log_gain = log_gain >> 1;
        error_gain = exp32(log_gain);

        std::unique_lock<std::mutex> locker(a_coeff_mutex, std::defer_lock);
        locker.lock();

        for (int i = 0; i < num_coeffs; i++)
            a32[i] = a32_temp[i];

        locker.unlock();
    }
    else
    {
        error_gain = 0;
    }

    acf_index++;
    if (acf_index >= num_acf)
        acf_index = 0;

    block_ready = false;
}

void TalkBox32::resetStates(void)
{
    voice_rms = 0;
    error_gain = 0;
    buffer_position = 0;
    block_ready = false;

    memory_hp[0] = memory_hp[1] = 0;

    for (int i=0; i<memory_rms_size; i++)
        memory_rms32[i] = 0;

    for (int i=0; i<num_coeffs + 1; i++)
        acf32_smooth[i] = 0;

    for (int i=0; i<num_coeffs; i++)
        memory_lpc[i] = 0;
}

void TalkBox32::setSmoothingTime(float tau)
{
    double alpha;

    if (tau > 0)
        alpha = 1 - (block_length / ( tau * fs ));
    else
        alpha = 0;

    if (alpha < 0)
        alpha = 0;

    acf_alpha0 = (int32_t) (alpha * 0x7FFFFFFF);
    acf_alpha1 = (int32_t) ((1-alpha) * 0x7FFFFFFF);
}

void TalkBox32::setGateLevel(float level)
{
    gate_level = (int32_t) (level * 0x7FFFFFFF);
}

void TalkBox32::setPreemphasis(float fcuttoff)
{
    double ftan = tan(M_PI * fcuttoff / fs);
    high_pass_coeff = (int32_t) ((ftan-1) / (ftan+1) * 0x7FFFFFFF);
}

int TalkBox32::getNumCoeffs(void)
{
    return num_coeffs;
}

void TalkBox32::getCoefficients(float all_pole_coefficients[])
{
    for (int i=0; i<num_coeffs; i++)
        all_pole_coefficients[i] = a32[i] / float(1 << fractional_digits);
}

float TalkBox32::getPreemphasis(void)
{
    return ( high_pass_coeff / float(0x7FFFFFFF) );
}

float TalkBox32::getErrorGain(void)
{
    return ( error_gain / float(0x7FFFFFFF) );
}

float TalkBox32::getVoiceGain(void)
{
    return ( voice_rms / float(0x7FFFFFFF) );
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
