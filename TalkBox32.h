#ifndef _TALK_BOX32
#define _TALK_BOX32

#include <stdint.h>
#include <mutex>

const int num_coeffs = 50;
const int block_length = 512;
const int num_acf = 4;
const int memory_rms_size = 4;
const int fractional_digits = 24;

class TalkBox32
{
protected:
    double fs;
    int32_t voice_rms;
    int32_t error_gain;
    int32_t buffer_position;
    int32_t input_buffer0[block_length];
    int32_t input_buffer1[block_length];
    int32_t *sample_buffer;
    int32_t *block_buffer;
    bool block_ready;
    int16_t n_shift_memory;
    int16_t n_shift_block;
    int32_t high_pass_coeff;
    int32_t memory_hp[2];
    int32_t memory_rms32[memory_rms_size];
    int32_t acf_alpha0;
    int32_t acf_alpha1;
    int32_t gate_level;
    int16_t acf_index;
    int32_t acf32[num_acf][num_coeffs + 1];
    int32_t acf32_smooth[num_coeffs + 1];
    int32_t a32_temp[num_coeffs];
    int32_t a32[num_coeffs];
    int32_t memory_lpc[num_coeffs];
    std::mutex a_coeff_mutex;

public:
    TalkBox32(double fs);
    ~TalkBox32(void);
    void process(int32_t samples[]);
    void calculateLPCcoefficients(void);
    void resetStates(void);
    void setSmoothingTime(float tau);
    void setGateLevel(float level);
    void setPreemphasis(float fcuttoff);
    int  getNumCoeffs(void);
    void getCoefficients(float all_pole_coefficients[]);
    float getPreemphasis(void);
    float getErrorGain(void);
    float getVoiceGain(void);
};

#endif  // _TALK_BOX32
