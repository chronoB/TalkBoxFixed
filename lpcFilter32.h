#ifndef _LPCFILTER32
#define _LPCFILTER32

#include <stdint.h>

int32_t lpcFilter32(int32_t inputSample, int32_t *a, int32_t *memory, int num_coeff, const int fractional_digits);

#endif  // _LPCFILTER32
