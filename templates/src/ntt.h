#ifndef NTT_H
#define NTT_H
#include <cmath>
#include <algorithm>
#include "hls_vector.h"
#include "hls_stream.h"
#include "ap_int.h"
#include <tapa.h>

using Data = {DATA_FORMAT};
constexpr int kVecLen = 64 / {DATA_BSIZE}; // / sizeof(Data)

using DataVec = tapa::vec_t<Data, kVecLen>;

template <typename T>
using bits = ap_uint<tapa::widthof<T>()>;

constexpr int log2(int x) {
    return (x <= 1) ? 0 : 1 + log2(x / 2);
}

constexpr int mod = {MOD};

// Number of coefficients
constexpr int n = {N};
constexpr int log2N = {log2N};

constexpr int B = {B};
constexpr int log2B = {log2B};

// WIDTH: number of coeffs processed in parallel (per NTT CORE)
constexpr int WIDTH = 2*B;
constexpr int DEPTH = n / WIDTH;
constexpr int logDEPTH = {logDEPTH};

constexpr int num_spat_stage = log2B + 1;
constexpr int num_temp_stage = log2N - (log2B + 1);

constexpr int NUM_CH = {NUM_CH};

constexpr int NUM_CORE = kVecLen*NUM_CH / (2*B);

constexpr int NUM_CH_PER_CORE = (2*B) / kVecLen; 
constexpr int NUM_CORE_PER_CH = kVecLen / (2*B);

#if {MCH} //(2*B) > kVecLen
#define MCH
#endif

constexpr int GROUP_NUM = {GROUP_NUM}; //std::min(NUM_CH, NUM_CH * kVecLen / (2*B));
#define GROUP_CH_NUM {GROUP_CH_NUM}
//constexpr int GROUP_CH_NUM = NUM_CH / GROUP_NUM;
constexpr int log_GROUP_CH_NUM = log2(GROUP_CH_NUM);
constexpr int GROUP_CORE_NUM = kVecLen*NUM_CH / (2*B*GROUP_NUM);
constexpr int log_GROUP_CORE_NUM = log2(GROUP_CORE_NUM);
using ReshapeDataVec = tapa::vec_t<Data, 2*B>;
constexpr int SAMPLE_FIFO_DEPTH_S = n / kVecLen;
constexpr int SAMPLE_FIFO_DEPTH_M = n / (2*B);

// Bit reversed array of twiddle factors
constexpr int tw_factors[n] = {{TW_FACTORS}};

constexpr int psi = tw_factors[n/2];

#endif // NTT_H
