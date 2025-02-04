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

constexpr int BU = {BU};
constexpr int log2BU = {log2BU};

// WIDTH: number of coeffs processed in parallel (per NTT CORE)
constexpr int WIDTH = 2*BU;
constexpr int DEPTH = n / WIDTH;
constexpr int logDEPTH = {logDEPTH};

constexpr int num_spat_stage = log2BU + 1;
constexpr int num_temp_stage = log2N - (log2BU + 1);

constexpr int CH = {CH};

constexpr int NUM_CORE = kVecLen*CH / (2*BU);

constexpr int NUM_CH_PER_CORE = (2*BU) / kVecLen; 
constexpr int NUM_CORE_PER_CH = kVecLen / (2*BU);

#if {MCH} //(2*BU) > kVecLen
#define MCH
#endif

constexpr int GROUP_NUM = {GROUP_NUM}; //std::min(CH, CH * kVecLen / (2*BU));
#define GROUP_CH_NUM {GROUP_CH_NUM}
//constexpr int GROUP_CH_NUM = CH / GROUP_NUM;
constexpr int log_GROUP_CH_NUM = log2(GROUP_CH_NUM);
constexpr int GROUP_CORE_NUM = kVecLen*CH / (2*BU*GROUP_NUM);
constexpr int log_GROUP_CORE_NUM = log2(GROUP_CORE_NUM);
using ReshapeDataVec = tapa::vec_t<Data, 2*BU>;
constexpr int SAMPLE_FIFO_DEPTH_S = n / kVecLen;
constexpr int SAMPLE_FIFO_DEPTH_M = n / (2*BU);

// Bit reversed array of twiddle factors
constexpr int tw_factors[n] = {{TW_FACTORS}};

constexpr int psi = tw_factors[n/2];

#endif // NTT_H
