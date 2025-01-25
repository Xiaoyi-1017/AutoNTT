#ifndef NTT_H
#define NTT_H
#include <cmath>
#include "hls_vector.h"
#include "hls_stream.h"
#include "ap_int.h"
#include <tapa.h>

using Data = {DATA_FORMAT};
constexpr int kVecLen = 64 / sizeof(Data);

using DataVec = tapa::vec_t<Data, kVecLen>;

template <typename T>
using bits = ap_uint<tapa::widthof<T>()>;

const int mod = {MOD};

// Number of coefficients
constexpr int n = {N};
constexpr int log2N = {log2N};

const int B = {B};
const int log2B = {log2B};

// WIDTH: number of coeffs processed in parallel (per NTT CORE)
const int WIDTH = 2*B;
const int DEPTH = n / WIDTH;
const int logDEPTH = {logDEPTH};

const int spatial_stages = log2B + 1;
const int temporal_stages = log2N - (log2B + 1);

const int NUM_CH = {NUM_CH};

const int NUM_CORE = kVecLen*NUM_CH / (2*B);

const int NUM_CH_PER_CORE = (2*B) / kVecLen; 
const int NUM_CORE_PER_CH = kVecLen / (2*B);


// Bit reversed array of twiddle factors
const int tw_factors[n] = {{TW_FACTORS}};

const int psi = tw_factors[n/2];

#endif // NTT_H
