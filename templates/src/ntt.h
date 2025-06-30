#ifndef NTT_H
#define NTT_H

// reduction mode selection: 0=Original, 1=Barrett, 2=HEAX MulRed
// {MR_mode}

// TW_MODE
// 0 = all pre-calculated tables;
// 1 = both tw_factors and tw_h on-the-fly;
// 2 = tw_factors pre-calculated, tw_h on-the-fly
#define TW_MODE 1

#include <cmath>
#include <algorithm>
#include "hls_vector.h"
#include "hls_stream.h"
#include "ap_int.h"
#include <tapa.h>

constexpr int K = {K};

using HostData = {DATA_FORMAT}; 
using Data = ap_uint<K>; // ap_uint<14|23|32>
using Data2 = ap_uint<K*2>;
 
constexpr int DataCHLen = 64 / {DATA_BSIZE}; // 64 / sizeof(Data)
constexpr int EffDataCHLen = {EffDataCHLen}; 

using DataVec = tapa::vec_t<HostData, DataCHLen>;

template <typename T>
using bits = ap_uint<tapa::widthof<T>()>;

constexpr int log2(int x) {
    return (x <= 1) ? 0 : 1 + log2(x / 2);
}

#define MOD {MOD} // 12289 | 8380417 | 3221225473, 167772161 | 469762049 | 998244353 | 2013265921

#if REDUCTION_MODE == 0
// Original modular reduction support is limited 
#if MOD == 12289
    #define USE_Q12289
#elif MOD == 8380417
    #define USE_Q8380417
#elif MOD == 3221225473
    #define USE_Q3221225473
#else
    #error "Unsupported mod value. Please define MOD as 12289 or 8380417."
#endif

#elif REDUCTION_MODE == 1
#if MOD == 3221225473
constexpr unsigned long long BARRETT_MU = 5726623059;
#else
constexpr unsigned long BARRETT_MU = (1ULL << (2 * K)) / MOD;
#endif

#elif REDUCTION_MODE == 2
using Data_W = ap_uint<K+2>;
using Data_W2 = ap_int<K*2+2>;
using Data3  = ap_int<K>;
using Data4  = ap_int<2*K>;

#endif

// Number of coefficients
constexpr int n = {N};
constexpr int logN = {logN};

constexpr int BU = {BU};
constexpr int logBU = {logBU};

// WIDTH: number of coeffs processed in parallel (per NTT CORE)
constexpr int WIDTH = 2*BU;
constexpr int DEPTH = n / WIDTH;
constexpr int logDEPTH = {logDEPTH};

constexpr int num_spat_stage = logBU + 1;
constexpr int num_temp_stage = logN - (logBU + 1);

constexpr int CH = {CH};

constexpr int NUM_CORE = EffDataCHLen*CH / (2*BU);

constexpr int GROUP_NUM = {GROUP_NUM}; //std::min(CH, NUM_CORE);
#define GROUP_CH_NUM {GROUP_CH_NUM}
//constexpr int GROUP_CH_NUM = CH / GROUP_NUM;
constexpr int GROUP_CORE_NUM = NUM_CORE / GROUP_NUM;

#if {MCH} // GROUP_CH_NUM > 1
#define MCH
#endif

constexpr int DRAM_CORE_WIDTH_RATIO = DataCHLen / (2*BU);
constexpr int CORE_DRAM_WIDTH_RATIO = (2*BU) / DataCHLen;

constexpr int POLY_FIFO_DEPTH_S = n / DataCHLen;
constexpr int POLY_FIFO_DEPTH_M = n / (2*BU);

constexpr HostData psi = {PSI};

#if TW_MODE == 0 || TW_MODE == 2
// Bit reversed array of twiddle factors
const Data tw_factors[n] = {{TW_FACTORS}};
#endif

// HEAX MulRed table
#if REDUCTION_MODE == 2
#if TW_MODE == 0
// const Data2 tw_h_factors[MOD] = ;
const Data2 tw_h_factors[n] = {{TW_H_FACTORS}};
#define GET_TW_H(idx) (tw_h_factors[idx])
#endif

#else
  #define GET_TW_H(idx) (Data2(0))
#endif
#endif

#if TW_MODE != 0
inline int bit_reverse(int x, int logN) {
#pragma HLS INLINE
    int y = 0;
    for(int b = 0; b < logN; ++b) {
        y = (y<<1) | (x & 1);
        x >>= 1;
    }
    return y;
}

// base^exp mod MOD
inline Data modexp(Data base, int exp) {
#pragma HLS INLINE
    Data result = 1;
    Data2 b = base;
    while (exp) {
        if (exp & 1) result = (Data)((Data2)result * b % MOD);
        b = (Data2)((Data2)b * b % MOD);
        exp >>= 1;
    }
    return result;
}

// twiddle_generator
inline Data compute_tw_factor(int idx) {
#pragma HLS INLINE
    // int natural_idx = bit_reverse(idx, logN);
    // return modexp(psi, natural_idx);
    return modexp(psi, idx);
}

inline Data2 compute_tw_h(int idx) {
#pragma HLS INLINE
#if TW_MODE == 1
    // All on-the-fly
    Data twf = compute_tw_factor(idx);
    Data2 tw_h_shifted = (Data2)twf << K;
    return (Data2)(tw_h_shifted / MOD);
#elif TW_MODE == 2
    // tw_factors pre-calculated, tw_h on-the-fly
    Data twf = tw_factors[idx]; 
    Data2 tw_h_shifted = (Data2)twf << K;  
    return (Data2)(tw_h_shifted / MOD);  
#else
    return tw_h_factors[idx];
#endif

}

inline Data2 compute_tw_h_from_twf(Data twf) {
#pragma HLS INLINE
    // Directly twf << K
    Data2 shifted = (Data2)twf << K;
    return (Data2)(shifted / MOD);
}

#if REDUCTION_MODE == 2
inline void reduce_mulred(Data coeff, Data tw_factor, Data2 tw_h_factor, Data &remainder){
#pragma HLS INLINE

        Data q = MOD;
        Data2 prod_z = Data2(coeff) * Data2(tw_factor);
        Data_W  z = Data_W(prod_z);
        Data2 prod_t_ = Data2(coeff * tw_h_factor);
        Data_W  t = Data_W(prod_t_ >> K);

        Data2 prod_z2 = Data2(t) * Data2(q);
        Data_W  z2 = Data_W(prod_z2);

        Data_W res = Data4(z) - Data4(z2);
	if (res >= q)
        res = res - q;

	remainder =  static_cast<Data>(res);
}
#endif

#if REDUCTION_MODE == 1
inline void reduce_barrett(Data coeff, Data tw_factor, Data &remainder){
#pragma HLS INLINE

	Data2 q = MOD;

        // Full 2K-bit product
        Data2 x = (Data2)coeff * (Data2)tw_factor;
    
        // BARRETT_MU = floor(2^(2K) / Q)

        // q1 = floor(x / 2^K)
        Data2 q1 = x >> K;
        // q2 = q1 * MU
        Data2 q2 = q1 * (Data2)BARRETT_MU;
        // q3 = floor(q2 / 2^K)
        Data2 q3 = q2 >> K;
        // r = x - q3 * Q
        Data2 r = x - q3 * q;
        // Correct at most twice
        if (r >= q) r -= q;
        if (r >= q) r -= q;

	remainder =  static_cast<Data>(r);
}
#endif

#if REDUCTION_MODE == 0
inline void reduce_original(Data coeff, Data tw_factor, Data &remainder){
#pragma HLS INLINE

	Data q = MOD;
	Data2 odd_cast = coeff;
	Data2 tw_factor_cast = tw_factor;

	Data2 z = odd_cast * tw_factor_cast;

#ifdef USE_Q12289 //Data is 14b

	// max(c) = 19 (maximum 5 bits are needed)
	ap_uint<5> c = z.range(27, 26) + z.range(25, 24) + z.range(23, 22)
		+ z.range(21, 20) + z.range(19, 18) + z.range(17, 16)
		+ z.range(15, 14);

	// max(d) = 14??h2FFF (max 14 bits are needed)
	ap_uint<14> d = z.range(27, 26) + z.range(27, 24) + z.range(27, 22)
		+ z.range(27, 20) + z.range(27, 18) + z.range(27, 16)
		+ z.range(27, 14);

	// max(e) = 6 (max 3 bits are needed)
	ap_uint<3> e = c.range(1, 0) + c.range(3, 2) + c[4];

	ap_uint<16> f = (ap_uint<16>) (e[2] + e.range(1, 0)) << 12;
	f = f -(e[2]+ c[4] + c.range(4, 2));

	ap_uint<14> g = z.range(13, 0);

        Data q_temp = (g > q) ? q : (Data)0;

	// Calculate y
	ap_uint<16> y = f + g - q_temp;

#elif defined(USE_Q8380417) //Data is 23b

	// c = z[45:43} + z[42:33} + z[32:23], max(c) = 2052 (12-bit integer)
	ap_uint<12> c = z.range(45, 43) + z.range(42, 33) + z.range(32, 23);
	// max(d) = 8380415 (23-bit integer)
	ap_uint<23> d = z.range(45, 43) + z.range(45, 33) + z.range(45, 23);

	// max(e) = 1024 (max 11 bits are needed)
	ap_uint<11> e = c.range(11, 10) + c.range(9, 0);


	// f is less than q
	ap_uint<23> f = (ap_uint<23>) (e[10] + e.range(9, 0)) << 13;
	f = f -(e[10]+ c.range(11, 10));

	// g could be bigger than q (23-bit integer) 
	ap_uint<23> g = z.range(22, 0);

        Data q_temp = (g > q) ? q : (Data)0;

	// Calculate y
	ap_uint<32> y = f + g - q_temp;

#elif defined(USE_Q3221225473) //Data is 32b

	// max(c) = 46 (6-bit integer)
	ap_uint<6> c = 	z.range(63, 62) + z.range(61, 60) + z.range(59, 58) + z.range(57, 56)
					+ z.range(55, 54) + z.range(53, 52) + z.range(51, 50) + z.range(49, 48)
					+ z.range(47, 46) + z.range(45, 44) + z.range(43, 42) + z.range(41, 40)
					+ z.range(39, 38) + z.range(37, 36) + z.range(35, 34) + z.range(33, 32);
	
	// max(d) = 8380415 (32-bit integer)
	ap_uint<32> d = z.range(63, 62) + z.range(63, 60) + z.range(63, 58) + z.range(63, 56)
					+ z.range(63, 54) + z.range(63, 52) + z.range(63, 50) + z.range(63, 48)
					+ z.range(63, 46) + z.range(63, 44) + z.range(63, 42) + z.range(63, 40)
					+ z.range(63, 38) + z.range(63, 36) + z.range(63, 34) + z.range(63, 32);

	// max(e) = 7 (max 3 bits are needed)
	ap_uint<3> e = c.range(5, 4) + c.range(3, 2) + c.range(1, 0);

	// max(f) = 4 (3-bit integer)
	ap_uint<3> f = e[2] + e.range(1, 0);

	// max(f[2] + f.range(1, 0)) = 3 -> g is less than q
	ap_uint<32> g = (ap_uint<32>) (f[2] + f.range(1, 0)) << 30 ;
	
	g = g -(f[2] + e[2] + c.range(5, 4) + c.range(5, 2));

	// h could be bigger than q (23-bit integer) 
	ap_uint<32> h = z.range(31, 0);

        Data q_temp = (h > q) ? q : (Data)0;

	// Calculate y
	ap_uint<33> y = g + h - q_temp;

#else
	#error "One of USE_Q12289, USE_Q8380417, or USE_Q3221225473 must be defined."
#endif

        Data q_temp2 = (y > q) ? q : (Data)0;
	y = y - q_temp2;

        Data q_temp3 = (y < d) ? q : (Data)0;
	y = y + q_temp3 - d;

	remainder =  static_cast<Data>(y);
}
#endif

#endif // NTT_H
