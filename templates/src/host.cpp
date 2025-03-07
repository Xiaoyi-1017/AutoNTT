#include <algorithm>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>
#include "ntt.h"

int bit_reverse(int i){
	ap_uint<logN> x = (ap_uint<logN>) i;
	ap_uint<logN> reversed_idx = x.reverse();
	return (int) reversed_idx;
}

/* Changed to prevent overflow */
int mod_power(HostData x, int exp, HostData mod){
	uint64_t result = 1;
	for(int i = 0; i < exp; i++){
		result *= x;
		result %= mod;
	}

	HostData result_output = static_cast<HostData>(result);

	return result_output;
}

void get_omega_mat(HostData Omega[n][n], HostData mod, HostData psi) {
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			int exp = (2*(i*j)+j) % (2*n);
			Omega[i][j] = mod_power(psi, exp, mod);
		}
	}
}

/* Changed to prevent overflow */

void sw_ntt(std::vector<HostData, tapa::aligned_allocator<HostData>> A, std::vector<HostData, tapa::aligned_allocator<HostData>> &out_sw, 
		HostData psi, HostData q, const int POLY_NUM){

	HostData omega[n][n];
	get_omega_mat(omega, q, psi);

#pragma omp parallel
	{
		int tid = omp_get_thread_num();
		if( tid == 0 ){
			int nthreads = omp_get_num_threads();
		}
	}

#pragma omp parallel for
	for(int i = 0; i < POLY_NUM; i++){
		for (int j = 0; j < n; ++j) {
			uint64_t sum = 0;
			for (int k = 0; k < n; ++k) {
				sum = sum + static_cast<uint64_t>(omega[j][k]) * static_cast<uint64_t>(A[n*i+k]); 
				sum %= q;
			}
			out_sw[n*i + j] = static_cast<HostData>(sum);
		}
	}
}

void bit_reverse_hw_out(std::vector<HostData, tapa::aligned_allocator<HostData>> out_hw, std::vector<HostData, tapa::aligned_allocator<HostData>> &out_hw_BR, const int POLY_NUM){

#pragma omp parallel 
	{
		int tid = omp_get_thread_num();
		if( tid == 0 ){
			int nthreads = omp_get_num_threads();
		}
	}

#pragma omp parallel for
	//Bit-reversing the HW output
	for(int i = 0; i < POLY_NUM; i++){
		for(int j = 0; j < n; j++){
			out_hw_BR[n*i + bit_reverse(j)] = out_hw[n*i + j];
		}
	}
}

void copy_input(std::vector<HostData, tapa::aligned_allocator<HostData>> input, 
		std::vector<std::vector<HostData, tapa::aligned_allocator<HostData> >, std::allocator<std::vector<HostData, tapa::aligned_allocator<HostData> > >> &HBM_CH, int POLY_NUM){

	for(int i = 0; i < POLY_NUM; i++){
		int index = (((i/GROUP_CH_NUM)*(2*GROUP_CH_NUM))%(2*CH)) + (i%GROUP_CH_NUM); 
		for(int j = 0; j < n; j++){
			//X[index][(i/CH)*n+j] = rand() % mod;
			HBM_CH[index][(i/CH)*n+j] = input[i*n + j];
		}
	}
}

void copy_output(std::vector<std::vector<HostData, tapa::aligned_allocator<HostData> >, std::allocator<std::vector<HostData, tapa::aligned_allocator<HostData> > >> HBM_CH, 
		std::vector<HostData, tapa::aligned_allocator<HostData>> &out_hw, int POLY_NUM){

	for(int i = 0; i < POLY_NUM; i++){
		int index = (((i/GROUP_CH_NUM)*(2*GROUP_CH_NUM))%(2*CH)) + (i%GROUP_CH_NUM) + GROUP_CH_NUM; 
		for(int j = 0; j < n; j++){
			out_hw[i*n + j] = HBM_CH[index][(i/CH)*n+j];
		}
	}
}


void ntt(tapa::mmaps<bits<DataVec>, 2*CH> hbm_ch, int poly_num);

DEFINE_string(bitstream, "", "path to bitstream file, run csim if empty");

int main(int argc, char* argv[]) {

	gflags::ParseCommandLineFlags(&argc, &argv, /*remove_flags=*/true);

	const int POLY_NUM = argc > 1 ? atoll(argv[1]) : 1000;

	const HostData mod = MOD;

	const int N = n;

	const int ADJ_FACTOR = CH > NUM_CORE ? CH : NUM_CORE;
	const int ADJ_POLY_NUM = ((POLY_NUM - 1) / ADJ_FACTOR + 1) * ADJ_FACTOR;

	std::vector<HostData, tapa::aligned_allocator<HostData>> input(n*ADJ_POLY_NUM);

	std::vector<HostData, tapa::aligned_allocator<HostData>> out_sw(n*ADJ_POLY_NUM);
	std::vector<HostData, tapa::aligned_allocator<HostData>> out_hw(n*ADJ_POLY_NUM);
	std::vector<HostData, tapa::aligned_allocator<HostData>> out_hw_BR(n*ADJ_POLY_NUM);

	std::vector<std::vector<HostData, tapa::aligned_allocator<HostData>>> HBM_CH(2*CH);

	int last_stride = n >> (logN - logBU -1);

	// Generate twiddle factors
	std::cout << "n: " << n << " mod: " << mod << std::endl;
	std::cout << "Number of polynomials: " << ADJ_POLY_NUM << " (adjusted from " << POLY_NUM << ")" << std::endl; 
	std::cout << "Number of NTT cores: " << NUM_CORE << std::endl;
	std::cout << "Number of butterfly units per NTT core per stage: " << BU << std::endl; 
	std::cout << "Number of total butterfly units: " << BU*NUM_CORE*logN << std::endl; 
	std::cout << "Number of input and output DRAM channels: " << CH << std::endl;
	std::cout << "DEPTH per polynomial: " << DEPTH << std::endl;

	srand(time(NULL));

	// Create the test HostData 
	for(int i = 0; i < ADJ_POLY_NUM; i++){
		for(int j = 0; j < n; j++){
			//   input[i*n + j] = rand() % mod;
			input[i*n + j] = (i*n + j) % mod;
		}
	}

	// Resize each channel into appropriate size
	int CHANNEL_SIZE = n*ADJ_POLY_NUM/CH;

	if( CHANNEL_SIZE * sizeof(HostData) > 256 * 1024 * 1024 ){
		std::cout << "Error: HostData size of each channel (" << CHANNEL_SIZE * sizeof(HostData) << ") exceeds 256 MB" << std::endl;
		exit(1);
	}  

	for (int cc = 0; cc < 2*CH; ++cc) {
		HBM_CH[cc].resize(CHANNEL_SIZE, 0);
	}

	copy_input(input, HBM_CH, ADJ_POLY_NUM);

	std::cout << "Test HostData generation DONE" << std::endl;

	sw_ntt(input, out_sw, psi, mod, ADJ_POLY_NUM);

	std::cout << "SW Computation Done. Launching FPGA Kernel..." << std::endl;

	int64_t kernel_time_ns =
		tapa::invoke(ntt, FLAGS_bitstream,
				tapa::read_write_mmaps<HostData, 2*CH>(HBM_CH).reinterpret<bits<DataVec>>(),
				ADJ_POLY_NUM);


	//std::cout << "kernel time: " << kernel_time_ns * 1e-9 << " s" << std::endl;
	std::cout << "Kernel time: " << (double)kernel_time_ns * 1e-6  << " ms" << std::endl;
	//std::cout << "kernel time: " << kernel_time_ns * 1e-3 << " us" << std::endl;

	std::cout << "Kernel throughput: " << (double)ADJ_POLY_NUM / kernel_time_ns * 1e3 << " M Polynomials/s" << std::endl; 
	std::cout << "Eff BW per HBM channel: " << (double)CHANNEL_SIZE*sizeof(HostData) / kernel_time_ns << " GB/s" << std::endl; 

	std::cout << "Copying output ...\n";

	copy_output(HBM_CH, out_hw, ADJ_POLY_NUM);    

	bit_reverse_hw_out(out_hw, out_hw_BR, ADJ_POLY_NUM);

	std::cout << "Done\n";

	// Compare the results of the Device to the simulation

	int err_cnt = 0;
	for(int i = 0; i<ADJ_POLY_NUM; i++){
		for(int j = 0; j< n; j++){
			//printf("Sample %d index %d - sw:%d, hw:%d\n", i, j, out_sw[i*n+j], out_hw_BR[i*n+j]); 
			if(out_sw[i*n+j] != out_hw_BR[i*n+j]) {
				err_cnt++;
				if(err_cnt < 10) printf("Error in polynomial %d index %d - sw:%d, hw:%d\n", i, j, out_sw[i*n+j], out_hw_BR[i*n+j]); 
			}
		}
	}

	if(err_cnt != 0){
		printf("FAILED! Error count : %d / %d\n", err_cnt, ADJ_POLY_NUM*n);
	}
	else{
		printf("PASSED!\n");
	}

	return EXIT_SUCCESS;
}
