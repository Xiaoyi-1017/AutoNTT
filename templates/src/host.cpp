#include <algorithm>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>
#include "ntt.h"


int bit_reverse(int i){
	ap_uint<log2N> x = (ap_uint<log2N>) i;
	ap_uint<log2N> reversed_idx = x.reverse();
	return (int) reversed_idx;
}

int mod_power(int x, int exp, int mod){
    int result = 1;
    for(int i = 0; i < exp; i++){
        result *= x;
        result %= mod;
    }
    return result;

}

void get_omega_mat(int Omega[n][n], int mod, int psi) {
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            int exp = (2*(i*j)+j) % (2*n);
            Omega[i][j] = mod_power(psi, exp, mod);
        }
    }
}


void sw_ntt(std::vector<Data, tapa::aligned_allocator<Data>> A, std::vector<Data, tapa::aligned_allocator<Data>> &out_sw, int psi, int q, const int SAMPLE_NUM){
  int omega[n][n];
  get_omega_mat(omega, mod, psi);
  #pragma omp parallel
	{
		int tid = omp_get_thread_num();
		if( tid == 0 ){
			int nthreads = omp_get_num_threads();
			// std::cout << "Running OpenMP with " << nthreads << " threads...\n";
		}
	}
  #pragma omp parallel for
  for(int i = 0; i < SAMPLE_NUM; i++){
    for (int j = 0; j < n; ++j) {
        int sum = 0;
        for (int k = 0; k < n; ++k) {
            sum = (sum + omega[j][k] * A[n*i+k]) % q;
        }
        out_sw[n*i + j] = sum;
    }
  }
}

void bit_reverse_hw_out(std::vector<Data, tapa::aligned_allocator<Data>> out_hw, std::vector<Data, tapa::aligned_allocator<Data>> &out_hw_BR, const int SAMPLE_NUM){
#pragma omp parallel 
{
    int tid = omp_get_thread_num();
    if( tid == 0 ){
            int nthreads = omp_get_num_threads();
            // std::cout << "Running OpenMP with " << nthreads << " threads...\n";
    }
}
#pragma omp parallel for
//Bit-reversing the HW output
    for(int i = 0; i < SAMPLE_NUM; i++){
      for(int j = 0; j < n; j++){
          out_hw_BR[n*i + bit_reverse(j)] = out_hw[n*i + j];
      }
    }
}

void copy_input(std::vector<Data, tapa::aligned_allocator<Data>> input, 
		std::vector<std::vector<int, tapa::aligned_allocator<int> >, std::allocator<std::vector<int, tapa::aligned_allocator<int> > >> &X, int SAMPLE_NUM){

	for(int i = 0; i < SAMPLE_NUM; i++){
		for(int j = 0; j < n; j++){
			//X[i%CH][(i/CH)*n+j] = rand() % mod;
			X[i%CH][(i/CH)*n+j] = input[i*n + j];
		}
	}
}

void copy_output(std::vector<std::vector<int, tapa::aligned_allocator<int> >, std::allocator<std::vector<int, tapa::aligned_allocator<int> > >> Y, 
		std::vector<Data, tapa::aligned_allocator<Data>> &out_hw, int SAMPLE_NUM){
	for(int i = 0; i < SAMPLE_NUM; i++){
		for(int j = 0; j < n; j++){
			out_hw[i*n + j] = Y[i%CH][(i/CH)*n+j];
		}
	}
}




void ntt(tapa::mmaps<bits<DataVec>, CH> x, tapa::mmaps<bits<DataVec>, CH> y, int SAMPLE_NUM);

DEFINE_string(bitstream, "", "path to bitstream file, run csim if empty");

int main(int argc, char* argv[]) {

    gflags::ParseCommandLineFlags(&argc, &argv, /*remove_flags=*/true);

    const int SAMPLE_NUM = argc > 1 ? atoll(argv[1]) : 1000;

    const int N = n;

    const int ADJ_FACTOR = CH > NUM_CORE ? CH : NUM_CORE;
    const int ADJ_SAMPLE_NUM = ((SAMPLE_NUM - 1) / ADJ_FACTOR + 1) * ADJ_FACTOR;

    std::vector<Data, tapa::aligned_allocator<Data>> input(n*ADJ_SAMPLE_NUM);
    
    std::vector<Data, tapa::aligned_allocator<Data>> out_sw(n*ADJ_SAMPLE_NUM);
    std::vector<Data, tapa::aligned_allocator<Data>> out_hw(n*ADJ_SAMPLE_NUM);
    std::vector<Data, tapa::aligned_allocator<Data>> out_hw_BR(n*ADJ_SAMPLE_NUM);

    std::vector<std::vector<Data, tapa::aligned_allocator<Data>>> X(CH);
    std::vector<std::vector<Data, tapa::aligned_allocator<Data>>> Y(CH);



    int last_stride = n >> (log2N - log2BU -1);

    // Generate twiddle factors
    std::cout << "n: " << n << " mod: " << mod << std::endl;
    std::cout << "Number OF samples: " << ADJ_SAMPLE_NUM << " (adjusted from " << SAMPLE_NUM << ")" << std::endl; 
    std::cout << "Number of NTT cores: " << NUM_CORE << std::endl;
    std::cout << "Number of butterfly units per NTT core: " << BU << std::endl; 
    std::cout << "Number of input and output DRAM channels: " << CH << std::endl;
    std::cout << "DEPTH per sample: " << DEPTH << std::endl;
  
    srand(time(NULL));

    // Create the test data 
    for(int i = 0; i < ADJ_SAMPLE_NUM; i++){
      for(int j = 0; j < n; j++){
        //   input[i*n + j] = rand() % mod;
          input[i*n + j] = (i*n + j) % mod;
      }
    }

    // Resize each channel into appropriate size
    int CHANNEL_SIZE = n*ADJ_SAMPLE_NUM/CH;

    if( CHANNEL_SIZE * sizeof(Data) > 256 * 1024 * 1024 ){
        std::cout << "Error: Data size of each channel (" << CHANNEL_SIZE * sizeof(Data) << ") exceeds 256 MB" << std::endl;
        exit(1);
    }  

    for (int cc = 0; cc < CH; ++cc) {
        X[cc].resize(CHANNEL_SIZE, 0);
        Y[cc].resize(CHANNEL_SIZE, 0);
    }

    copy_input(input, X, ADJ_SAMPLE_NUM);
    
    std::cout << "Test data generation DONE" << std::endl;

    sw_ntt(input, out_sw, psi, mod, ADJ_SAMPLE_NUM);

    std::cout << "SW Computation Done. Launching FPGA Kernel..." << std::endl;

    int64_t kernel_time_ns =
      tapa::invoke(ntt, FLAGS_bitstream,
                   tapa::read_only_mmaps<Data, CH>(X).reinterpret<bits<DataVec>>(),
                   tapa::write_only_mmaps<Data, CH>(Y).reinterpret<bits<DataVec>>(), 
		ADJ_SAMPLE_NUM);


    //std::cout << "kernel time: " << kernel_time_ns * 1e-9 << " s" << std::endl;
    std::cout << "Kernel time: " << (double)kernel_time_ns * 1e-6  << " ms" << std::endl;
    //std::cout << "kernel time: " << kernel_time_ns * 1e-3 << " us" << std::endl;

    std::cout << "Kernel throughput: " << (double)ADJ_SAMPLE_NUM / kernel_time_ns * 1e3 << " Msamples/s" << std::endl; 
    std::cout << "Eff BW per HBM channel: " << (double)CHANNEL_SIZE*sizeof(Data) / kernel_time_ns << " GB/s" << std::endl; 
    
    std::cout << "Copying output ...\n";
  
    copy_output(Y, out_hw, ADJ_SAMPLE_NUM);    

    bit_reverse_hw_out(out_hw, out_hw_BR, ADJ_SAMPLE_NUM);

    std::cout << "Done\n";
        

    // Compare the results of the Device to the simulation

	int err_cnt = 0;
	for(int i = 0; i<ADJ_SAMPLE_NUM; i++){
		for(int j = 0; j< n; j++){
			//printf("Sample %d index %d - sw:%d, hw:%d\n", i, j, out_sw[i*n+j], out_hw_BR[i*n+j]); 
			if(out_sw[i*n+j] != out_hw_BR[i*n+j]) {
				err_cnt++;
				if(err_cnt < 10) printf("Error in sample %d index %d - sw:%d, hw:%d\n", i, j, out_sw[i*n+j], out_hw_BR[i*n+j]); 
			}
		}
	}

	if(err_cnt != 0){
		printf("FAILED! Error count : %d / %d\n", err_cnt, ADJ_SAMPLE_NUM*n);
	}
	else{
		printf("PASSED!\n");
	}

	return EXIT_SUCCESS;
}
