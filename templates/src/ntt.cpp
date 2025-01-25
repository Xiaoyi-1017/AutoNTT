#include <iostream>
#include "ntt.h"

void reduce(int x, int &remainder){
#pragma HLS INLINE
	ap_uint <16> q = 12289;
	ap_uint <32> z = static_cast<ap_uint<32>>(x);

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

	// Calculate y
	ap_uint<32> y = f + g - q*(g > q);

	y = y- q*(y > q);

	y = y + q * (y < d) - d;

	remainder =  static_cast<int>(y);
}

void butterfly(int even, int odd, int tw_factor, int *out_even, int* out_odd){
#pragma HLS INLINE
	int multiplied = odd * tw_factor;
	int t; reduce(multiplied, t);

	// Cooley-Tukey butterfly
	int coeff_odd = even - t;
	int coeff_even = even + t;

	// Save them to dereferenced variables out_odd, and out_even
	*out_odd = (coeff_odd < 0) ? (coeff_odd + mod) : coeff_odd;
	*out_even = (coeff_even >= mod) ? (coeff_even - mod) : coeff_even;
}

void bf_unit(int id, int stage, tapa::istream<int>& input0, tapa::istream<int>& input1,
		tapa::ostream<int>& output0, tapa::ostream<int>& output1, int SAMPLES){

	int stride = n >> (stage + 1); 

	int even, odd;
	int out_even, out_odd;
	ap_uint<logDEPTH> read_i = 0;

BF_LOOP:
	for(;;){
#pragma HLS PIPELINE II = 1
		if(!input0.empty() && !input1.empty()){

			even = input0.read();
			odd = input1.read();

			int tw_idx = (read_i*B)/stride +  n/(2 * stride);

			butterfly(even, odd,  tw_factors[tw_idx], &out_even, &out_odd);

			output0.write(out_even);
			output1.write(out_odd);                

			read_i++;
		}
	}
}


void mem_stage(const int stage, tapa::istream<int>& input_stream0, tapa::istream<int>& input_stream1, 
		tapa::ostream<int>& output_stream0, tapa::ostream<int>& output_stream1, int SAMPLES)
{
	const int stride = DEPTH >> (stage+1);
	const int shift = temporal_stages - (stage+1);
	const int mask = (1<<shift) -1;

	// memory for entry with EVEN indices
	int mem0[2][DEPTH/2]; 
	int mem1[2][DEPTH/2];

	// memory for entry with ODD indices
	int mem2[2][DEPTH/2]; 
	int mem3[2][DEPTH/2];

	int input_data0; int input_data1;    
	int output_data0; int output_data1;

	bool input_suc0 = false; bool input_suc1 = false;
	bool output_suc0 = false; bool output_suc1 = false;

	bool read_exist = 0;
	ap_uint<1> rd_s = 0;

MEM_STAGE_LOOP:
	for(;;){
		bool read_done = 0;
		bool write_done = !read_exist; // write can proceed only when read data exists
		read_exist = 0;

		//switch read/write sample
		ap_uint<1> wr_s = rd_s;
		rd_s ^= (ap_uint<1>)1;

		//memory read/write data count
		unsigned int read_i = 0; unsigned int write_i = 0;    

		bool sample_done = false;

MEM_STAGE_SAMPLE_LOOP:
		while(!sample_done){
#pragma HLS pipeline II = 1;
#pragma HLS dependence variable=mem0 type=inter false
#pragma HLS dependence variable=mem1 type=inter false
#pragma HLS dependence variable=mem2 type=inter false
#pragma HLS dependence variable=mem3 type=inter false

			// Indexing
			int read_mem_idx = (read_i >> shift) % 2;
			int read_chunk = (read_i >> shift) / 2;
			int read_idx = (read_i & mask);

			int raddr = read_chunk*stride + read_idx;

			int write_mem_idx = (write_i >> shift) % 2;
			int write_chunk = (write_i >> shift) / 2;
			int write_idx = (write_i & mask);

			int waddr = (write_chunk*stride) +  write_idx;

			if(read_done == 0 && !input_stream0.empty() && !input_stream1.empty()){
				input_data0 = input_stream0.read();
				input_data1 = input_stream1.read();

				if(read_mem_idx ==  0){    
					mem0[rd_s][raddr] = input_data0;
					mem1[rd_s][raddr] = input_data1;
				}
				else{
					mem2[rd_s][raddr] = input_data0;
					mem3[rd_s][raddr] = input_data1;
				}

				read_i++;

				if(read_i == DEPTH){
					read_done = true;
					read_exist = true;
				}
			}

			if(write_done == 0){
				if (write_mem_idx == 0){
					output_data0 = mem0[wr_s][waddr];
					output_data1 = mem2[wr_s][waddr];
				}
				else{
					output_data0 = mem1[wr_s][waddr];
					output_data1 = mem3[wr_s][waddr];
				}

				output_stream0.write(output_data0);
				output_stream1.write(output_data1);

				write_i++;

				if(write_i == DEPTH){
					write_done = true;
				}
			}

			//Go to the next sample when all existing data has been written out
			//and when either reading data has finished or no data has been read
			if( (read_done == true || read_i == 0) && write_done == true){
				sample_done = true;
			}
		}
	}
}

void ntt_temporal_stage(int stage, tapa::istreams<int, B>& input_stream0, tapa::istreams<int, B>& input_stream1,
		tapa::ostreams<int, B>& output_stream0, tapa::ostreams<int, B>& output_stream1, int SAMPLES){

	tapa::streams<int, B> inter_streams0;
	tapa::streams<int, B> inter_streams1;

	tapa::task()
		.invoke<tapa::detach, B>(bf_unit, tapa::seq(), stage, input_stream0, input_stream1, inter_streams0, inter_streams1, SAMPLES)
		.invoke<tapa::detach, B>(mem_stage, stage, inter_streams0, inter_streams1, output_stream0, output_stream1, SAMPLES);
}

void ntt_spatial_stages(tapa::istreams<int, B>& input_stream0, tapa::istreams<int, B>& input_stream1,
		tapa::ostreams<int, B>& output_stream0, tapa::ostreams<int, B>& output_stream1, int SAMPLES){

	const int num_of_stages = log2B + 1;
	int mem[num_of_stages+1][WIDTH];

	for(;;){
NTT_SPATIAL_LOOP:
		for(int i = 0; i < DEPTH; i++){
#pragma HLS PIPELINE

INPUT_LOOP:
			for(int j = 0; j < B; j++){
#pragma HLS UNROLL
				mem[0][2*j] = input_stream0[j].read();
				mem[0][2*j+1] = input_stream1[j].read();
			}

STAGE_LOOP:
			for(int s = 0; s < num_of_stages; s++){
#pragma HLS unroll 

				int current_stage = s + log2N - log2B -1;
				int stride = n >> ((current_stage) + 1);
				int next_stride = stride >> 1;

				// For division
				int shift = num_of_stages - (s + 1);

				// For remainder
				int mask = (1 << shift) -1;
				int next_mask = ( 1 << (shift-1) ) -1;

BUTTERFLY_LOOP:
				for(int idx = 0; idx < B; idx++){
#pragma HLS unroll
					int out_even; int out_odd;

					int j = idx >> shift;
					int k = idx & mask;

					int tw_idx = (i*B+idx)/stride + n / (2 * stride);

					int ind_even = (next_stride == 0) ? ( j << (shift+1) ) 
						: ( j << (shift+1) ) + (k & next_mask) * 2 + (k >> (shift-1));
					int ind_odd = ind_even + stride;

					butterfly(mem[s][2*idx], mem[s][2*idx+1], tw_factors[tw_idx], &out_even, &out_odd);
					mem[s+1][ind_even] = out_even;
					mem[s+1][ind_odd] = out_odd;                   
				}
			}

OUTPUT_LOOP:
			for(int j = 0; j < B; j++){
#pragma HLS UNROLL
				output_stream0[j] << mem[num_of_stages][2*j];
				output_stream1[j] << mem[num_of_stages][2*j+1];
			}
		}
	}
}

void input_module(int id, tapa::mmap<bits<DataVec>> x_in, tapa::ostreams<int, kVecLen>& input_stream0, int SAMPLES){

    for(int sample = 0; sample < SAMPLES; sample++){
        for(int i = 0; i < DEPTH; i++){
        #pragma HLS PIPELINE II = 1
            DataVec in = tapa::bit_cast<DataVec>(x_in[sample*DEPTH  + i]);
            
            for(int k = 0; k < kVecLen; k++){
            #pragma HLS UNROLL    
                input_stream0[k] << in[k];
            }
        }
    }
}

void output_module(int id, tapa::istreams<int, kVecLen/2>& output_stream0, tapa::istreams<int, kVecLen/2>& output_stream1, 
                    tapa::mmap<bits<DataVec>> y_out, int SAMPLES){
    
    for(int sample = 0; sample < SAMPLES; sample++){
        for(int i = 0; i < DEPTH; i++){
        #pragma HLS PIPELINE II = 1
            DataVec out;
            for(int k = 0; k < (kVecLen/2); k++){
            #pragma HLS UNROLL    
                out[2*k] = output_stream0[k].read();
                out[2*k+1] = output_stream1[k].read();
            }
            y_out[sample * DEPTH + i] =  tapa::bit_cast<bits<DataVec>>(out);
        }
    }
}

void input_mem_stage_int(int id, tapa::istream<int>& x_0, tapa::istream<int>& x_1,
		tapa::ostream<int>& even_y, tapa::ostream<int>& odd_y, int SAMPLES){

	// memory for entry with EVEN indices
	int mem0[2][DEPTH/2]; 
	int mem1[2][DEPTH/2];

	// memory for entry with ODD indices
	int mem2[2][DEPTH/2]; 
	int mem3[2][DEPTH/2];

	bool read_exist = 0;
	ap_uint<1> rd_s = 0;

INPUT_MEM_STAGE_LOOP:
	for(;;){
		bool read_done = 0;
		bool write_done = !read_exist; // write can proceed only when read data exists
		read_exist = 0;

		//switch read/write sample
		ap_uint<1> wr_s = rd_s;
		rd_s ^= (ap_uint<1>)1;

		//memory read/write data count
		unsigned int read_i = 0; unsigned int write_i = 0;    

		bool sample_done = false;

INPUT_MEM_STAGE_SAMPLE_LOOP:
		while(!sample_done){
#pragma HLS pipeline II = 1;
#pragma HLS dependence variable=mem0 type=inter false
#pragma HLS dependence variable=mem1 type=inter false
#pragma HLS dependence variable=mem2 type=inter false
#pragma HLS dependence variable=mem3 type=inter false

			if(read_done == false && !x_0.empty() && !x_1.empty()){
				if(read_i < DEPTH/2){
					// For even inputs 0 ~ n/2-1      
					mem0[rd_s][read_i] = x_0.read();
					mem1[rd_s][read_i] = x_1.read();
				}
				else{
					// For odd inputs n/2 ~ n-1
					mem2[rd_s][read_i-DEPTH/2] = x_0.read();
					mem3[rd_s][read_i-DEPTH/2] = x_1.read();
				}

				read_i++;

				if(read_i == DEPTH){
					read_done = true;
					read_exist = true;
				} 
			}

			if(write_done == false){
				// For even inputs 0 ~ n/2-1
				if (write_i % 2 == 0){ // 0, 512 - 2, 514 ...
					even_y.write(mem0[wr_s][write_i/2]);
					odd_y.write(mem2[wr_s][write_i/2]);
				}
				else{            // 1, 513 - 3, 515 ...
					even_y.write(mem1[wr_s][write_i/2]);
					odd_y.write(mem3[wr_s][write_i/2]);
				}

				write_i++;

				if(write_i == DEPTH){
					write_done = true;
				} 
			}

			//Go to the next sample when all existing data has been written out
			//and when either reading data has finished or no data has been read
			if( (read_done == true || read_i == 0) && write_done == true){
				sample_done = true;
			}
		}
	}
}

void input_mem_stage(tapa::istreams<int, B>& x_0, tapa::istreams<int, B>& x_1,
		tapa::ostreams<int, B>& y_0, tapa::ostreams<int, B>& y_1, int SAMPLES){

	tapa::task()
		.invoke<tapa::detach, B>(input_mem_stage_int, tapa::seq(), x_0, x_1, y_0, y_1, SAMPLES);
}

void divide_streams(tapa::istreams<int, 2*B> input_streams, tapa::ostreams<int, B> outputstreams_0, tapa::ostreams<int, B> outputstreams_1, int SAMPLES){

	for(;;){
		for(int i = 0; i < DEPTH; i++){
#pragma HLS pipeline II = 1;

			for(int k = 0; k < B; k++){
#pragma HLS UNROLL    
				outputstreams_0[k] << input_streams[k].read();
			}

			for(int k = 0; k < B; k++){
#pragma HLS UNROLL    
				outputstreams_1[k] << input_streams[k+B].read();
			}
		}
	}
}

void ntt_core(int id, tapa::istreams<int, 2*B> inputstreams, tapa::ostreams<int, B> outputstreams_0, tapa::ostreams<int, B> outputstreams_1, int SAMPLES){

	tapa::streams<int, B*(temporal_stages+2)> even_streams("even_streams");
	tapa::streams<int, B*(temporal_stages+2)> odd_streams("odd_streams");

	tapa::task()
		.invoke<tapa::detach>(divide_streams, inputstreams, even_streams, odd_streams, SAMPLES)
		.invoke<tapa::detach>(input_mem_stage, even_streams, odd_streams, even_streams, odd_streams, SAMPLES)
		.invoke<tapa::detach, temporal_stages>(ntt_temporal_stage, tapa::seq(), even_streams, odd_streams, even_streams, odd_streams, SAMPLES)
		.invoke<tapa::detach>(ntt_spatial_stages, even_streams, odd_streams, outputstreams_0, outputstreams_1, SAMPLES); 
}

void ntt(tapa::mmaps<bits<DataVec>, NUM_CH> x, tapa::mmaps<bits<DataVec>, NUM_CH> y, int SAMPLES){

	tapa::streams<int, NUM_CH*kVecLen> inputstreams("inputstreams");
	tapa::streams<int, NUM_CH*kVecLen/2> outputstreams_0("outputstreams_0");
	tapa::streams<int, NUM_CH*kVecLen/2> outputstreams_1("outputstreams_1");

	tapa::task()
		.invoke<tapa::join, NUM_CH>(input_module, tapa::seq(), x, inputstreams, SAMPLES)
		.invoke<tapa::detach, NUM_CORE>(ntt_core, tapa::seq(), inputstreams, outputstreams_0, outputstreams_1, SAMPLES)
		.invoke<tapa::join, NUM_CH>(output_module, tapa::seq(), outputstreams_0, outputstreams_1, y, SAMPLES);
}

