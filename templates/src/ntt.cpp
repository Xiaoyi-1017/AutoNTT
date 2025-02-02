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


void mem_stage(const int stage, tapa::istream<Data>& input_stream0, tapa::istream<Data>& input_stream1, 
		tapa::ostream<Data>& output_stream0, tapa::ostream<Data>& output_stream1, int SAMPLES)
{
	const int stride = DEPTH >> (stage+1);
	const int shift = num_temp_stage - (stage+1);
	const int mask = (1<<shift) -1;

	// memory for entry with EVEN indices
	Data mem0[2][DEPTH/2]; 
	Data mem1[2][DEPTH/2];

	// memory for entry with ODD indices
	Data mem2[2][DEPTH/2]; 
	Data mem3[2][DEPTH/2];

	//double buffer index
	ap_uint<1> rd_s = 0;
	ap_uint<1> wr_s = 1;

	//memory read/write data count
	unsigned int read_i = 0; unsigned int write_i = 0;    

	bool read_exist = false;
	bool read_done = false;
	bool write_done = false;

MEM_STAGE_LOOP:
	for(;;){
#pragma HLS pipeline II = 1
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

		if(read_exist == true && write_done == false){
			Data output_data0;
			Data output_data1;

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
				write_i = 0;
				write_done = true;
			}
		}

		if(read_done == false && !input_stream0.empty() && !input_stream1.empty()){
			Data input_data0 = input_stream0.read();
			Data input_data1 = input_stream1.read();

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
				read_i = 0;
				read_done = true;
			}
		}

		//Go to the next sample when all existing data has been written out
		//and when either reading data has finished or no data has been read
		if( (read_done == true || read_i == 0) && (read_exist == false || write_done == true) ){
			if(read_done == true){
				read_exist = true;
			}
			else{
				read_exist = false;
			}
			read_done = false;
			write_done = false;

			wr_s = rd_s;
			rd_s ^= (ap_uint<1>)1;
		}
	}
}

void temporal_stage(int stage, tapa::istreams<Data, B>& input_stream0, tapa::istreams<Data, B>& input_stream1,
		tapa::ostreams<Data, B>& output_stream0, tapa::ostreams<Data, B>& output_stream1, int SAMPLES){

	tapa::streams<Data, B> inter_streams0;
	tapa::streams<Data, B> inter_streams1;

	tapa::task()
		.invoke<tapa::detach, B>(bf_unit, tapa::seq(), stage, input_stream0, input_stream1, inter_streams0, inter_streams1, SAMPLES)
		.invoke<tapa::detach, B>(mem_stage, stage, inter_streams0, inter_streams1, output_stream0, output_stream1, SAMPLES);
}

void spatial_stages(tapa::istreams<Data, B>& input_stream0, tapa::istreams<Data, B>& input_stream1,
		tapa::ostreams<Data, B>& output_stream0, tapa::ostreams<Data, B>& output_stream1, int SAMPLES){

	const int num_of_stages = log2B + 1;
	int mem[num_of_stages+1][WIDTH];
	ap_uint<logDEPTH> i = 0;

NTT_SPATIAL_LOOP:
	for(;;){
#pragma HLS PIPELINE II = 1
		bool do_read = 1;
		for(int j = 0; j < B; j++){
#pragma HLS UNROLL
			do_read &= !input_stream0[j].empty() & !input_stream1[j].empty();
		}

		//for(int i = 0; i < DEPTH; i++){
		if( do_read == true ){
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

			i++;
		}
	}
}

void input_mem_stage_int(int id, tapa::istream<Data>& x_0, tapa::istream<Data>& x_1,
		tapa::ostream<Data>& even_y, tapa::ostream<Data>& odd_y, int SAMPLES){

	// memory for entry with EVEN indices
	Data mem0[2][DEPTH/2]; 
	Data mem1[2][DEPTH/2];

	// memory for entry with ODD indices
	Data mem2[2][DEPTH/2]; 
	Data mem3[2][DEPTH/2];

	//double buffer index
	ap_uint<1> rd_s = 0;
	ap_uint<1> wr_s = 1;

	//memory read/write data count
	unsigned int read_i = 0; unsigned int write_i = 0;    

	bool read_exist = false;
	bool read_done = false;
	bool write_done = false;

INPUT_MEM_STAGE_LOOP:
	for(;;){
#pragma HLS pipeline II = 1
#pragma HLS dependence variable=mem0 type=inter false
#pragma HLS dependence variable=mem1 type=inter false
#pragma HLS dependence variable=mem2 type=inter false
#pragma HLS dependence variable=mem3 type=inter false

		if(read_exist == true && write_done == false){
			Data output_data0;
			Data output_data1;

			// For even inputs 0 ~ n/2-1
			if (write_i % 2 == 0){ // 0, 512 - 2, 514 ...
				output_data0 = mem0[wr_s][write_i/2];
				output_data1 = mem2[wr_s][write_i/2];
			}
			else{            // 1, 513 - 3, 515 ...
				output_data0 = mem1[wr_s][write_i/2];
				output_data1 = mem3[wr_s][write_i/2];
			}

			even_y.write(output_data0);
			odd_y.write(output_data1);

			write_i++;

			if(write_i == DEPTH){
				write_i = 0;
				write_done = true;
			} 
		}

		if(read_done == false && !x_0.empty() && !x_1.empty()){
			Data output_data0 = x_0.read();
			Data output_data1 = x_1.read();

			if(read_i < DEPTH/2){
				// For even inputs 0 ~ n/2-1      
				mem0[rd_s][read_i] = output_data0;
				mem1[rd_s][read_i] = output_data1;
			}
			else{
				// For odd inputs n/2 ~ n-1
				mem2[rd_s][read_i-DEPTH/2] = output_data0;
				mem3[rd_s][read_i-DEPTH/2] = output_data1;
			}

			read_i++;

			if(read_i == DEPTH){
				read_i = 0;
				read_done = true;
			} 
		}

		//Go to the next sample when all existing data has been written out
		//and when either reading data has finished or no data has been read
		if( (read_done == true || read_i == 0) && (read_exist == false || write_done == true) ){
			if(read_done == true){
				read_exist = true;
			}
			else{
				read_exist = false;
			}
			read_done = false;
			write_done = false;

			wr_s = rd_s;
			rd_s ^= (ap_uint<1>)1;
		}
	}
}

void input_mem_stage(tapa::istreams<Data, B>& x_0, tapa::istreams<Data, B>& x_1,
		tapa::ostreams<Data, B>& y_0, tapa::ostreams<Data, B>& y_1, int SAMPLES){

	tapa::task()
		.invoke<tapa::detach, B>(input_mem_stage_int, tapa::seq(), x_0, x_1, y_0, y_1, SAMPLES);
}




#ifdef MCH

void read_dram_m(tapa::mmap<bits<DataVec>> x, tapa::ostreams<Data, 2*B> & dramrd_streams, int sample_num){

	for(int sample = 0; sample < sample_num/NUM_CH; sample++){
        	for(int i = 0; i < n/kVecLen; i++){
#pragma HLS PIPELINE II = 1
			DataVec data = tapa::bit_cast<DataVec>( x[sample * (n/kVecLen) + i] ); 
			for(int d = 0; d < kVecLen; d++){
#pragma HLS UNROLL
				dramrd_streams[(i % GROUP_CH_NUM)*kVecLen+d].write( data[d] );			
			}
		}	
	}
}

void read_collect_2_core_m(tapa::istream<Data> & dramrd_streams0, tapa::istream<Data> & dramrd_streams1, tapa::ostream<Data> & core_istream){

READ_COLLECT_LOOP:
	for(;;){		
		for(int ch = 0; ch < GROUP_CH_NUM; ch++){
			for(int i = 0; i < n/(2*B); i++){
#pragma HLS PIPELINE II = 1
				Data data;
				if( ch == 0 ){
					data = dramrd_streams0.read();
				}
				else{
					data = dramrd_streams1.read();
				}

				core_istream.write( data );
			}
		}
	}
}

void read_collect_4_core_m(tapa::istream<Data> & dramrd_streams0, tapa::istream<Data> & dramrd_streams1, tapa::istream<Data> & dramrd_streams2, tapa::istream<Data> & dramrd_streams3, tapa::ostream<Data> & core_istream){

READ_COLLECT_LOOP:
	for(;;){		
		for(int ch = 0; ch < GROUP_CH_NUM; ch++){
			for(int i = 0; i < n/(2*B); i++){
#pragma HLS PIPELINE II = 1
				Data data;
				if( ch == 0 ){
					data = dramrd_streams0.read();
				}
				else if( ch == 1 ){
					data = dramrd_streams1.read();
				}
				else if( ch == 2 ){
					data = dramrd_streams2.read();
				}
				else{
					data = dramrd_streams3.read();
				}

				core_istream.write( data );
			}
		}
	}
}

void read_collect_8_core_m(tapa::istream<Data> & dramrd_streams0, tapa::istream<Data> & dramrd_streams1, tapa::istream<Data> & dramrd_streams2, tapa::istream<Data> & dramrd_streams3, tapa::istream<Data> & dramrd_streams4, tapa::istream<Data> & dramrd_streams5, tapa::istream<Data> & dramrd_streams6, tapa::istream<Data> & dramrd_streams7, tapa::ostream<Data> & core_istream){

READ_COLLECT_LOOP:
	for(;;){		
		for(int ch = 0; ch < GROUP_CH_NUM; ch++){
			for(int i = 0; i < n/(2*B); i++){
#pragma HLS PIPELINE II = 1
				Data data;
				if( ch == 0 ){
					data = dramrd_streams0.read();
				}
				else if( ch == 1 ){
					data = dramrd_streams1.read();
				}
				else if( ch == 2 ){
					data = dramrd_streams2.read();
				}
				else if( ch == 3 ){
					data = dramrd_streams3.read();
				}
				else if( ch == 4 ){
					data = dramrd_streams4.read();
				}
				else if( ch == 5 ){
					data = dramrd_streams5.read();
				}
				else if( ch == 6 ){
					data = dramrd_streams6.read();
				}
				else{
					data = dramrd_streams7.read();
				}

				core_istream.write( data );
			}
		}
	}
}

void read_collect_2_m(tapa::istreams<Data, 2*B> & dramrd_streams0, tapa::istreams<Data, 2*B> & dramrd_streams1, tapa::ostreams<Data, B> & core_istreams_e, tapa::ostreams<Data, B> & core_istreams_o ){

	tapa::task()
		.invoke<tapa::detach, B>(read_collect_2_core_m, dramrd_streams0, dramrd_streams1, core_istreams_e)
		.invoke<tapa::detach, B>(read_collect_2_core_m, dramrd_streams0, dramrd_streams1, core_istreams_o)
	;
}

void read_collect_4_m(tapa::istreams<Data, 2*B> & dramrd_streams0, tapa::istreams<Data, 2*B> & dramrd_streams1, tapa::istreams<Data, 2*B> & dramrd_streams2, tapa::istreams<Data, 2*B> & dramrd_streams3, tapa::ostreams<Data, B> & core_istreams_e, tapa::ostreams<Data, B> & core_istreams_o ){

	tapa::task()
		.invoke<tapa::detach, B>(read_collect_4_core_m, dramrd_streams0, dramrd_streams1, dramrd_streams2, dramrd_streams3, core_istreams_e)
		.invoke<tapa::detach, B>(read_collect_4_core_m, dramrd_streams0, dramrd_streams1, dramrd_streams2, dramrd_streams3, core_istreams_o)
	;
}

void read_collect_8_m(tapa::istreams<Data, 2*B> & dramrd_streams0, tapa::istreams<Data, 2*B> & dramrd_streams1, tapa::istreams<Data, 2*B> & dramrd_streams2, tapa::istreams<Data, 2*B> & dramrd_streams3, tapa::istreams<Data, 2*B> & dramrd_streams4, tapa::istreams<Data, 2*B> & dramrd_streams5, tapa::istreams<Data, 2*B> & dramrd_streams6, tapa::istreams<Data, 2*B> & dramrd_streams7, tapa::ostreams<Data, B> & core_istreams_e, tapa::ostreams<Data, B> & core_istreams_o ){

	tapa::task()
		.invoke<tapa::detach, B>(read_collect_8_core_m, dramrd_streams0, dramrd_streams1, dramrd_streams2, dramrd_streams3, dramrd_streams4, dramrd_streams5, dramrd_streams6, dramrd_streams7, core_istreams_e)
		.invoke<tapa::detach, B>(read_collect_8_core_m, dramrd_streams0, dramrd_streams1, dramrd_streams2, dramrd_streams3, dramrd_streams4, dramrd_streams5, dramrd_streams6, dramrd_streams7, core_istreams_o)
	;
}

void read_collect_m(tapa::istreams<Data, 2*B*GROUP_CH_NUM> & dramrd_streams, tapa::ostreams<Data, B> & core_istreams_e, tapa::ostreams<Data, B> & core_istreams_o ){

	tapa::task()
#if GROUP_CH_NUM == 2
		.invoke<tapa::detach, 1>(read_collect_2_m, dramrd_streams, dramrd_streams, core_istreams_e, core_istreams_o)
#elif GROUP_CH_NUM == 4
		.invoke<tapa::detach, 1>(read_collect_4_m, dramrd_streams, dramrd_streams, dramrd_streams, dramrd_streams, core_istreams_e, core_istreams_o)
#elif GROUP_CH_NUM == 8
		.invoke<tapa::detach, 1>(read_collect_8_m, dramrd_streams, dramrd_streams, dramrd_streams, dramrd_streams, dramrd_streams, dramrd_streams, dramrd_streams, dramrd_streams, core_istreams_e, core_istreams_o)
#endif
	;
}

void write_dist_2_core_m(tapa::ostream<Data> & dramwr_streams0, tapa::ostream<Data> & dramwr_streams1, tapa::istream<Data> & core_ostream){

WRITE_DIST_LOOP:
	for(;;){		
		for(int ch = 0; ch < GROUP_CH_NUM; ch++){
			for(int i = 0; i < n/(2*B); i++){
#pragma HLS PIPELINE II = 1
				Data data = core_ostream.read();
				if( ch == 0 ){
					dramwr_streams0.write( data );
				}
				else{
					dramwr_streams1.write( data );
				}
			}
		}
	}
}

void write_dist_4_core_m(tapa::ostream<Data> & dramwr_streams0, tapa::ostream<Data> & dramwr_streams1, tapa::ostream<Data> & dramwr_streams2, tapa::ostream<Data> & dramwr_streams3, tapa::istream<Data> & core_ostream){

WRITE_DIST_LOOP:
	for(;;){		
		for(int ch = 0; ch < GROUP_CH_NUM; ch++){
			for(int i = 0; i < n/(2*B); i++){
#pragma HLS PIPELINE II = 1
				Data data = core_ostream.read();
				if( ch == 0 ){
					dramwr_streams0.write( data );
				}
				else if( ch == 1 ){
					dramwr_streams1.write( data );
				}
				else if( ch == 2 ){
					dramwr_streams2.write( data );
				}
				else{
					dramwr_streams3.write( data );
				}
			}
		}
	}
}

void write_dist_8_core_m(tapa::ostream<Data> & dramwr_streams0, tapa::ostream<Data> & dramwr_streams1, tapa::ostream<Data> & dramwr_streams2, tapa::ostream<Data> & dramwr_streams3, tapa::ostream<Data> & dramwr_streams4, tapa::ostream<Data> & dramwr_streams5, tapa::ostream<Data> & dramwr_streams6, tapa::ostream<Data> & dramwr_streams7, tapa::istream<Data> & core_ostream){

WRITE_DIST_LOOP:
	for(;;){		
		for(int ch = 0; ch < GROUP_CH_NUM; ch++){
			for(int i = 0; i < n/(2*B); i++){
#pragma HLS PIPELINE II = 1
				Data data = core_ostream.read();
				if( ch == 0 ){
					dramwr_streams0.write( data );
				}
				else if( ch == 1 ){
					dramwr_streams1.write( data );
				}
				else if( ch == 2 ){
					dramwr_streams2.write( data );
				}
				else if( ch == 3 ){
					dramwr_streams3.write( data );
				}
				else if( ch == 4 ){
					dramwr_streams4.write( data );
				}
				else if( ch == 5 ){
					dramwr_streams5.write( data );
				}
				else if( ch == 6 ){
					dramwr_streams6.write( data );
				}
				else{
					dramwr_streams7.write( data );
				}
			}
		}
	}
}

void write_dist_2_int_m(tapa::ostreams<Data, 2> & dramwr_streams0, tapa::ostreams<Data, 2> & dramwr_streams1, tapa::istream<Data> & core_ostreams_e, tapa::istream<Data> & core_ostreams_o ){

	tapa::task()
		.invoke<tapa::detach>(write_dist_2_core_m, dramwr_streams0, dramwr_streams1, core_ostreams_e)
		.invoke<tapa::detach>(write_dist_2_core_m, dramwr_streams0, dramwr_streams1, core_ostreams_o)
	;
}

void write_dist_2_m(tapa::ostreams<Data, 2*B> & dramwr_streams0, tapa::ostreams<Data, 2*B> & dramwr_streams1, tapa::istreams<Data, B> & core_ostreams_e, tapa::istreams<Data, B> & core_ostreams_o ){

	tapa::task()
		.invoke<tapa::detach, B>(write_dist_2_int_m, dramwr_streams0, dramwr_streams1, core_ostreams_e, core_ostreams_o)
	;
}

void write_dist_4_int_m(tapa::ostreams<Data, 2> & dramwr_streams0, tapa::ostreams<Data, 2> & dramwr_streams1, tapa::ostreams<Data, 2> & dramwr_streams2, tapa::ostreams<Data, 2> & dramwr_streams3, tapa::istream<Data> & core_ostreams_e, tapa::istream<Data> & core_ostreams_o ){

	tapa::task()
		.invoke<tapa::detach>(write_dist_4_core_m, dramwr_streams0, dramwr_streams1, dramwr_streams2, dramwr_streams3, core_ostreams_e)
		.invoke<tapa::detach>(write_dist_4_core_m, dramwr_streams0, dramwr_streams1, dramwr_streams2, dramwr_streams3, core_ostreams_o)
	;
}

void write_dist_4_m(tapa::ostreams<Data, 2*B> & dramwr_streams0, tapa::ostreams<Data, 2*B> & dramwr_streams1, tapa::ostreams<Data, 2*B> & dramwr_streams2, tapa::ostreams<Data, 2*B> & dramwr_streams3, tapa::istreams<Data, B> & core_ostreams_e, tapa::istreams<Data, B> & core_ostreams_o ){

	tapa::task()
		.invoke<tapa::detach, B>(write_dist_4_int_m, dramwr_streams0, dramwr_streams1, dramwr_streams2, dramwr_streams3, core_ostreams_e, core_ostreams_o)
	;
}

void write_dist_8_int_m(tapa::ostreams<Data, 2> & dramwr_streams0, tapa::ostreams<Data, 2> & dramwr_streams1, tapa::ostreams<Data, 2> & dramwr_streams2, tapa::ostreams<Data, 2> & dramwr_streams3, tapa::ostreams<Data, 2> & dramwr_streams4, tapa::ostreams<Data, 2> & dramwr_streams5, tapa::ostreams<Data, 2> & dramwr_streams6, tapa::ostreams<Data, 2> & dramwr_streams7, tapa::istream<Data> & core_ostreams_e, tapa::istream<Data> & core_ostreams_o ){

	tapa::task()
		.invoke<tapa::detach>(write_dist_8_core_m, dramwr_streams0, dramwr_streams1, dramwr_streams2, dramwr_streams3, dramwr_streams4, dramwr_streams5, dramwr_streams6, dramwr_streams7, core_ostreams_e)
		.invoke<tapa::detach>(write_dist_8_core_m, dramwr_streams0, dramwr_streams1, dramwr_streams2, dramwr_streams3, dramwr_streams4, dramwr_streams5, dramwr_streams6, dramwr_streams7, core_ostreams_o)
	;
}

void write_dist_8_m(tapa::ostreams<Data, 2*B> & dramwr_streams0, tapa::ostreams<Data, 2*B> & dramwr_streams1, tapa::ostreams<Data, 2*B> & dramwr_streams2, tapa::ostreams<Data, 2*B> & dramwr_streams3, tapa::ostreams<Data, 2*B> & dramwr_streams4, tapa::ostreams<Data, 2*B> & dramwr_streams5, tapa::ostreams<Data, 2*B> & dramwr_streams6, tapa::ostreams<Data, 2*B> & dramwr_streams7, tapa::istreams<Data, B> & core_ostreams_e, tapa::istreams<Data, B> & core_ostreams_o ){

	tapa::task()
		.invoke<tapa::detach, B>(write_dist_8_int_m, dramwr_streams0, dramwr_streams1, dramwr_streams2, dramwr_streams3, dramwr_streams4, dramwr_streams5, dramwr_streams6, dramwr_streams7, core_ostreams_e, core_ostreams_o)
	;
}


void write_dist_m(tapa::ostreams<Data, 2*B*GROUP_CH_NUM> & dramwr_streams, tapa::istreams<Data, B> & core_ostreams_e, tapa::istreams<Data, B> & core_ostreams_o ){

	tapa::task()
#if GROUP_CH_NUM == 2
		.invoke<tapa::detach, 1>(write_dist_2_m, dramwr_streams, dramwr_streams, core_ostreams_e, core_ostreams_o)
#elif GROUP_CH_NUM == 4
		.invoke<tapa::detach, 1>(write_dist_4_m, dramwr_streams, dramwr_streams, dramwr_streams, dramwr_streams, core_ostreams_e, core_ostreams_o)
#elif GROUP_CH_NUM == 8
		.invoke<tapa::detach, 1>(write_dist_8_m, dramwr_streams, dramwr_streams, dramwr_streams, dramwr_streams, dramwr_streams, dramwr_streams, dramwr_streams, dramwr_streams, core_ostreams_e, core_ostreams_o)
#endif
	;
}

void write_dram_m(tapa::mmap<bits<DataVec>> y, tapa::istreams<Data, 2*B> & dramwr_streams, int sample_num){

	for(int sample = 0; sample < sample_num/NUM_CH; sample++){
        	for(int i = 0; i < n/kVecLen; i++){
#pragma HLS PIPELINE II = 1
			DataVec data;
			for(int d = 0; d < kVecLen; d++){
#pragma HLS UNROLL
				data[d] = dramwr_streams[(i % GROUP_CH_NUM)*kVecLen+d].read();			
			}
            		y[sample * (n/kVecLen) + i] =  tapa::bit_cast<bits<DataVec>>(data);
		}	
	}
}

#else

void read_dram_s(tapa::mmap<bits<DataVec>> x, tapa::ostreams<DataVec, GROUP_CORE_NUM> & dramrd_streams, int sample_num){

	for(int sample = 0; sample < sample_num/NUM_CH; sample++){
	       	for(int i = 0; i < n/kVecLen; i++){
#pragma HLS PIPELINE II = 1
			DataVec data = tapa::bit_cast<DataVec>( x[sample * (n/kVecLen) + i] );
			dramrd_streams[sample%GROUP_CORE_NUM].write(data);
		}
	}
}

void read_dist_s(tapa::istream<DataVec> & dramrd_stream, tapa::ostreams<Data, B> & core_istreams_e, tapa::ostreams<Data, B> & core_istreams_o){
	for(;;){
	        //for(int i = 0; i < n/kVecLen; i++){
#pragma HLS PIPELINE II = GROUP_CORE_NUM
		if( !dramrd_stream.empty() ){
			DataVec data = dramrd_stream.read();
			for(int core = 0; core < GROUP_CORE_NUM; core++){
				for(int d = 0; d < B; d++){
#pragma HLS UNROLL
					core_istreams_e[d].write( data[core*2*B+d] );			
				}
				for(int d = 0; d < B; d++){
#pragma HLS UNROLL
					core_istreams_o[d].write( data[core*2*B+B+d] );			
				}
			}
		}
		//}
	}
}

void write_reshape_s(tapa::ostream<DataVec> & dramwr_stream, tapa::istreams<Data, B> & core_ostreams_e, tapa::istreams<Data, B> & core_ostreams_o ){

WRITE_RESHAPE_S_LOOP:
	for(;;){
#pragma HLS PIPELINE II = GROUP_CORE_NUM 
		DataVec data;
		for(int c = 0; c < GROUP_CORE_NUM; c++){
			for(int d = 0; d < B; d++){
				data[c*2*B+2*d] = core_ostreams_e[d].read();
			}
			for(int d = 0; d < B; d++){
				data[c*2*B+2*d+1] = core_ostreams_o[d].read();
			}
		}
		dramwr_stream.write(data);
	}
}

void write_dram_s(tapa::mmap<bits<DataVec>> y, tapa::istreams<DataVec, GROUP_CORE_NUM> & dramwr_streams, int sample_num){

	for(int sample = 0; sample < sample_num/NUM_CH; sample++){
	       	for(int i = 0; i < n/kVecLen; i++){
#pragma HLS PIPELINE II = 1
			DataVec data = dramwr_streams[sample%GROUP_CORE_NUM].read();
            		y[sample * (n/kVecLen) + i] =  tapa::bit_cast<bits<DataVec>>(data);
		}
	}
}

#endif


void ntt_core(int id, tapa::istreams<Data, B> core_istreams_e, tapa::istreams<Data, B> core_istreams_o, tapa::ostreams<Data, B> core_ostreams_e, tapa::ostreams<Data, B> core_ostreams_o, int SAMPLES){

	tapa::streams<Data, B*(num_temp_stage+1), 2> even_streams("even_streams");
	tapa::streams<Data, B*(num_temp_stage+1), 2> odd_streams("odd_streams");

	tapa::task()
		.invoke<tapa::detach>(input_mem_stage, core_istreams_e, core_istreams_o, even_streams, odd_streams, SAMPLES)
		.invoke<tapa::detach, num_temp_stage>(temporal_stage, tapa::seq(), even_streams, odd_streams, even_streams, odd_streams, SAMPLES)
		.invoke<tapa::detach>(spatial_stages, even_streams, odd_streams, core_ostreams_e, core_ostreams_o, SAMPLES); 
}


void ntt(tapa::mmaps<bits<DataVec>, NUM_CH> x, tapa::mmaps<bits<DataVec>, NUM_CH> y, int sample_num){

#ifdef MCH
	tapa::streams<Data, 2*B*NUM_CH, SAMPLE_FIFO_DEPTH_M> dramrd_streams("dramrd_streams");
	tapa::streams<Data, 2*B*NUM_CH, SAMPLE_FIFO_DEPTH_M> dramwr_streams("dramwr_streams");
#else
	tapa::streams<DataVec, NUM_CORE, SAMPLE_FIFO_DEPTH_S> dramrd_streams("dramrd_streams");
	tapa::streams<DataVec, NUM_CORE, SAMPLE_FIFO_DEPTH_S> dramwr_streams("dramwr_streams");
#endif
	tapa::streams<Data, B*NUM_CORE, 2> core_istreams_e("core_istreams_e");
	tapa::streams<Data, B*NUM_CORE, 2> core_istreams_o("core_istreams_o");

	tapa::streams<Data, B*NUM_CORE, 2> core_ostreams_e("core_ostreams_e");
	tapa::streams<Data, B*NUM_CORE, 2> core_ostreams_o("core_ostreams_o");

	tapa::task()
#ifdef MCH
		.invoke<tapa::join, NUM_CH>(read_dram_m, x, dramrd_streams, sample_num)
		.invoke<tapa::detach, GROUP_NUM>(read_collect_m, dramrd_streams, core_istreams_e, core_istreams_o)
#else
		.invoke<tapa::join, NUM_CH>(read_dram_s, x, dramrd_streams, sample_num)
	 	.invoke<tapa::detach, NUM_CORE>(read_dist_s, dramrd_streams, core_istreams_e, core_istreams_o)
#endif

		.invoke<tapa::detach, NUM_CORE>(ntt_core, tapa::seq(), core_istreams_e, core_istreams_o, core_ostreams_e, core_ostreams_o, sample_num)

#ifdef MCH
		.invoke<tapa::detach, GROUP_NUM>(write_dist_m, dramwr_streams, core_ostreams_e, core_ostreams_o)
		.invoke<tapa::join, NUM_CH>(write_dram_m, y, dramwr_streams, sample_num)
#else
		.invoke<tapa::detach, NUM_CORE>(write_reshape_s, dramwr_streams, core_ostreams_e, core_ostreams_o)
		.invoke<tapa::join, NUM_CH>(write_dram_s, y, dramwr_streams, sample_num)
#endif
	;

}
