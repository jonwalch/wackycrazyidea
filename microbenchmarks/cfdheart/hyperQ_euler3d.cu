// Copyright 2009, Andrew Corrigan, acorriga@gmu.edu
// This code is from the AIAA-2009-4001 paper

//#include <cutil.h>
#include <helper_cuda.h>
#include <helper_timer.h>
#include <iostream>
#include <fstream>


//heartwall
#include <avilib.h>
#include <avimod.h>
#include "define.c"
params_common_change common_change;
__constant__ params_common_change d_common_change;
params_common common;
__constant__ params_common d_common;
params_unique unique[ALL_POINTS];								// cannot determine size dynamically so choose more than usually needed
__constant__ params_unique d_unique[ALL_POINTS];

 
/*
 * Options 
 * 
 */ 
#define GAMMA 1.4f
#define iterations 2000
#ifndef block_length
	#define block_length 256
#endif

#define NDIM 3
#define NNB 4

#define RK 3	// 3rd order RK
#define ff_mach 1.2f
#define deg_angle_of_attack 0.0f



#define VAR_DENSITY 0
#define VAR_MOMENTUM  1
#define VAR_DENSITY_ENERGY (VAR_MOMENTUM+NDIM)
#define NVAR (VAR_DENSITY_ENERGY+1)


/*
 * Generic functions
 */
template <typename T>
T* alloc(int N)
{
	T* t;
	checkCudaErrors(cudaMalloc((void**)&t, sizeof(T)*N));
	return t;
}

template <typename T>
void dealloc(T* array)
{
	checkCudaErrors(cudaFree((void*)array));
}

template <typename T>
void copy(T* dst, T* src, int N)
{
	checkCudaErrors(cudaMemcpy((void*)dst, (void*)src, N*sizeof(T), cudaMemcpyDeviceToDevice));
}

template <typename T>
void upload(T* dst, T* src, int N)
{
	checkCudaErrors(cudaMemcpy((void*)dst, (void*)src, N*sizeof(T), cudaMemcpyHostToDevice));
}

/*
 * Element-based Cell-centered FVM solver functions
 */
__constant__ float ff_variable[NVAR];
__constant__ float3 ff_flux_contribution_momentum_x[1];
__constant__ float3 ff_flux_contribution_momentum_y[1];
__constant__ float3 ff_flux_contribution_momentum_z[1];
__constant__ float3 ff_flux_contribution_density_energy[1];

__global__ void cuda_initialize_variables(int nelr, float* variables)
{
	const int i = (blockDim.x*blockIdx.x + threadIdx.x);
	for(int j = 0; j < NVAR; j++)
		variables[i + j*nelr] = ff_variable[j];
}
void initialize_variables(int nelr, float* variables)
{
	dim3 Dg(nelr / block_length), Db(block_length);
	cuda_initialize_variables<<<Dg, Db>>>(nelr, variables);
	getLastCudaError("initialize_variables failed");
}

__device__ __host__ inline void compute_flux_contribution(float& density, float3& momentum, float& density_energy, float& pressure, float3& velocity, float3& fc_momentum_x, float3& fc_momentum_y, float3& fc_momentum_z, float3& fc_density_energy)
{
	fc_momentum_x.x = velocity.x*momentum.x + pressure;
	fc_momentum_x.y = velocity.x*momentum.y;
	fc_momentum_x.z = velocity.x*momentum.z;
	
	
	fc_momentum_y.x = fc_momentum_x.y;
	fc_momentum_y.y = velocity.y*momentum.y + pressure;
	fc_momentum_y.z = velocity.y*momentum.z;

	fc_momentum_z.x = fc_momentum_x.z;
	fc_momentum_z.y = fc_momentum_y.z;
	fc_momentum_z.z = velocity.z*momentum.z + pressure;

	float de_p = density_energy+pressure;
	fc_density_energy.x = velocity.x*de_p;
	fc_density_energy.y = velocity.y*de_p;
	fc_density_energy.z = velocity.z*de_p;
}

__device__ inline void compute_velocity(float& density, float3& momentum, float3& velocity)
{
	velocity.x = momentum.x / density;
	velocity.y = momentum.y / density;
	velocity.z = momentum.z / density;
}
	
__device__ inline float compute_speed_sqd(float3& velocity)
{
	return velocity.x*velocity.x + velocity.y*velocity.y + velocity.z*velocity.z;
}

__device__ inline float compute_pressure(float& density, float& density_energy, float& speed_sqd)
{
	return (float(GAMMA)-float(1.0f))*(density_energy - float(0.5f)*density*speed_sqd);
}

__device__ inline float compute_speed_of_sound(float& density, float& pressure)
{
	return sqrtf(float(GAMMA)*pressure/density);
}


//===============================================================================================================================================================================================================
//===============================================================================================================================================================================================================
//	KERNEL FUNCTION
//===============================================================================================================================================================================================================
//===============================================================================================================================================================================================================

__device__ void heartwall_kernel() {

	//======================================================================================================================================================
	//	COMMON VARIABLES
	//======================================================================================================================================================

	/*__shared__ volatile int smem[1024];*/

	fp* d_in;
	int rot_row;
	int rot_col;
	int in2_rowlow;
	int in2_collow;
	int ic;
	int jc;
	int jp1;
	int ja1, ja2;
	int ip1;
	int ia1, ia2;
	int ja, jb;
	int ia, ib;
	float s;
	int i;
	int j;
	int row;
	int col;
	int ori_row;
	int ori_col;
	int position;
	float sum;
	int pos_ori;
	float temp;
	float temp2;
	int location;
	int cent;
	int tMask_row; 
	int tMask_col;
	float largest_value_current = 0;
	float largest_value = 0;
	int largest_coordinate_current = 0;
	int largest_coordinate = 0;
	float fin_max_val = 0;
	int fin_max_coo = 0;
	int largest_row;
	int largest_col;
	int offset_row;
	int offset_col;
	__shared__ float in_partial_sum[51];															// WATCH THIS !!! HARDCODED VALUE
	__shared__ float in_sqr_partial_sum[51];															// WATCH THIS !!! HARDCODED VALUE
	__shared__ float in_final_sum;
	__shared__ float in_sqr_final_sum;
	float mean;
	float mean_sqr;
	float variance;
	float deviation;
	__shared__ float denomT;
	__shared__ float par_max_val[131];															// WATCH THIS !!! HARDCODED VALUE
	__shared__ int par_max_coo[131];															// WATCH THIS !!! HARDCODED VALUE
	int pointer;
	__shared__ float d_in_mod_temp[2601];
	int ori_pointer;
	int loc_pointer;

	//======================================================================================================================================================
	//	THREAD PARAMETERS
	//======================================================================================================================================================

	int bx = blockIdx.x;																// get current horizontal block index (0-n)
	int tx = threadIdx.x;																// get current horizontal thread index (0-n)
	int ei_new;

	/*smem[threadIdx.x] = 0;*/

	//===============================================================================================================================================================================================================
	//===============================================================================================================================================================================================================
	//	GENERATE TEMPLATE
	//===============================================================================================================================================================================================================
	//===============================================================================================================================================================================================================

	// generate templates based on the first frame only
	if(d_common_change.frame_no == 0){

		//======================================================================================================================================================
		// GET POINTER TO TEMPLATE FOR THE POINT
		//======================================================================================================================================================

		// pointers to: current template for current point
		d_in = &d_unique[bx].d_T[d_unique[bx].in_pointer];

		//======================================================================================================================================================
		//	UPDATE ROW LOC AND COL LOC
		//======================================================================================================================================================

		// uptade temporary endo/epi row/col coordinates (in each block corresponding to point, narrow work to one thread)
		ei_new = tx;
		if(ei_new == 0){

			// update temporary row/col coordinates
			pointer = d_unique[bx].point_no*d_common.no_frames+d_common_change.frame_no;
			d_unique[bx].d_tRowLoc[pointer] = d_unique[bx].d_Row[d_unique[bx].point_no];
			d_unique[bx].d_tColLoc[pointer] = d_unique[bx].d_Col[d_unique[bx].point_no];

		}

		//======================================================================================================================================================
		//	CREATE TEMPLATES
		//======================================================================================================================================================

		// work
		ei_new = tx;
		while(ei_new < d_common.in_elem){

			// figure out row/col location in new matrix
			row = (ei_new+1) % d_common.in_rows - 1;												// (0-n) row
			col = (ei_new+1) / d_common.in_rows + 1 - 1;											// (0-n) column
			if((ei_new+1) % d_common.in_rows == 0){
				row = d_common.in_rows - 1;
				col = col-1;
			}

			// figure out row/col location in corresponding new template area in image and give to every thread (get top left corner and progress down and right)
			ori_row = d_unique[bx].d_Row[d_unique[bx].point_no] - 25 + row - 1;
			ori_col = d_unique[bx].d_Col[d_unique[bx].point_no] - 25 + col - 1;
			ori_pointer = ori_col*d_common.frame_rows+ori_row;

			// update template
			d_in[col*d_common.in_rows+row] = d_common_change.d_frame[ori_pointer];

			// go for second round
			ei_new = ei_new + NUMBER_THREADS;

		}

	}

	//===============================================================================================================================================================================================================
	//===============================================================================================================================================================================================================
	//	PROCESS POINTS
	//===============================================================================================================================================================================================================
	//===============================================================================================================================================================================================================

	// process points in all frames except for the first one
	if(d_common_change.frame_no != 0){

		//======================================================================================================================================================
		//	SELECTION
		//======================================================================================================================================================

		in2_rowlow = d_unique[bx].d_Row[d_unique[bx].point_no] - d_common.sSize;													// (1 to n+1)
		in2_collow = d_unique[bx].d_Col[d_unique[bx].point_no] - d_common.sSize;

		// work
		ei_new = tx;
		while(ei_new < d_common.in2_elem){

			// figure out row/col location in new matrix
			row = (ei_new+1) % d_common.in2_rows - 1;												// (0-n) row
			col = (ei_new+1) / d_common.in2_rows + 1 - 1;											// (0-n) column
			if((ei_new+1) % d_common.in2_rows == 0){
				row = d_common.in2_rows - 1;
				col = col-1;
			}

			// figure out corresponding location in old matrix and copy values to new matrix
			ori_row = row + in2_rowlow - 1;
			ori_col = col + in2_collow - 1;
			d_unique[bx].d_in2[ei_new] = d_common_change.d_frame[ori_col*d_common.frame_rows+ori_row];

			// go for second round
			ei_new = ei_new + NUMBER_THREADS;

		}

		//======================================================================================================================================================
		//	SYNCHRONIZE THREADS
		//======================================================================================================================================================

		__syncthreads();

		//======================================================================================================================================================
		//	CONVOLUTION
		//======================================================================================================================================================

		//====================================================================================================
		//	ROTATION
		//====================================================================================================

		// variables
		d_in = &d_unique[bx].d_T[d_unique[bx].in_pointer];

		// work
		ei_new = tx;
		while(ei_new < d_common.in_elem){

			// figure out row/col location in padded array
			row = (ei_new+1) % d_common.in_rows - 1;												// (0-n) row
			col = (ei_new+1) / d_common.in_rows + 1 - 1;											// (0-n) column
			if((ei_new+1) % d_common.in_rows == 0){
				row = d_common.in_rows - 1;
				col = col-1;
			}
		
			// execution
			rot_row = (d_common.in_rows-1) - row;
			rot_col = (d_common.in_rows-1) - col;
			d_in_mod_temp[ei_new] = d_in[rot_col*d_common.in_rows+rot_row];

			// go for second round
			ei_new = ei_new + NUMBER_THREADS;

		}

		//====================================================================================================
		//	SYNCHRONIZE THREADS
		//====================================================================================================

		__syncthreads();

		//====================================================================================================
		//	ACTUAL CONVOLUTION
		//====================================================================================================

		// work
		ei_new = tx;
		while(ei_new < d_common.conv_elem){

			// figure out row/col location in array
			ic = (ei_new+1) % d_common.conv_rows;												// (1-n)
			jc = (ei_new+1) / d_common.conv_rows + 1;											// (1-n)
			if((ei_new+1) % d_common.conv_rows == 0){
				ic = d_common.conv_rows;
				jc = jc-1;
			}

			//
			j = jc + d_common.joffset;
			jp1 = j + 1;
			if(d_common.in2_cols < jp1){
				ja1 = jp1 - d_common.in2_cols;
			}
			else{
				ja1 = 1;
			}
			if(d_common.in_cols < j){
				ja2 = d_common.in_cols;
			}
			else{
				ja2 = j;
			}

			i = ic + d_common.ioffset;
			ip1 = i + 1;
			
			if(d_common.in2_rows < ip1){
				ia1 = ip1 - d_common.in2_rows;
			}
			else{
				ia1 = 1;
			}
			if(d_common.in_rows < i){
				ia2 = d_common.in_rows;
			}
			else{
				ia2 = i;
			}

			s = 0;

			for(ja=ja1; ja<=ja2; ja++){
				jb = jp1 - ja;
				for(ia=ia1; ia<=ia2; ia++){
					ib = ip1 - ia;
					s = s + d_in_mod_temp[d_common.in_rows*(ja-1)+ia-1] * d_unique[bx].d_in2[d_common.in2_rows*(jb-1)+ib-1];
				}
			}

			//d_unique[bx].d_conv[d_common.conv_rows*(jc-1)+ic-1] = s;
			d_unique[bx].d_conv[ei_new] = s;

			// go for second round
			ei_new = ei_new + NUMBER_THREADS;

		}

		//======================================================================================================================================================
		//	SYNCHRONIZE THREADS
		//======================================================================================================================================================

		__syncthreads();

		//======================================================================================================================================================
		//	CUMULATIVE SUM
		//======================================================================================================================================================

		//====================================================================================================
		//	PAD ARRAY, VERTICAL CUMULATIVE SUM
		//====================================================================================================

		//==================================================
		//	PADD ARRAY
		//==================================================

		// work
		ei_new = tx;
		while(ei_new < d_common.in2_pad_cumv_elem){

			// figure out row/col location in padded array
			row = (ei_new+1) % d_common.in2_pad_cumv_rows - 1;												// (0-n) row
			col = (ei_new+1) / d_common.in2_pad_cumv_rows + 1 - 1;											// (0-n) column
			if((ei_new+1) % d_common.in2_pad_cumv_rows == 0){
				row = d_common.in2_pad_cumv_rows - 1;
				col = col-1;
			}

			// execution
			if(	row > (d_common.in2_pad_add_rows-1) &&														// do if has numbers in original array
				row < (d_common.in2_pad_add_rows+d_common.in2_rows) && 
				col > (d_common.in2_pad_add_cols-1) && 
				col < (d_common.in2_pad_add_cols+d_common.in2_cols)){
				ori_row = row - d_common.in2_pad_add_rows;
				ori_col = col - d_common.in2_pad_add_cols;
				d_unique[bx].d_in2_pad_cumv[ei_new] = d_unique[bx].d_in2[ori_col*d_common.in2_rows+ori_row];
			}
			else{																			// do if otherwise
				d_unique[bx].d_in2_pad_cumv[ei_new] = 0;
			}

			// go for second round
			ei_new = ei_new + NUMBER_THREADS;

		}

		//==================================================
		//	SYNCHRONIZE THREADS
		//==================================================

		__syncthreads();

		//==================================================
		//	VERTICAL CUMULATIVE SUM
		//==================================================

		//work
		ei_new = tx;
		while(ei_new < d_common.in2_pad_cumv_cols){

			// figure out column position
			pos_ori = ei_new*d_common.in2_pad_cumv_rows;

			// variables
			sum = 0;
			
			// loop through all rows
			for(position = pos_ori; position < pos_ori+d_common.in2_pad_cumv_rows; position = position + 1){
				d_unique[bx].d_in2_pad_cumv[position] = d_unique[bx].d_in2_pad_cumv[position] + sum;
				sum = d_unique[bx].d_in2_pad_cumv[position];
			}

			// go for second round
			ei_new = ei_new + NUMBER_THREADS;

		}

		//====================================================================================================
		//	SYNCHRONIZE THREADS
		//====================================================================================================

		__syncthreads();

		//====================================================================================================
		//	SELECTION
		//====================================================================================================

		// work
		ei_new = tx;
		while(ei_new < d_common.in2_pad_cumv_sel_elem){

			// figure out row/col location in new matrix
			row = (ei_new+1) % d_common.in2_pad_cumv_sel_rows - 1;												// (0-n) row
			col = (ei_new+1) / d_common.in2_pad_cumv_sel_rows + 1 - 1;											// (0-n) column
			if((ei_new+1) % d_common.in2_pad_cumv_sel_rows == 0){
				row = d_common.in2_pad_cumv_sel_rows - 1;
				col = col-1;
			}

			// figure out corresponding location in old matrix and copy values to new matrix
			ori_row = row + d_common.in2_pad_cumv_sel_rowlow - 1;
			ori_col = col + d_common.in2_pad_cumv_sel_collow - 1;
			d_unique[bx].d_in2_pad_cumv_sel[ei_new] = d_unique[bx].d_in2_pad_cumv[ori_col*d_common.in2_pad_cumv_rows+ori_row];

			// go for second round
			ei_new = ei_new + NUMBER_THREADS;

		}

		//====================================================================================================
		//	SYNCHRONIZE THREADS
		//====================================================================================================

		__syncthreads();

		//====================================================================================================
		//	SELECTION 2, SUBTRACTION, HORIZONTAL CUMULATIVE SUM
		//====================================================================================================

		//==================================================
		//	SELECTION 2
		//==================================================

		// work
		ei_new = tx;
		while(ei_new < d_common.in2_sub_cumh_elem){

			// figure out row/col location in new matrix
			row = (ei_new+1) % d_common.in2_sub_cumh_rows - 1;												// (0-n) row
			col = (ei_new+1) / d_common.in2_sub_cumh_rows + 1 - 1;											// (0-n) column
			if((ei_new+1) % d_common.in2_sub_cumh_rows == 0){
				row = d_common.in2_sub_cumh_rows - 1;
				col = col-1;
			}

			// figure out corresponding location in old matrix and copy values to new matrix
			ori_row = row + d_common.in2_pad_cumv_sel2_rowlow - 1;
			ori_col = col + d_common.in2_pad_cumv_sel2_collow - 1;
			d_unique[bx].d_in2_sub_cumh[ei_new] = d_unique[bx].d_in2_pad_cumv[ori_col*d_common.in2_pad_cumv_rows+ori_row];

			// go for second round
			ei_new = ei_new + NUMBER_THREADS;

		}

		//==================================================
		//	SYNCHRONIZE THREADS
		//==================================================

		__syncthreads();

		//==================================================
		//	SUBTRACTION
		//==================================================
		
		// work
		ei_new = tx;
		while(ei_new < d_common.in2_sub_cumh_elem){

			// subtract
			d_unique[bx].d_in2_sub_cumh[ei_new] = d_unique[bx].d_in2_pad_cumv_sel[ei_new] - d_unique[bx].d_in2_sub_cumh[ei_new];

			// go for second round
			ei_new = ei_new + NUMBER_THREADS;

		}

		//==================================================
		//	SYNCHRONIZE THREADS
		//==================================================

		__syncthreads();

		//==================================================
		//	HORIZONTAL CUMULATIVE SUM
		//==================================================

		// work
		ei_new = tx;
		while(ei_new < d_common.in2_sub_cumh_rows){

			// figure out row position
			pos_ori = ei_new;

			// variables
			sum = 0;

			// loop through all rows
			for(position = pos_ori; position < pos_ori+d_common.in2_sub_cumh_elem; position = position + d_common.in2_sub_cumh_rows){
				d_unique[bx].d_in2_sub_cumh[position] = d_unique[bx].d_in2_sub_cumh[position] + sum;
				sum = d_unique[bx].d_in2_sub_cumh[position];
			}

			// go for second round
			ei_new = ei_new + NUMBER_THREADS;

		}

		//====================================================================================================
		//	SYNCHRONIZE THREADS
		//====================================================================================================

		__syncthreads();

		//====================================================================================================
		//	SELECTION
		//====================================================================================================

		// work
		ei_new = tx;
		while(ei_new < d_common.in2_sub_cumh_sel_elem){

			// figure out row/col location in new matrix
			row = (ei_new+1) % d_common.in2_sub_cumh_sel_rows - 1;												// (0-n) row
			col = (ei_new+1) / d_common.in2_sub_cumh_sel_rows + 1 - 1;											// (0-n) column
			if((ei_new+1) % d_common.in2_sub_cumh_sel_rows == 0){
				row = d_common.in2_sub_cumh_sel_rows - 1;
				col = col - 1;
			}

			// figure out corresponding location in old matrix and copy values to new matrix
			ori_row = row + d_common.in2_sub_cumh_sel_rowlow - 1;
			ori_col = col + d_common.in2_sub_cumh_sel_collow - 1;
			d_unique[bx].d_in2_sub_cumh_sel[ei_new] = d_unique[bx].d_in2_sub_cumh[ori_col*d_common.in2_sub_cumh_rows+ori_row];

			// go for second round
			ei_new = ei_new + NUMBER_THREADS;

		}

		//====================================================================================================
		//	SYNCHRONIZE THREADS
		//====================================================================================================

		__syncthreads();

		//====================================================================================================
		//	SELECTION 2, SUBTRACTION
		//====================================================================================================

		//==================================================
		//	SELECTION 2
		//==================================================

		// work
		ei_new = tx;
		while(ei_new < d_common.in2_sub2_elem){

			// figure out row/col location in new matrix
			row = (ei_new+1) % d_common.in2_sub2_rows - 1;												// (0-n) row
			col = (ei_new+1) / d_common.in2_sub2_rows + 1 - 1;											// (0-n) column
			if((ei_new+1) % d_common.in2_sub2_rows == 0){
				row = d_common.in2_sub2_rows - 1;
				col = col-1;
			}

			// figure out corresponding location in old matrix and copy values to new matrix
			ori_row = row + d_common.in2_sub_cumh_sel2_rowlow - 1;
			ori_col = col + d_common.in2_sub_cumh_sel2_collow - 1;
			d_unique[bx].d_in2_sub2[ei_new] = d_unique[bx].d_in2_sub_cumh[ori_col*d_common.in2_sub_cumh_rows+ori_row];

			// go for second round
			ei_new = ei_new + NUMBER_THREADS;

		}

		//==================================================
		//	SYNCHRONIZE THREADS
		//==================================================

		__syncthreads();

		//==================================================
		//	SUBTRACTION
		//==================================================

		// work
		ei_new = tx;
		while(ei_new < d_common.in2_sub2_elem){

			// subtract
			d_unique[bx].d_in2_sub2[ei_new] = d_unique[bx].d_in2_sub_cumh_sel[ei_new] - d_unique[bx].d_in2_sub2[ei_new];

			// go for second round
			ei_new = ei_new + NUMBER_THREADS;

		}

		//======================================================================================================================================================
		//	SYNCHRONIZE THREADS
		//======================================================================================================================================================

		__syncthreads();

		//======================================================================================================================================================
		//	CUMULATIVE SUM 2
		//======================================================================================================================================================

		//====================================================================================================
		//	MULTIPLICATION
		//====================================================================================================

		// work
		ei_new = tx;
		while(ei_new < d_common.in2_sqr_elem){

			temp = d_unique[bx].d_in2[ei_new];
			d_unique[bx].d_in2_sqr[ei_new] = temp * temp;

			// go for second round
			ei_new = ei_new + NUMBER_THREADS;

		}

		//====================================================================================================
		//	SYNCHRONIZE THREADS
		//====================================================================================================

		__syncthreads();

		//====================================================================================================
		//	PAD ARRAY, VERTICAL CUMULATIVE SUM
		//====================================================================================================

		//==================================================
		//	PAD ARRAY
		//==================================================

		// work
		ei_new = tx;
		while(ei_new < d_common.in2_pad_cumv_elem){

			// figure out row/col location in padded array
			row = (ei_new+1) % d_common.in2_pad_cumv_rows - 1;												// (0-n) row
			col = (ei_new+1) / d_common.in2_pad_cumv_rows + 1 - 1;											// (0-n) column
			if((ei_new+1) % d_common.in2_pad_cumv_rows == 0){
				row = d_common.in2_pad_cumv_rows - 1;
				col = col-1;
			}

			// execution
			if(	row > (d_common.in2_pad_add_rows-1) &&													// do if has numbers in original array
				row < (d_common.in2_pad_add_rows+d_common.in2_sqr_rows) && 
				col > (d_common.in2_pad_add_cols-1) && 
				col < (d_common.in2_pad_add_cols+d_common.in2_sqr_cols)){
				ori_row = row - d_common.in2_pad_add_rows;
				ori_col = col - d_common.in2_pad_add_cols;
				d_unique[bx].d_in2_pad_cumv[ei_new] = d_unique[bx].d_in2_sqr[ori_col*d_common.in2_sqr_rows+ori_row];
			}
			else{																							// do if otherwise
				d_unique[bx].d_in2_pad_cumv[ei_new] = 0;
			}

			// go for second round
			ei_new = ei_new + NUMBER_THREADS;

		}

		//==================================================
		//	SYNCHRONIZE THREADS
		//==================================================

		__syncthreads();

		//==================================================
		//	VERTICAL CUMULATIVE SUM
		//==================================================

		//work
		ei_new = tx;
		while(ei_new < d_common.in2_pad_cumv_cols){

			// figure out column position
			pos_ori = ei_new*d_common.in2_pad_cumv_rows;

			// variables
			sum = 0;
			
			// loop through all rows
			for(position = pos_ori; position < pos_ori+d_common.in2_pad_cumv_rows; position = position + 1){
				d_unique[bx].d_in2_pad_cumv[position] = d_unique[bx].d_in2_pad_cumv[position] + sum;
				sum = d_unique[bx].d_in2_pad_cumv[position];
			}

			// go for second round
			ei_new = ei_new + NUMBER_THREADS;

		}

		//====================================================================================================
		//	SYNCHRONIZE THREADS
		//====================================================================================================

		__syncthreads();

		//====================================================================================================
		//	SELECTION
		//====================================================================================================

		// work
		ei_new = tx;
		while(ei_new < d_common.in2_pad_cumv_sel_elem){

			// figure out row/col location in new matrix
			row = (ei_new+1) % d_common.in2_pad_cumv_sel_rows - 1;												// (0-n) row
			col = (ei_new+1) / d_common.in2_pad_cumv_sel_rows + 1 - 1;											// (0-n) column
			if((ei_new+1) % d_common.in2_pad_cumv_sel_rows == 0){
				row = d_common.in2_pad_cumv_sel_rows - 1;
				col = col-1;
			}

			// figure out corresponding location in old matrix and copy values to new matrix
			ori_row = row + d_common.in2_pad_cumv_sel_rowlow - 1;
			ori_col = col + d_common.in2_pad_cumv_sel_collow - 1;
			d_unique[bx].d_in2_pad_cumv_sel[ei_new] = d_unique[bx].d_in2_pad_cumv[ori_col*d_common.in2_pad_cumv_rows+ori_row];

			// go for second round
			ei_new = ei_new + NUMBER_THREADS;

		}

		//====================================================================================================
		//	SYNCHRONIZE THREADS
		//====================================================================================================

		__syncthreads();

		//====================================================================================================
		//	SELECTION 2, SUBTRACTION, HORIZONTAL CUMULATIVE SUM
		//====================================================================================================

		//==================================================
		//	SELECTION 2
		//==================================================

		// work
		ei_new = tx;
		while(ei_new < d_common.in2_sub_cumh_elem){

			// figure out row/col location in new matrix
			row = (ei_new+1) % d_common.in2_sub_cumh_rows - 1;												// (0-n) row
			col = (ei_new+1) / d_common.in2_sub_cumh_rows + 1 - 1;											// (0-n) column
			if((ei_new+1) % d_common.in2_sub_cumh_rows == 0){
				row = d_common.in2_sub_cumh_rows - 1;
				col = col-1;
			}

			// figure out corresponding location in old matrix and copy values to new matrix
			ori_row = row + d_common.in2_pad_cumv_sel2_rowlow - 1;
			ori_col = col + d_common.in2_pad_cumv_sel2_collow - 1;
			d_unique[bx].d_in2_sub_cumh[ei_new] = d_unique[bx].d_in2_pad_cumv[ori_col*d_common.in2_pad_cumv_rows+ori_row];

			// go for second round
			ei_new = ei_new + NUMBER_THREADS;

		}

		//==================================================
		//	SYNCHRONIZE THREADS
		//==================================================

		__syncthreads();

		//==================================================
		//	SUBTRACTION
		//==================================================

		// work
		ei_new = tx;
		while(ei_new < d_common.in2_sub_cumh_elem){

			// subtract
			d_unique[bx].d_in2_sub_cumh[ei_new] = d_unique[bx].d_in2_pad_cumv_sel[ei_new] - d_unique[bx].d_in2_sub_cumh[ei_new];

			// go for second round
			ei_new = ei_new + NUMBER_THREADS;

		}

		//==================================================
		//	HORIZONTAL CUMULATIVE SUM
		//==================================================

		// work
		ei_new = tx;
		while(ei_new < d_common.in2_sub_cumh_rows){

			// figure out row position
			pos_ori = ei_new;

			// variables
			sum = 0;

			// loop through all rows
			for(position = pos_ori; position < pos_ori+d_common.in2_sub_cumh_elem; position = position + d_common.in2_sub_cumh_rows){
				d_unique[bx].d_in2_sub_cumh[position] = d_unique[bx].d_in2_sub_cumh[position] + sum;
				sum = d_unique[bx].d_in2_sub_cumh[position];
			}

			// go for second round
			ei_new = ei_new + NUMBER_THREADS;

		}

		//====================================================================================================
		//	SYNCHRONIZE THREADS
		//====================================================================================================

		__syncthreads();

		//====================================================================================================
		//	SELECTION
		//====================================================================================================

		// work
		ei_new = tx;
		while(ei_new < d_common.in2_sub_cumh_sel_elem){

			// figure out row/col location in new matrix
			row = (ei_new+1) % d_common.in2_sub_cumh_sel_rows - 1;												// (0-n) row
			col = (ei_new+1) / d_common.in2_sub_cumh_sel_rows + 1 - 1;											// (0-n) column
			if((ei_new+1) % d_common.in2_sub_cumh_sel_rows == 0){
				row = d_common.in2_sub_cumh_sel_rows - 1;
				col = col - 1;
			}

			// figure out corresponding location in old matrix and copy values to new matrix
			ori_row = row + d_common.in2_sub_cumh_sel_rowlow - 1;
			ori_col = col + d_common.in2_sub_cumh_sel_collow - 1;
			d_unique[bx].d_in2_sub_cumh_sel[ei_new] = d_unique[bx].d_in2_sub_cumh[ori_col*d_common.in2_sub_cumh_rows+ori_row];

			// go for second round
			ei_new = ei_new + NUMBER_THREADS;

		}

		//====================================================================================================
		//	SYNCHRONIZE THREADS
		//====================================================================================================

		__syncthreads();

		//====================================================================================================
		//	SELECTION 2, SUBTRACTION
		//====================================================================================================

		//==================================================
		//	SELECTION 2
		//==================================================

		// work
		ei_new = tx;
		while(ei_new < d_common.in2_sub2_elem){

			// figure out row/col location in new matrix
			row = (ei_new+1) % d_common.in2_sub2_rows - 1;												// (0-n) row
			col = (ei_new+1) / d_common.in2_sub2_rows + 1 - 1;											// (0-n) column
			if((ei_new+1) % d_common.in2_sub2_rows == 0){
				row = d_common.in2_sub2_rows - 1;
				col = col-1;
			}

			// figure out corresponding location in old matrix and copy values to new matrix
			ori_row = row + d_common.in2_sub_cumh_sel2_rowlow - 1;
			ori_col = col + d_common.in2_sub_cumh_sel2_collow - 1;
			d_unique[bx].d_in2_sqr_sub2[ei_new] = d_unique[bx].d_in2_sub_cumh[ori_col*d_common.in2_sub_cumh_rows+ori_row];

			// go for second round
			ei_new = ei_new + NUMBER_THREADS;

		}

		//==================================================
		//	SYNCHRONIZE THREADS
		//==================================================

		__syncthreads();

		//==================================================
		//	SUBTRACTION
		//==================================================

		// work
		ei_new = tx;
		while(ei_new < d_common.in2_sub2_elem){

			// subtract
			d_unique[bx].d_in2_sqr_sub2[ei_new] = d_unique[bx].d_in2_sub_cumh_sel[ei_new] - d_unique[bx].d_in2_sqr_sub2[ei_new];

			// go for second round
			ei_new = ei_new + NUMBER_THREADS;

		}

		//======================================================================================================================================================
		//	SYNCHRONIZE THREADS
		//======================================================================================================================================================

		__syncthreads();

		//======================================================================================================================================================
		//	FINAL
		//======================================================================================================================================================

		//====================================================================================================
		//	DENOMINATOR A		SAVE RESULT IN CUMULATIVE SUM A2
		//====================================================================================================

		// work
		ei_new = tx;
		while(ei_new < d_common.in2_sub2_elem){

			temp = d_unique[bx].d_in2_sub2[ei_new];
			temp2 = d_unique[bx].d_in2_sqr_sub2[ei_new] - (temp * temp / d_common.in_elem);
			if(temp2 < 0){
				temp2 = 0;
			}
			d_unique[bx].d_in2_sqr_sub2[ei_new] = sqrt(temp2);
			

			// go for second round
			ei_new = ei_new + NUMBER_THREADS;

		}

		//====================================================================================================
		//	SYNCHRONIZE THREADS
		//====================================================================================================

		__syncthreads();

		//====================================================================================================
		//	MULTIPLICATION
		//====================================================================================================

		// work
		ei_new = tx;
		while(ei_new < d_common.in_sqr_elem){

			temp = d_in[ei_new];
			d_unique[bx].d_in_sqr[ei_new] = temp * temp;

			// go for second round
			ei_new = ei_new + NUMBER_THREADS;

		}

		//====================================================================================================
		//	SYNCHRONIZE THREADS
		//====================================================================================================

		__syncthreads();

		//====================================================================================================
		//	IN SUM
		//====================================================================================================

		// work
		ei_new = tx;
		while(ei_new < d_common.in_cols){

			sum = 0;
			for(i = 0; i < d_common.in_rows; i++){

				sum = sum + d_in[ei_new*d_common.in_rows+i];

			}
			in_partial_sum[ei_new] = sum;

			// go for second round
			ei_new = ei_new + NUMBER_THREADS;

		}

		//====================================================================================================
		//	SYNCHRONIZE THREADS
		//====================================================================================================

		__syncthreads();

		//====================================================================================================
		//	IN_SQR SUM
		//====================================================================================================

		ei_new = tx;
		while(ei_new < d_common.in_sqr_rows){
				
			sum = 0;
			for(i = 0; i < d_common.in_sqr_cols; i++){

				sum = sum + d_unique[bx].d_in_sqr[ei_new+d_common.in_sqr_rows*i];

			}
			in_sqr_partial_sum[ei_new] = sum;

			// go for second round
			ei_new = ei_new + NUMBER_THREADS;

		}

		//====================================================================================================
		//	SYNCHRONIZE THREADS
		//====================================================================================================

		__syncthreads();

		//====================================================================================================
		//	FINAL SUMMATION
		//====================================================================================================

		if(tx == 0){

			in_final_sum = 0;
			for(i = 0; i<d_common.in_cols; i++){
				in_final_sum = in_final_sum + in_partial_sum[i];
			}

		}else if(tx == 1){

			in_sqr_final_sum = 0;
			for(i = 0; i<d_common.in_sqr_cols; i++){
				in_sqr_final_sum = in_sqr_final_sum + in_sqr_partial_sum[i];
			}

		}

		//====================================================================================================
		//	SYNCHRONIZE THREADS
		//====================================================================================================

		__syncthreads();

		//====================================================================================================
		//	DENOMINATOR T
		//====================================================================================================

		if(tx == 0){

			mean = in_final_sum / d_common.in_elem;													// gets mean (average) value of element in ROI
			mean_sqr = mean * mean;
			variance  = (in_sqr_final_sum / d_common.in_elem) - mean_sqr;							// gets variance of ROI
			deviation = sqrt(variance);																// gets standard deviation of ROI

			denomT = sqrt(float(d_common.in_elem-1))*deviation;

		}

		//====================================================================================================
		//	SYNCHRONIZE THREADS
		//====================================================================================================

		__syncthreads();

		//====================================================================================================
		//	DENOMINATOR		SAVE RESULT IN CUMULATIVE SUM A2
		//====================================================================================================

		// work
		ei_new = tx;
		while(ei_new < d_common.in2_sub2_elem){

			d_unique[bx].d_in2_sqr_sub2[ei_new] = d_unique[bx].d_in2_sqr_sub2[ei_new] * denomT;

			// go for second round
			ei_new = ei_new + NUMBER_THREADS;

		}

		//====================================================================================================
		//	SYNCHRONIZE THREADS
		//====================================================================================================

		__syncthreads();

		//====================================================================================================
		//	NUMERATOR	SAVE RESULT IN CONVOLUTION
		//====================================================================================================

		// work
		ei_new = tx;
		while(ei_new < d_common.conv_elem){

			d_unique[bx].d_conv[ei_new] = d_unique[bx].d_conv[ei_new] - d_unique[bx].d_in2_sub2[ei_new] * in_final_sum / d_common.in_elem;

			// go for second round
			ei_new = ei_new + NUMBER_THREADS;

		}

		//====================================================================================================
		//	SYNCHRONIZE THREADS
		//====================================================================================================

		__syncthreads();

		//====================================================================================================
		//	CORRELATION	SAVE RESULT IN CUMULATIVE SUM A2
		//====================================================================================================

		// work
		ei_new = tx;
		while(ei_new < d_common.in2_sub2_elem){

			d_unique[bx].d_in2_sqr_sub2[ei_new] = d_unique[bx].d_conv[ei_new] / d_unique[bx].d_in2_sqr_sub2[ei_new];

			// go for second round
			ei_new = ei_new + NUMBER_THREADS;

		}

		//======================================================================================================================================================
		//	SYNCHRONIZE THREADS
		//======================================================================================================================================================

		__syncthreads();

		//======================================================================================================================================================
		//	TEMPLATE MASK CREATE
		//======================================================================================================================================================

		cent = d_common.sSize + d_common.tSize + 1;
		if(d_common_change.frame_no == 0){
			tMask_row = cent + d_unique[bx].d_Row[d_unique[bx].point_no] - d_unique[bx].d_Row[d_unique[bx].point_no] - 1;
			tMask_col = cent + d_unique[bx].d_Col[d_unique[bx].point_no] - d_unique[bx].d_Col[d_unique[bx].point_no] - 1;
		}
		else{
			pointer = d_common_change.frame_no-1+d_unique[bx].point_no*d_common.no_frames;
			tMask_row = cent + d_unique[bx].d_tRowLoc[pointer] - d_unique[bx].d_Row[d_unique[bx].point_no] - 1;
			tMask_col = cent + d_unique[bx].d_tColLoc[pointer] - d_unique[bx].d_Col[d_unique[bx].point_no] - 1;
		}


		//work
		ei_new = tx;
		while(ei_new < d_common.tMask_elem){

			location = tMask_col*d_common.tMask_rows + tMask_row;

			if(ei_new==location){
				d_unique[bx].d_tMask[ei_new] = 1;
			}
			else{
				d_unique[bx].d_tMask[ei_new] = 0;
			}

			//go for second round
			ei_new = ei_new + NUMBER_THREADS;

		}

		//======================================================================================================================================================
		//	SYNCHRONIZE THREADS
		//======================================================================================================================================================

		__syncthreads();

		//======================================================================================================================================================
		//	MASK CONVOLUTION
		//======================================================================================================================================================

		// work
		ei_new = tx;
		while(ei_new < d_common.mask_conv_elem){

			// figure out row/col location in array
			ic = (ei_new+1) % d_common.mask_conv_rows;												// (1-n)
			jc = (ei_new+1) / d_common.mask_conv_rows + 1;											// (1-n)
			if((ei_new+1) % d_common.mask_conv_rows == 0){
				ic = d_common.mask_conv_rows;
				jc = jc-1;
			}

			//
			j = jc + d_common.mask_conv_joffset;
			jp1 = j + 1;
			if(d_common.mask_cols < jp1){
				ja1 = jp1 - d_common.mask_cols;
			}
			else{
				ja1 = 1;
			}
			if(d_common.tMask_cols < j){
				ja2 = d_common.tMask_cols;
			}
			else{
				ja2 = j;
			}

			i = ic + d_common.mask_conv_ioffset;
			ip1 = i + 1;
			
			if(d_common.mask_rows < ip1){
				ia1 = ip1 - d_common.mask_rows;
			}
			else{
				ia1 = 1;
			}
			if(d_common.tMask_rows < i){
				ia2 = d_common.tMask_rows;
			}
			else{
				ia2 = i;
			}

			s = 0;

			for(ja=ja1; ja<=ja2; ja++){
				jb = jp1 - ja;
				for(ia=ia1; ia<=ia2; ia++){
					ib = ip1 - ia;
					s = s + d_unique[bx].d_tMask[d_common.tMask_rows*(ja-1)+ia-1] * 1;
				}
			}

			// //d_unique[bx].d_mask_conv[d_common.mask_conv_rows*(jc-1)+ic-1] = s;
			d_unique[bx].d_mask_conv[ei_new] = d_unique[bx].d_in2_sqr_sub2[ei_new] * s;

			// go for second round
			ei_new = ei_new + NUMBER_THREADS;

		}

		//======================================================================================================================================================
		//	SYNCHRONIZE THREADS
		//======================================================================================================================================================

		__syncthreads();

		//======================================================================================================================================================
		//	MAXIMUM VALUE
		//======================================================================================================================================================

		//====================================================================================================
		//	INITIAL SEARCH
		//====================================================================================================

		ei_new = tx;
		while(ei_new < d_common.mask_conv_rows){

			for(i=0; i<d_common.mask_conv_cols; i++){
				largest_coordinate_current = ei_new*d_common.mask_conv_rows+i;
				largest_value_current = abs(d_unique[bx].d_mask_conv[largest_coordinate_current]);
				if(largest_value_current > largest_value){
					largest_coordinate = largest_coordinate_current;
					largest_value = largest_value_current;
				}
			}
			par_max_coo[ei_new] = largest_coordinate;
			par_max_val[ei_new] = largest_value;

			// go for second round
			ei_new = ei_new + NUMBER_THREADS;

		}

		//====================================================================================================
		//	SYNCHRONIZE THREADS
		//====================================================================================================

		__syncthreads();

		//====================================================================================================
		//	FINAL SEARCH
		//====================================================================================================

		if(tx == 0){

			for(i = 0; i < d_common.mask_conv_rows; i++){
				if(par_max_val[i] > fin_max_val){
					fin_max_val = par_max_val[i];
					fin_max_coo = par_max_coo[i];
				}
			}

			// convert coordinate to row/col form
			largest_row = (fin_max_coo+1) % d_common.mask_conv_rows - 1;											// (0-n) row
			largest_col = (fin_max_coo+1) / d_common.mask_conv_rows;												// (0-n) column
			if((fin_max_coo+1) % d_common.mask_conv_rows == 0){
				largest_row = d_common.mask_conv_rows - 1;
				largest_col = largest_col - 1;
			}

			// calculate offset
			largest_row = largest_row + 1;																	// compensate to match MATLAB format (1-n)
			largest_col = largest_col + 1;																	// compensate to match MATLAB format (1-n)
			offset_row = largest_row - d_common.in_rows - (d_common.sSize - d_common.tSize);
			offset_col = largest_col - d_common.in_cols - (d_common.sSize - d_common.tSize);
			pointer = d_common_change.frame_no+d_unique[bx].point_no*d_common.no_frames;
			d_unique[bx].d_tRowLoc[pointer] = d_unique[bx].d_Row[d_unique[bx].point_no] + offset_row;
			d_unique[bx].d_tColLoc[pointer] = d_unique[bx].d_Col[d_unique[bx].point_no] + offset_col;

		}

		//======================================================================================================================================================
		//	SYNCHRONIZE THREADS
		//======================================================================================================================================================

		__syncthreads();

	}
	
	//===============================================================================================================================================================================================================
	//===============================================================================================================================================================================================================
	//	COORDINATE AND TEMPLATE UPDATE
	//===============================================================================================================================================================================================================
	//===============================================================================================================================================================================================================

	// time19 = clock();

	// if the last frame in the bath, update template
	if(d_common_change.frame_no != 0 && (d_common_change.frame_no)%10 == 0){

		// update coordinate
		loc_pointer = d_unique[bx].point_no*d_common.no_frames+d_common_change.frame_no;
		d_unique[bx].d_Row[d_unique[bx].point_no] = d_unique[bx].d_tRowLoc[loc_pointer];
		d_unique[bx].d_Col[d_unique[bx].point_no] = d_unique[bx].d_tColLoc[loc_pointer];

		// work
		ei_new = tx;
		while(ei_new < d_common.in_elem){

			// figure out row/col location in new matrix
			row = (ei_new+1) % d_common.in_rows - 1;												// (0-n) row
			col = (ei_new+1) / d_common.in_rows + 1 - 1;											// (0-n) column
			if((ei_new+1) % d_common.in_rows == 0){
				row = d_common.in_rows - 1;
				col = col-1;
			}

			// figure out row/col location in corresponding new template area in image and give to every thread (get top left corner and progress down and right)
			ori_row = d_unique[bx].d_Row[d_unique[bx].point_no] - 25 + row - 1;
			ori_col = d_unique[bx].d_Col[d_unique[bx].point_no] - 25 + col - 1;
			ori_pointer = ori_col*d_common.frame_rows+ori_row;

			// update template
			d_in[ei_new] = d_common.alpha*d_in[ei_new] + (1.00-d_common.alpha)*d_common_change.d_frame[ori_pointer];

			// go for second round
			ei_new = ei_new + NUMBER_THREADS;

		}

	}

}

	//===============================================================================================================================================================================================================
	//===============================================================================================================================================================================================================
	//	END OF FUNCTION
	//===============================================================================================================================================================================================================
	//===============================================================================================================================================================================================================

int heartmain(int argc, char ** argv) {
	//cudaDeviceSetCacheConfig(cudaFuncCachePreferL1);
	//======================================================================================================================================================
	//	VARIABLES
	//======================================================================================================================================================

	// CUDA kernel execution parameters
	dim3 threads;
	dim3 blocks;

	// counter
	int i;
	int frames_processed;

	// frames
	char* video_file_name;
	avi_t* frames;
	fp* frame;

	//======================================================================================================================================================
	// 	FRAME
	//======================================================================================================================================================

	if(argc!=3){
		printf("ERROR: usage: heartwall <inputfile> <num of frames>\n");
		exit(1);
	}
	
	// open movie file
 	video_file_name = argv[1];
	frames = (avi_t*)AVI_open_input_file(video_file_name, 1);														// added casting
	if (frames == NULL)  {
		   AVI_print_error((char *) "Error with AVI_open_input_file");
		   return -1;
	}

	// common
	common.no_frames = AVI_video_frames(frames);
	common.frame_rows = AVI_video_height(frames);
	common.frame_cols = AVI_video_width(frames);
	common.frame_elem = common.frame_rows * common.frame_cols;
	common.frame_mem = sizeof(fp) * common.frame_elem;

	// pointers
	cudaMalloc((void **)&common_change.d_frame, common.frame_mem);

	//======================================================================================================================================================
	// 	CHECK INPUT ARGUMENTS
	//======================================================================================================================================================
	
	frames_processed = atoi(argv[2]);
		if(frames_processed<0 || frames_processed>common.no_frames){
			printf("ERROR: %d is an incorrect number of frames specified, select in the range of 0-%d\n", frames_processed, common.no_frames);
			return 0;
	}
	

	//======================================================================================================================================================
	//	HARDCODED INPUTS FROM MATLAB
	//======================================================================================================================================================

	//====================================================================================================
	//	CONSTANTS
	//====================================================================================================

	common.sSize = 40;
	common.tSize = 25;
	common.maxMove = 10;
	common.alpha = 0.87;

	//====================================================================================================
	//	ENDO POINTS
	//====================================================================================================

	common.endoPoints = ENDO_POINTS;
	common.endo_mem = sizeof(int) * common.endoPoints;

	common.endoRow = (int *)malloc(common.endo_mem);
	common.endoRow[ 0] = 369;
	common.endoRow[ 1] = 400;
	common.endoRow[ 2] = 429;
	common.endoRow[ 3] = 452;
	common.endoRow[ 4] = 476;
	common.endoRow[ 5] = 486;
	common.endoRow[ 6] = 479;
	common.endoRow[ 7] = 458;
	common.endoRow[ 8] = 433;
	common.endoRow[ 9] = 404;
	common.endoRow[10] = 374;
	common.endoRow[11] = 346;
	common.endoRow[12] = 318;
	common.endoRow[13] = 294;
	common.endoRow[14] = 277;
	common.endoRow[15] = 269;
	common.endoRow[16] = 275;
	common.endoRow[17] = 287;
	common.endoRow[18] = 311;
	common.endoRow[19] = 339;
	cudaMalloc((void **)&common.d_endoRow, common.endo_mem);
	cudaMemcpy(common.d_endoRow, common.endoRow, common.endo_mem, cudaMemcpyHostToDevice);

	common.endoCol = (int *)malloc(common.endo_mem);
	common.endoCol[ 0] = 408;
	common.endoCol[ 1] = 406;
	common.endoCol[ 2] = 397;
	common.endoCol[ 3] = 383;
	common.endoCol[ 4] = 354;
	common.endoCol[ 5] = 322;
	common.endoCol[ 6] = 294;
	common.endoCol[ 7] = 270;
	common.endoCol[ 8] = 250;
	common.endoCol[ 9] = 237;
	common.endoCol[10] = 235;
	common.endoCol[11] = 241;
	common.endoCol[12] = 254;
	common.endoCol[13] = 273;
	common.endoCol[14] = 300;
	common.endoCol[15] = 328;
	common.endoCol[16] = 356;
	common.endoCol[17] = 383;
	common.endoCol[18] = 401;
	common.endoCol[19] = 411;
	cudaMalloc((void **)&common.d_endoCol, common.endo_mem);
	cudaMemcpy(common.d_endoCol, common.endoCol, common.endo_mem, cudaMemcpyHostToDevice);

	common.tEndoRowLoc = (int *)malloc(common.endo_mem * common.no_frames);
	cudaMalloc((void **)&common.d_tEndoRowLoc, common.endo_mem * common.no_frames);

	common.tEndoColLoc = (int *)malloc(common.endo_mem * common.no_frames);
	cudaMalloc((void **)&common.d_tEndoColLoc, common.endo_mem * common.no_frames);

	//====================================================================================================
	//	EPI POINTS
	//====================================================================================================

	common.epiPoints = EPI_POINTS;
	common.epi_mem = sizeof(int) * common.epiPoints;

	common.epiRow = (int *)malloc(common.epi_mem);
	common.epiRow[ 0] = 390;
	common.epiRow[ 1] = 419;
	common.epiRow[ 2] = 448;
	common.epiRow[ 3] = 474;
	common.epiRow[ 4] = 501;
	common.epiRow[ 5] = 519;
	common.epiRow[ 6] = 535;
	common.epiRow[ 7] = 542;
	common.epiRow[ 8] = 543;
	common.epiRow[ 9] = 538;
	common.epiRow[10] = 528;
	common.epiRow[11] = 511;
	common.epiRow[12] = 491;
	common.epiRow[13] = 466;
	common.epiRow[14] = 438;
	common.epiRow[15] = 406;
	common.epiRow[16] = 376;
	common.epiRow[17] = 347;
	common.epiRow[18] = 318;
	common.epiRow[19] = 291;
	common.epiRow[20] = 275;
	common.epiRow[21] = 259;
	common.epiRow[22] = 256;
	common.epiRow[23] = 252;
	common.epiRow[24] = 252;
	common.epiRow[25] = 257;
	common.epiRow[26] = 266;
	common.epiRow[27] = 283;
	common.epiRow[28] = 305;
	common.epiRow[29] = 331;
	common.epiRow[30] = 360;
	cudaMalloc((void **)&common.d_epiRow, common.epi_mem);
	cudaMemcpy(common.d_epiRow, common.epiRow, common.epi_mem, cudaMemcpyHostToDevice);

	common.epiCol = (int *)malloc(common.epi_mem);
	common.epiCol[ 0] = 457;
	common.epiCol[ 1] = 454;
	common.epiCol[ 2] = 446;
	common.epiCol[ 3] = 431;
	common.epiCol[ 4] = 411;
	common.epiCol[ 5] = 388;
	common.epiCol[ 6] = 361;
	common.epiCol[ 7] = 331;
	common.epiCol[ 8] = 301;
	common.epiCol[ 9] = 273;
	common.epiCol[10] = 243;
	common.epiCol[11] = 218;
	common.epiCol[12] = 196;
	common.epiCol[13] = 178;
	common.epiCol[14] = 166;
	common.epiCol[15] = 157;
	common.epiCol[16] = 155;
	common.epiCol[17] = 165;
	common.epiCol[18] = 177;
	common.epiCol[19] = 197;
	common.epiCol[20] = 218;
	common.epiCol[21] = 248;
	common.epiCol[22] = 276;
	common.epiCol[23] = 304;
	common.epiCol[24] = 333;
	common.epiCol[25] = 361;
	common.epiCol[26] = 391;
	common.epiCol[27] = 415;
	common.epiCol[28] = 434;
	common.epiCol[29] = 448;
	common.epiCol[30] = 455;
	cudaMalloc((void **)&common.d_epiCol, common.epi_mem);
	cudaMemcpy(common.d_epiCol, common.epiCol, common.epi_mem, cudaMemcpyHostToDevice);

	common.tEpiRowLoc = (int *)malloc(common.epi_mem * common.no_frames);
	cudaMalloc((void **)&common.d_tEpiRowLoc, common.epi_mem * common.no_frames);

	common.tEpiColLoc = (int *)malloc(common.epi_mem * common.no_frames);
	cudaMalloc((void **)&common.d_tEpiColLoc, common.epi_mem * common.no_frames);

	//====================================================================================================
	//	ALL POINTS
	//====================================================================================================

	common.allPoints = ALL_POINTS;

	//======================================================================================================================================================
	// 	TEMPLATE SIZES
	//======================================================================================================================================================

	// common
	common.in_rows = common.tSize + 1 + common.tSize;
	common.in_cols = common.in_rows;
	common.in_elem = common.in_rows * common.in_cols;
	common.in_mem = sizeof(fp) * common.in_elem;

	//======================================================================================================================================================
	// 	CREATE ARRAY OF TEMPLATES FOR ALL POINTS
	//======================================================================================================================================================

	// common
	cudaMalloc((void **)&common.d_endoT, common.in_mem * common.endoPoints);
	cudaMalloc((void **)&common.d_epiT, common.in_mem * common.epiPoints);

	//======================================================================================================================================================
	//	SPECIFIC TO ENDO OR EPI TO BE SET HERE
	//======================================================================================================================================================

	for(i=0; i<common.endoPoints; i++){
		unique[i].point_no = i;
		unique[i].d_Row = common.d_endoRow;
		unique[i].d_Col = common.d_endoCol;
		unique[i].d_tRowLoc = common.d_tEndoRowLoc;
		unique[i].d_tColLoc = common.d_tEndoColLoc;
		unique[i].d_T = common.d_endoT;
	}
	for(i=common.endoPoints; i<common.allPoints; i++){
		unique[i].point_no = i-common.endoPoints;
		unique[i].d_Row = common.d_epiRow;
		unique[i].d_Col = common.d_epiCol;
		unique[i].d_tRowLoc = common.d_tEpiRowLoc;
		unique[i].d_tColLoc = common.d_tEpiColLoc;
		unique[i].d_T = common.d_epiT;
	}

	//======================================================================================================================================================
	// 	RIGHT TEMPLATE 	FROM 	TEMPLATE ARRAY
	//======================================================================================================================================================

	// pointers
	for(i=0; i<common.allPoints; i++){
		unique[i].in_pointer = unique[i].point_no * common.in_elem;
	}

	//======================================================================================================================================================
	// 	AREA AROUND POINT		FROM	FRAME
	//======================================================================================================================================================

	// common
	common.in2_rows = 2 * common.sSize + 1;
	common.in2_cols = 2 * common.sSize + 1;
	common.in2_elem = common.in2_rows * common.in2_cols;
	common.in2_mem = sizeof(float) * common.in2_elem;

	// pointers
	for(i=0; i<common.allPoints; i++){
		cudaMalloc((void **)&unique[i].d_in2, common.in2_mem);
	}

	//======================================================================================================================================================
	// 	CONVOLUTION
	//======================================================================================================================================================

	// common
	common.conv_rows = common.in_rows + common.in2_rows - 1;												// number of rows in I
	common.conv_cols = common.in_cols + common.in2_cols - 1;												// number of columns in I
	common.conv_elem = common.conv_rows * common.conv_cols;												// number of elements
	common.conv_mem = sizeof(float) * common.conv_elem;
	common.ioffset = 0;
	common.joffset = 0;

	// pointers
	for(i=0; i<common.allPoints; i++){
		cudaMalloc((void **)&unique[i].d_conv, common.conv_mem);
	}

	//======================================================================================================================================================
	// 	CUMULATIVE SUM
	//======================================================================================================================================================

	//====================================================================================================
	// 	PADDING OF ARRAY, VERTICAL CUMULATIVE SUM
	//====================================================================================================

	// common
	common.in2_pad_add_rows = common.in_rows;
	common.in2_pad_add_cols = common.in_cols;

	common.in2_pad_cumv_rows = common.in2_rows + 2*common.in2_pad_add_rows;
	common.in2_pad_cumv_cols = common.in2_cols + 2*common.in2_pad_add_cols;
	common.in2_pad_cumv_elem = common.in2_pad_cumv_rows * common.in2_pad_cumv_cols;
	common.in2_pad_cumv_mem = sizeof(float) * common.in2_pad_cumv_elem;

	// pointers
	for(i=0; i<common.allPoints; i++){
		cudaMalloc((void **)&unique[i].d_in2_pad_cumv, common.in2_pad_cumv_mem);
	}

	//====================================================================================================
	// 	SELECTION
	//====================================================================================================

	// common
	common.in2_pad_cumv_sel_rowlow = 1 + common.in_rows;													// (1 to n+1)
	common.in2_pad_cumv_sel_rowhig = common.in2_pad_cumv_rows - 1;
	common.in2_pad_cumv_sel_collow = 1;
	common.in2_pad_cumv_sel_colhig = common.in2_pad_cumv_cols;
	common.in2_pad_cumv_sel_rows = common.in2_pad_cumv_sel_rowhig - common.in2_pad_cumv_sel_rowlow + 1;
	common.in2_pad_cumv_sel_cols = common.in2_pad_cumv_sel_colhig - common.in2_pad_cumv_sel_collow + 1;
	common.in2_pad_cumv_sel_elem = common.in2_pad_cumv_sel_rows * common.in2_pad_cumv_sel_cols;
	common.in2_pad_cumv_sel_mem = sizeof(float) * common.in2_pad_cumv_sel_elem;

	// pointers
	for(i=0; i<common.allPoints; i++){
		cudaMalloc((void **)&unique[i].d_in2_pad_cumv_sel, common.in2_pad_cumv_sel_mem);
	}

	//====================================================================================================
	// 	SELECTION	2, SUBTRACTION, HORIZONTAL CUMULATIVE SUM
	//====================================================================================================

	// common
	common.in2_pad_cumv_sel2_rowlow = 1;
	common.in2_pad_cumv_sel2_rowhig = common.in2_pad_cumv_rows - common.in_rows - 1;
	common.in2_pad_cumv_sel2_collow = 1;
	common.in2_pad_cumv_sel2_colhig = common.in2_pad_cumv_cols;
	common.in2_sub_cumh_rows = common.in2_pad_cumv_sel2_rowhig - common.in2_pad_cumv_sel2_rowlow + 1;
	common.in2_sub_cumh_cols = common.in2_pad_cumv_sel2_colhig - common.in2_pad_cumv_sel2_collow + 1;
	common.in2_sub_cumh_elem = common.in2_sub_cumh_rows * common.in2_sub_cumh_cols;
	common.in2_sub_cumh_mem = sizeof(float) * common.in2_sub_cumh_elem;

	// pointers
	for(i=0; i<common.allPoints; i++){
		cudaMalloc((void **)&unique[i].d_in2_sub_cumh, common.in2_sub_cumh_mem);
	}

	//====================================================================================================
	// 	SELECTION
	//====================================================================================================

	// common
	common.in2_sub_cumh_sel_rowlow = 1;
	common.in2_sub_cumh_sel_rowhig = common.in2_sub_cumh_rows;
	common.in2_sub_cumh_sel_collow = 1 + common.in_cols;
	common.in2_sub_cumh_sel_colhig = common.in2_sub_cumh_cols - 1;
	common.in2_sub_cumh_sel_rows = common.in2_sub_cumh_sel_rowhig - common.in2_sub_cumh_sel_rowlow + 1;
	common.in2_sub_cumh_sel_cols = common.in2_sub_cumh_sel_colhig - common.in2_sub_cumh_sel_collow + 1;
	common.in2_sub_cumh_sel_elem = common.in2_sub_cumh_sel_rows * common.in2_sub_cumh_sel_cols;
	common.in2_sub_cumh_sel_mem = sizeof(float) * common.in2_sub_cumh_sel_elem;

	// pointers
	for(i=0; i<common.allPoints; i++){
		cudaMalloc((void **)&unique[i].d_in2_sub_cumh_sel, common.in2_sub_cumh_sel_mem);
	}

	//====================================================================================================
	//	SELECTION 2, SUBTRACTION
	//====================================================================================================

	// common
	common.in2_sub_cumh_sel2_rowlow = 1;
	common.in2_sub_cumh_sel2_rowhig = common.in2_sub_cumh_rows;
	common.in2_sub_cumh_sel2_collow = 1;
	common.in2_sub_cumh_sel2_colhig = common.in2_sub_cumh_cols - common.in_cols - 1;
	common.in2_sub2_rows = common.in2_sub_cumh_sel2_rowhig - common.in2_sub_cumh_sel2_rowlow + 1;
	common.in2_sub2_cols = common.in2_sub_cumh_sel2_colhig - common.in2_sub_cumh_sel2_collow + 1;
	common.in2_sub2_elem = common.in2_sub2_rows * common.in2_sub2_cols;
	common.in2_sub2_mem = sizeof(float) * common.in2_sub2_elem;

	// pointers
	for(i=0; i<common.allPoints; i++){
		cudaMalloc((void **)&unique[i].d_in2_sub2, common.in2_sub2_mem);
	}

	//======================================================================================================================================================
	//	CUMULATIVE SUM 2
	//======================================================================================================================================================

	//====================================================================================================
	//	MULTIPLICATION
	//====================================================================================================

	// common
	common.in2_sqr_rows = common.in2_rows;
	common.in2_sqr_cols = common.in2_cols;
	common.in2_sqr_elem = common.in2_elem;
	common.in2_sqr_mem = common.in2_mem;

	// pointers
	for(i=0; i<common.allPoints; i++){
		cudaMalloc((void **)&unique[i].d_in2_sqr, common.in2_sqr_mem);
	}

	//====================================================================================================
	//	SELECTION 2, SUBTRACTION
	//====================================================================================================

	// common
	common.in2_sqr_sub2_rows = common.in2_sub2_rows;
	common.in2_sqr_sub2_cols = common.in2_sub2_cols;
	common.in2_sqr_sub2_elem = common.in2_sub2_elem;
	common.in2_sqr_sub2_mem = common.in2_sub2_mem;

	// pointers
	for(i=0; i<common.allPoints; i++){
		cudaMalloc((void **)&unique[i].d_in2_sqr_sub2, common.in2_sqr_sub2_mem);
	}

	//======================================================================================================================================================
	//	FINAL
	//======================================================================================================================================================

	// common
	common.in_sqr_rows = common.in_rows;
	common.in_sqr_cols = common.in_cols;
	common.in_sqr_elem = common.in_elem;
	common.in_sqr_mem = common.in_mem;

	// pointers
	for(i=0; i<common.allPoints; i++){
		cudaMalloc((void **)&unique[i].d_in_sqr, common.in_sqr_mem);
	}

	//======================================================================================================================================================
	//	TEMPLATE MASK CREATE
	//======================================================================================================================================================

	// common
	common.tMask_rows = common.in_rows + (common.sSize+1+common.sSize) - 1;
	common.tMask_cols = common.tMask_rows;
	common.tMask_elem = common.tMask_rows * common.tMask_cols;
	common.tMask_mem = sizeof(float) * common.tMask_elem;

	// pointers
	for(i=0; i<common.allPoints; i++){
		cudaMalloc((void **)&unique[i].d_tMask, common.tMask_mem);
	}

	//======================================================================================================================================================
	//	POINT MASK INITIALIZE
	//======================================================================================================================================================

	// common
	common.mask_rows = common.maxMove;
	common.mask_cols = common.mask_rows;
	common.mask_elem = common.mask_rows * common.mask_cols;
	common.mask_mem = sizeof(float) * common.mask_elem;

	//======================================================================================================================================================
	//	MASK CONVOLUTION
	//======================================================================================================================================================

	// common
	common.mask_conv_rows = common.tMask_rows;												// number of rows in I
	common.mask_conv_cols = common.tMask_cols;												// number of columns in I
	common.mask_conv_elem = common.mask_conv_rows * common.mask_conv_cols;												// number of elements
	common.mask_conv_mem = sizeof(float) * common.mask_conv_elem;
	common.mask_conv_ioffset = (common.mask_rows-1)/2;
	if((common.mask_rows-1) % 2 > 0.5){
		common.mask_conv_ioffset = common.mask_conv_ioffset + 1;
	}
	common.mask_conv_joffset = (common.mask_cols-1)/2;
	if((common.mask_cols-1) % 2 > 0.5){
		common.mask_conv_joffset = common.mask_conv_joffset + 1;
	}

	// pointers
	for(i=0; i<common.allPoints; i++){
		cudaMalloc((void **)&unique[i].d_mask_conv, common.mask_conv_mem);
	}

	//======================================================================================================================================================
	//	KERNEL
	//======================================================================================================================================================

	//====================================================================================================
	//	THREAD BLOCK
	//====================================================================================================

	// All kernels operations within kernel use same max size of threads. Size of block size is set to the size appropriate for max size operation (on padded matrix). Other use subsets of that.
	threads.x = NUMBER_THREADS;											// define the number of threads in the block
	threads.y = 1;
	blocks.x = common.allPoints;							// define the number of blocks in the grid
	blocks.y = 1;

	//====================================================================================================
	//	COPY ARGUMENTS
	//====================================================================================================

	cudaMemcpyToSymbol(d_common, &common, sizeof(params_common));
	cudaMemcpyToSymbol(d_unique, &unique, sizeof(params_unique)*ALL_POINTS);

	//====================================================================================================
	//	PRINT FRAME PROGRESS START
	//====================================================================================================

	//printf("frame progress: ");
	//fflush(NULL);

	//====================================================================================================
	//	LAUNCH
	//====================================================================================================

    cudaThreadSynchronize();
    
    ///
	common_change.frame_no=0;
	frame = get_frame(	frames,						// pointer to video file
										common_change.frame_no,				// number of frame that needs to be returned
										0,								// cropped?
										0,								// scaled?
										1);							// converted
	cudaMemcpy(common_change.d_frame, frame, common.frame_mem, cudaMemcpyHostToDevice);
	cudaMemcpyToSymbol(d_common_change, &common_change, sizeof(params_common_change));
	
    ///
    
    /*long long start_time; // = get_time();
    long long end_time; // = get_time();
    double time_elapsed = 0.0; // = elapsed_time(start_time, end_time);

	
		printf("blocks: %d, %d, %d, threads: %d, %d, %d\n", blocks.x, blocks.y, blocks.z, threads.x, threads.y, threads.z);

	/*for(common_change.frame_no=0; common_change.frame_no<frames_processed; common_change.frame_no++){

		// Extract a cropped version of the first frame from the video file
		frame = get_frame(	frames,						// pointer to video file
										common_change.frame_no,				// number of frame that needs to be returned
										0,								// cropped?
										0,								// scaled?
										1);							// converted

		// copy frame to GPU memory
		cudaMemcpy(common_change.d_frame, frame, common.frame_mem, cudaMemcpyHostToDevice);
		cudaMemcpyToSymbol(d_common_change, &common_change, sizeof(params_common_change));




		cudaThreadSynchronize();
		start_time = get_time();

		// launch GPU kernel
		heartwall_kernel<<<blocks, threads>>>();

		cudaThreadSynchronize();
		
		end_time = get_time();

		cudaError_t errmsg = cudaGetLastError();

    if ( cudaSuccess != errmsg ) 
		{
				printf("Error msg: %d, cudaGetErrorString: %s\n", errmsg, cudaGetErrorString(errmsg));
        printf( "Kernel Execution Error!\n" );
		}

		time_elapsed += elapsed_time(start_time, end_time);


		// free frame after each loop iteration, since AVI library allocates memory for every frame fetched
		free(frame);

		// print frame progress
		printf("%d (time %lf)", common_change.frame_no, time_elapsed);
		fflush(NULL);

	}

	//====================================================================================================
	//	PRINT FRAME PROGRESS END
	//====================================================================================================

	printf("\n");
  printf("HIRREG time: %lf\n", time_elapsed);
	fflush(NULL);*/
}

__device__ void cfdkernel(int nelr, int* elements_surrounding_elements, float* normals, float* variables, float* fluxes)
{
	const float smoothing_coefficient = float(0.2f);
	const int i = (blockDim.x*blockIdx.x + threadIdx.x);

	int j, nb;
	float3 normal; float normal_len;
	float factor;

	float density_i = variables[i + VAR_DENSITY*nelr];
	float3 momentum_i;
	momentum_i.x = variables[i + (VAR_MOMENTUM+0)*nelr];
	momentum_i.y = variables[i + (VAR_MOMENTUM+1)*nelr];
	momentum_i.z = variables[i + (VAR_MOMENTUM+2)*nelr];

	float density_energy_i = variables[i + VAR_DENSITY_ENERGY*nelr];

	float3 velocity_i;             				compute_velocity(density_i, momentum_i, velocity_i);
	float speed_sqd_i                          = compute_speed_sqd(velocity_i);
	float speed_i                              = sqrtf(speed_sqd_i);
	float pressure_i                           = compute_pressure(density_i, density_energy_i, speed_sqd_i);
	float speed_of_sound_i                     = compute_speed_of_sound(density_i, pressure_i);
	float3 flux_contribution_i_momentum_x, flux_contribution_i_momentum_y, flux_contribution_i_momentum_z;
	float3 flux_contribution_i_density_energy;	
	compute_flux_contribution(density_i, momentum_i, density_energy_i, pressure_i, velocity_i, flux_contribution_i_momentum_x, flux_contribution_i_momentum_y, flux_contribution_i_momentum_z, flux_contribution_i_density_energy);

	float flux_i_density = float(0.0f);
	float3 flux_i_momentum;
	flux_i_momentum.x = float(0.0f);
	flux_i_momentum.y = float(0.0f);
	flux_i_momentum.z = float(0.0f);
	float flux_i_density_energy = float(0.0f);
	
	float3 velocity_nb;
	float density_nb, density_energy_nb;
	float3 momentum_nb;
	float3 flux_contribution_nb_momentum_x, flux_contribution_nb_momentum_y, flux_contribution_nb_momentum_z;
	float3 flux_contribution_nb_density_energy;	
	float speed_sqd_nb, speed_of_sound_nb, pressure_nb;

	#pragma unroll
	for(j = 0; j < NNB; j++)
	{
		nb = elements_surrounding_elements[i + j*nelr];
		normal.x = normals[i + (j + 0*NNB)*nelr];
		normal.y = normals[i + (j + 1*NNB)*nelr];
		normal.z = normals[i + (j + 2*NNB)*nelr];
		normal_len = sqrtf(normal.x*normal.x + normal.y*normal.y + normal.z*normal.z);
	
		if(nb >= 0) 	// a legitimate neighbor
		{
			density_nb = variables[nb + VAR_DENSITY*nelr];
			momentum_nb.x = variables[nb + (VAR_MOMENTUM+0)*nelr];
			momentum_nb.y = variables[nb + (VAR_MOMENTUM+1)*nelr];
			momentum_nb.z = variables[nb + (VAR_MOMENTUM+2)*nelr];
			density_energy_nb = variables[nb + VAR_DENSITY_ENERGY*nelr];
												compute_velocity(density_nb, momentum_nb, velocity_nb);
			speed_sqd_nb                      = compute_speed_sqd(velocity_nb);
			pressure_nb                       = compute_pressure(density_nb, density_energy_nb, speed_sqd_nb);
			speed_of_sound_nb                 = compute_speed_of_sound(density_nb, pressure_nb);
				                            compute_flux_contribution(density_nb, momentum_nb, density_energy_nb, pressure_nb, velocity_nb, flux_contribution_nb_momentum_x, flux_contribution_nb_momentum_y, flux_contribution_nb_momentum_z, flux_contribution_nb_density_energy);
		
			// artificial viscosity
			factor = -normal_len*smoothing_coefficient*float(0.5f)*(speed_i + sqrtf(speed_sqd_nb) + speed_of_sound_i + speed_of_sound_nb);
			flux_i_density += factor*(density_i-density_nb);
			flux_i_density_energy += factor*(density_energy_i-density_energy_nb);
			flux_i_momentum.x += factor*(momentum_i.x-momentum_nb.x);
			flux_i_momentum.y += factor*(momentum_i.y-momentum_nb.y);
			flux_i_momentum.z += factor*(momentum_i.z-momentum_nb.z);

			// accumulate cell-centered fluxes
			factor = float(0.5f)*normal.x;
			flux_i_density += factor*(momentum_nb.x+momentum_i.x);
			flux_i_density_energy += factor*(flux_contribution_nb_density_energy.x+flux_contribution_i_density_energy.x);
			flux_i_momentum.x += factor*(flux_contribution_nb_momentum_x.x+flux_contribution_i_momentum_x.x);
			flux_i_momentum.y += factor*(flux_contribution_nb_momentum_y.x+flux_contribution_i_momentum_y.x);
			flux_i_momentum.z += factor*(flux_contribution_nb_momentum_z.x+flux_contribution_i_momentum_z.x);
		
			factor = float(0.5f)*normal.y;
			flux_i_density += factor*(momentum_nb.y+momentum_i.y);
			flux_i_density_energy += factor*(flux_contribution_nb_density_energy.y+flux_contribution_i_density_energy.y);
			flux_i_momentum.x += factor*(flux_contribution_nb_momentum_x.y+flux_contribution_i_momentum_x.y);
			flux_i_momentum.y += factor*(flux_contribution_nb_momentum_y.y+flux_contribution_i_momentum_y.y);
			flux_i_momentum.z += factor*(flux_contribution_nb_momentum_z.y+flux_contribution_i_momentum_z.y);
		
			factor = float(0.5f)*normal.z;
			flux_i_density += factor*(momentum_nb.z+momentum_i.z);
			flux_i_density_energy += factor*(flux_contribution_nb_density_energy.z+flux_contribution_i_density_energy.z);
			flux_i_momentum.x += factor*(flux_contribution_nb_momentum_x.z+flux_contribution_i_momentum_x.z);
			flux_i_momentum.y += factor*(flux_contribution_nb_momentum_y.z+flux_contribution_i_momentum_y.z);
			flux_i_momentum.z += factor*(flux_contribution_nb_momentum_z.z+flux_contribution_i_momentum_z.z);
		}
		else if(nb == -1)	// a wing boundary
		{
			flux_i_momentum.x += normal.x*pressure_i;
			flux_i_momentum.y += normal.y*pressure_i;
			flux_i_momentum.z += normal.z*pressure_i;
		}
		else if(nb == -2) // a far field boundary
		{
			factor = float(0.5f)*normal.x;
			flux_i_density += factor*(ff_variable[VAR_MOMENTUM+0]+momentum_i.x);
			flux_i_density_energy += factor*(ff_flux_contribution_density_energy[0].x+flux_contribution_i_density_energy.x);
			flux_i_momentum.x += factor*(ff_flux_contribution_momentum_x[0].x + flux_contribution_i_momentum_x.x);
			flux_i_momentum.y += factor*(ff_flux_contribution_momentum_y[0].x + flux_contribution_i_momentum_y.x);
			flux_i_momentum.z += factor*(ff_flux_contribution_momentum_z[0].x + flux_contribution_i_momentum_z.x);
		
			factor = float(0.5f)*normal.y;
			flux_i_density += factor*(ff_variable[VAR_MOMENTUM+1]+momentum_i.y);
			flux_i_density_energy += factor*(ff_flux_contribution_density_energy[0].y+flux_contribution_i_density_energy.y);
			flux_i_momentum.x += factor*(ff_flux_contribution_momentum_x[0].y + flux_contribution_i_momentum_x.y);
			flux_i_momentum.y += factor*(ff_flux_contribution_momentum_y[0].y + flux_contribution_i_momentum_y.y);
			flux_i_momentum.z += factor*(ff_flux_contribution_momentum_z[0].y + flux_contribution_i_momentum_z.y);

			factor = float(0.5f)*normal.z;
			flux_i_density += factor*(ff_variable[VAR_MOMENTUM+2]+momentum_i.z);
			flux_i_density_energy += factor*(ff_flux_contribution_density_energy[0].z+flux_contribution_i_density_energy.z);
			flux_i_momentum.x += factor*(ff_flux_contribution_momentum_x[0].z + flux_contribution_i_momentum_x.z);
			flux_i_momentum.y += factor*(ff_flux_contribution_momentum_y[0].z + flux_contribution_i_momentum_y.z);
			flux_i_momentum.z += factor*(ff_flux_contribution_momentum_z[0].z + flux_contribution_i_momentum_z.z);

		}
	}

	fluxes[i + VAR_DENSITY*nelr] = flux_i_density;
	fluxes[i + (VAR_MOMENTUM+0)*nelr] = flux_i_momentum.x;
	fluxes[i + (VAR_MOMENTUM+1)*nelr] = flux_i_momentum.y;
	fluxes[i + (VAR_MOMENTUM+2)*nelr] = flux_i_momentum.z;
	fluxes[i + VAR_DENSITY_ENERGY*nelr] = flux_i_density_energy;
}

__global__ void kernel(int nelr, int* elements_surrounding_elements, float* normals, float* variables, float* fluxes)
{
	if(threadIdx.x % 3 == 1) {
		const float smoothing_coefficient = float(0.2f);
		const int i = (blockDim.x*blockIdx.x + threadIdx.x);
	
		int j, nb;
		float3 normal; float normal_len;
		float factor;
	
		float density_i = variables[i + VAR_DENSITY*nelr];
		float3 momentum_i;
		momentum_i.x = variables[i + (VAR_MOMENTUM+0)*nelr];
		momentum_i.y = variables[i + (VAR_MOMENTUM+1)*nelr];
		momentum_i.z = variables[i + (VAR_MOMENTUM+2)*nelr];

		float density_energy_i = variables[i + VAR_DENSITY_ENERGY*nelr];

		float3 velocity_i;             				compute_velocity(density_i, momentum_i, velocity_i);
		float speed_sqd_i                          = compute_speed_sqd(velocity_i);
		float speed_i                              = sqrtf(speed_sqd_i);
		float pressure_i                           = compute_pressure(density_i, density_energy_i, speed_sqd_i);
		float speed_of_sound_i                     = compute_speed_of_sound(density_i, pressure_i);
		float3 flux_contribution_i_momentum_x, flux_contribution_i_momentum_y, flux_contribution_i_momentum_z;
		float3 flux_contribution_i_density_energy;	
		compute_flux_contribution(density_i, momentum_i, density_energy_i, pressure_i, velocity_i, flux_contribution_i_momentum_x, flux_contribution_i_momentum_y, flux_contribution_i_momentum_z, flux_contribution_i_density_energy);
	
		float flux_i_density = float(0.0f);
		float3 flux_i_momentum;
		flux_i_momentum.x = float(0.0f);
		flux_i_momentum.y = float(0.0f);
		flux_i_momentum.z = float(0.0f);
		float flux_i_density_energy = float(0.0f);
		
		float3 velocity_nb;
		float density_nb, density_energy_nb;
		float3 momentum_nb;
		float3 flux_contribution_nb_momentum_x, flux_contribution_nb_momentum_y, flux_contribution_nb_momentum_z;
		float3 flux_contribution_nb_density_energy;	
		float speed_sqd_nb, speed_of_sound_nb, pressure_nb;
	
		#pragma unroll
		for(j = 0; j < NNB; j++)
		{
			nb = elements_surrounding_elements[i + j*nelr];
			normal.x = normals[i + (j + 0*NNB)*nelr];
			normal.y = normals[i + (j + 1*NNB)*nelr];
			normal.z = normals[i + (j + 2*NNB)*nelr];
			normal_len = sqrtf(normal.x*normal.x + normal.y*normal.y + normal.z*normal.z);
		
			if(nb >= 0) 	// a legitimate neighbor
			{
				density_nb = variables[nb + VAR_DENSITY*nelr];
				momentum_nb.x = variables[nb + (VAR_MOMENTUM+0)*nelr];
				momentum_nb.y = variables[nb + (VAR_MOMENTUM+1)*nelr];
				momentum_nb.z = variables[nb + (VAR_MOMENTUM+2)*nelr];
				density_energy_nb = variables[nb + VAR_DENSITY_ENERGY*nelr];
													compute_velocity(density_nb, momentum_nb, velocity_nb);
				speed_sqd_nb                      = compute_speed_sqd(velocity_nb);
				pressure_nb                       = compute_pressure(density_nb, density_energy_nb, speed_sqd_nb);
				speed_of_sound_nb                 = compute_speed_of_sound(density_nb, pressure_nb);
					                            compute_flux_contribution(density_nb, momentum_nb, density_energy_nb, pressure_nb, velocity_nb, flux_contribution_nb_momentum_x, flux_contribution_nb_momentum_y, flux_contribution_nb_momentum_z, flux_contribution_nb_density_energy);
			
				// artificial viscosity
				factor = -normal_len*smoothing_coefficient*float(0.5f)*(speed_i + sqrtf(speed_sqd_nb) + speed_of_sound_i + speed_of_sound_nb);
				flux_i_density += factor*(density_i-density_nb);
				flux_i_density_energy += factor*(density_energy_i-density_energy_nb);
				flux_i_momentum.x += factor*(momentum_i.x-momentum_nb.x);
				flux_i_momentum.y += factor*(momentum_i.y-momentum_nb.y);
				flux_i_momentum.z += factor*(momentum_i.z-momentum_nb.z);

				// accumulate cell-centered fluxes
				factor = float(0.5f)*normal.x;
				flux_i_density += factor*(momentum_nb.x+momentum_i.x);
				flux_i_density_energy += factor*(flux_contribution_nb_density_energy.x+flux_contribution_i_density_energy.x);
				flux_i_momentum.x += factor*(flux_contribution_nb_momentum_x.x+flux_contribution_i_momentum_x.x);
				flux_i_momentum.y += factor*(flux_contribution_nb_momentum_y.x+flux_contribution_i_momentum_y.x);
				flux_i_momentum.z += factor*(flux_contribution_nb_momentum_z.x+flux_contribution_i_momentum_z.x);
			
				factor = float(0.5f)*normal.y;
				flux_i_density += factor*(momentum_nb.y+momentum_i.y);
				flux_i_density_energy += factor*(flux_contribution_nb_density_energy.y+flux_contribution_i_density_energy.y);
				flux_i_momentum.x += factor*(flux_contribution_nb_momentum_x.y+flux_contribution_i_momentum_x.y);
				flux_i_momentum.y += factor*(flux_contribution_nb_momentum_y.y+flux_contribution_i_momentum_y.y);
				flux_i_momentum.z += factor*(flux_contribution_nb_momentum_z.y+flux_contribution_i_momentum_z.y);
			
				factor = float(0.5f)*normal.z;
				flux_i_density += factor*(momentum_nb.z+momentum_i.z);
				flux_i_density_energy += factor*(flux_contribution_nb_density_energy.z+flux_contribution_i_density_energy.z);
				flux_i_momentum.x += factor*(flux_contribution_nb_momentum_x.z+flux_contribution_i_momentum_x.z);
				flux_i_momentum.y += factor*(flux_contribution_nb_momentum_y.z+flux_contribution_i_momentum_y.z);
				flux_i_momentum.z += factor*(flux_contribution_nb_momentum_z.z+flux_contribution_i_momentum_z.z);
			}
			else if(nb == -1)	// a wing boundary
			{
				flux_i_momentum.x += normal.x*pressure_i;
				flux_i_momentum.y += normal.y*pressure_i;
				flux_i_momentum.z += normal.z*pressure_i;
			}
			else if(nb == -2) // a far field boundary
			{
				factor = float(0.5f)*normal.x;
				flux_i_density += factor*(ff_variable[VAR_MOMENTUM+0]+momentum_i.x);
				flux_i_density_energy += factor*(ff_flux_contribution_density_energy[0].x+flux_contribution_i_density_energy.x);
				flux_i_momentum.x += factor*(ff_flux_contribution_momentum_x[0].x + flux_contribution_i_momentum_x.x);
				flux_i_momentum.y += factor*(ff_flux_contribution_momentum_y[0].x + flux_contribution_i_momentum_y.x);
				flux_i_momentum.z += factor*(ff_flux_contribution_momentum_z[0].x + flux_contribution_i_momentum_z.x);
			
				factor = float(0.5f)*normal.y;
				flux_i_density += factor*(ff_variable[VAR_MOMENTUM+1]+momentum_i.y);
				flux_i_density_energy += factor*(ff_flux_contribution_density_energy[0].y+flux_contribution_i_density_energy.y);
				flux_i_momentum.x += factor*(ff_flux_contribution_momentum_x[0].y + flux_contribution_i_momentum_x.y);
				flux_i_momentum.y += factor*(ff_flux_contribution_momentum_y[0].y + flux_contribution_i_momentum_y.y);
				flux_i_momentum.z += factor*(ff_flux_contribution_momentum_z[0].y + flux_contribution_i_momentum_z.y);

				factor = float(0.5f)*normal.z;
				flux_i_density += factor*(ff_variable[VAR_MOMENTUM+2]+momentum_i.z);
				flux_i_density_energy += factor*(ff_flux_contribution_density_energy[0].z+flux_contribution_i_density_energy.z);
				flux_i_momentum.x += factor*(ff_flux_contribution_momentum_x[0].z + flux_contribution_i_momentum_x.z);
				flux_i_momentum.y += factor*(ff_flux_contribution_momentum_y[0].z + flux_contribution_i_momentum_y.z);
				flux_i_momentum.z += factor*(ff_flux_contribution_momentum_z[0].z + flux_contribution_i_momentum_z.z);

			}
		}

		fluxes[i + VAR_DENSITY*nelr] = flux_i_density;
		fluxes[i + (VAR_MOMENTUM+0)*nelr] = flux_i_momentum.x;
		fluxes[i + (VAR_MOMENTUM+1)*nelr] = flux_i_momentum.y;
		fluxes[i + (VAR_MOMENTUM+2)*nelr] = flux_i_momentum.z;
		fluxes[i + VAR_DENSITY_ENERGY*nelr] = flux_i_density_energy;
	} else {
		heartwall_kernel();
	}
}

#define N 2

__global__ void Kernel1(int nelr, int* elements_surrounding_elements, float* normals, float* variables, float* fluxes)
{
    cfdkernel(nelr, elements_surrounding_elements, normals, variables, fluxes);
}

__global__ void Kernel2()
{
    heartwall_kernel();
}
void run_kernels(int nelr, int* elements_surrounding_elements, float* normals, float* variables, float* fluxes)
{
	dim3 Dg(nelr / block_length), Db(block_length);
	int numBlocks = nelr/block_length;
	if(numBlocks > ALL_POINTS){ 
		numBlocks = ALL_POINTS;
	}
    cudaStream_t streams[N];
    for(int i = 0; i < N; i++)
        cudaStreamCreate(&streams[i]);
	//kernel<<<numBlocks,Db>>>(nelr, elements_surrounding_elements, normals, variables, fluxes);
    Kernel1<<<numBlocks/3,Db, 0, streams[0]>>>(nelr, elements_surrounding_elements, normals, variables, fluxes);
    Kernel2<<<numBlocks-numBlocks/3,Db, 0, streams[1]>>>();
	getLastCudaError("kernel failed");
}

/*
 * Main function
 */
int main(int argc, char** argv)
{
	//cudaDeviceSetCacheConfig(cudaFuncCachePreferL1);
	
	if (argc < 2)
	{
		std::cout << "specify data file name" << std::endl;
		return 0;
	}
	const char* data_file_name = argv[1];
	
	cudaDeviceProp prop;
	int dev;
	
	checkCudaErrors(cudaSetDevice(0));
	checkCudaErrors(cudaGetDevice(&dev));
	checkCudaErrors(cudaGetDeviceProperties(&prop, dev));
	
	printf("Name:                     %s\n", prop.name);

	// set far field conditions and load them into constant memory on the gpu
	{
		float h_ff_variable[NVAR];
		const float angle_of_attack = float(3.1415926535897931 / 180.0f) * float(deg_angle_of_attack);
		
		h_ff_variable[VAR_DENSITY] = float(1.4);
		
		float ff_pressure = float(1.0f);
		float ff_speed_of_sound = sqrt(GAMMA*ff_pressure / h_ff_variable[VAR_DENSITY]);
		float ff_speed = float(ff_mach)*ff_speed_of_sound;
		
		float3 ff_velocity;
		ff_velocity.x = ff_speed*float(cos((float)angle_of_attack));
		ff_velocity.y = ff_speed*float(sin((float)angle_of_attack));
		ff_velocity.z = 0.0f;
		
		h_ff_variable[VAR_MOMENTUM+0] = h_ff_variable[VAR_DENSITY] * ff_velocity.x;
		h_ff_variable[VAR_MOMENTUM+1] = h_ff_variable[VAR_DENSITY] * ff_velocity.y;
		h_ff_variable[VAR_MOMENTUM+2] = h_ff_variable[VAR_DENSITY] * ff_velocity.z;
				
		h_ff_variable[VAR_DENSITY_ENERGY] = h_ff_variable[VAR_DENSITY]*(float(0.5f)*(ff_speed*ff_speed)) + (ff_pressure / float(GAMMA-1.0f));

		float3 h_ff_momentum;
		h_ff_momentum.x = *(h_ff_variable+VAR_MOMENTUM+0);
		h_ff_momentum.y = *(h_ff_variable+VAR_MOMENTUM+1);
		h_ff_momentum.z = *(h_ff_variable+VAR_MOMENTUM+2);
		float3 h_ff_flux_contribution_momentum_x;
		float3 h_ff_flux_contribution_momentum_y;
		float3 h_ff_flux_contribution_momentum_z;
		float3 h_ff_flux_contribution_density_energy;
		compute_flux_contribution(h_ff_variable[VAR_DENSITY], h_ff_momentum, h_ff_variable[VAR_DENSITY_ENERGY], ff_pressure, ff_velocity, h_ff_flux_contribution_momentum_x, h_ff_flux_contribution_momentum_y, h_ff_flux_contribution_momentum_z, h_ff_flux_contribution_density_energy);

		// copy far field conditions to the gpu
		checkCudaErrors( cudaMemcpyToSymbol(ff_variable,          h_ff_variable,          NVAR*sizeof(float)) );
		checkCudaErrors( cudaMemcpyToSymbol(ff_flux_contribution_momentum_x, &h_ff_flux_contribution_momentum_x, sizeof(float3)) );
		checkCudaErrors( cudaMemcpyToSymbol(ff_flux_contribution_momentum_y, &h_ff_flux_contribution_momentum_y, sizeof(float3)) );
		checkCudaErrors( cudaMemcpyToSymbol(ff_flux_contribution_momentum_z, &h_ff_flux_contribution_momentum_z, sizeof(float3)) );
		
		checkCudaErrors( cudaMemcpyToSymbol(ff_flux_contribution_density_energy, &h_ff_flux_contribution_density_energy, sizeof(float3)) );		
	}
	int nel;
	int nelr;
	
	// read in domain geometry
	float* areas;
	int* elements_surrounding_elements;
	float* normals;
	{
		std::ifstream file(data_file_name);
	
		file >> nel;
		nelr = block_length*((nel / block_length )+ std::min(1, nel % block_length));

		float* h_areas = new float[nelr];
		int* h_elements_surrounding_elements = new int[nelr*NNB];
		float* h_normals = new float[nelr*NDIM*NNB];

				
		// read in data
		for(int i = 0; i < nel; i++)
		{
			file >> h_areas[i];
			for(int j = 0; j < NNB; j++)
			{
				file >> h_elements_surrounding_elements[i + j*nelr];
				if(h_elements_surrounding_elements[i+j*nelr] < 0) h_elements_surrounding_elements[i+j*nelr] = -1;
				h_elements_surrounding_elements[i + j*nelr]--; //it's coming in with Fortran numbering				
				
				for(int k = 0; k < NDIM; k++)
				{
					file >> h_normals[i + (j + k*NNB)*nelr];
					h_normals[i + (j + k*NNB)*nelr] = -h_normals[i + (j + k*NNB)*nelr];
				}
			}
		}
		
		// fill in remaining data
		int last = nel-1;
		for(int i = nel; i < nelr; i++)
		{
			h_areas[i] = h_areas[last];
			for(int j = 0; j < NNB; j++)
			{
				// duplicate the last element
				h_elements_surrounding_elements[i + j*nelr] = h_elements_surrounding_elements[last + j*nelr];	
				for(int k = 0; k < NDIM; k++) h_normals[last + (j + k*NNB)*nelr] = h_normals[last + (j + k*NNB)*nelr];
			}
		}
		
		areas = alloc<float>(nelr);
		upload<float>(areas, h_areas, nelr);

		elements_surrounding_elements = alloc<int>(nelr*NNB);
		upload<int>(elements_surrounding_elements, h_elements_surrounding_elements, nelr*NNB);

		normals = alloc<float>(nelr*NDIM*NNB);
		upload<float>(normals, h_normals, nelr*NDIM*NNB);
				
		delete[] h_areas;
		delete[] h_elements_surrounding_elements;
		delete[] h_normals;
	}

	// Create arrays and set initial conditions
	float* variables = alloc<float>(nelr*NVAR);
	initialize_variables(nelr, variables);

	float* old_variables = alloc<float>(nelr*NVAR);   	
	float* fluxes = alloc<float>(nelr*NVAR);
	float* step_factors = alloc<float>(nelr); 

	// make sure all memory is floatly allocated before we start timing
	initialize_variables(nelr, old_variables);
	initialize_variables(nelr, fluxes);
	cudaMemset( (void*) step_factors, 0, sizeof(float)*nelr );
	// make sure CUDA isn't still doing something before we start timing
	cudaThreadSynchronize();

	// these need to be computed the first time in order to compute time step
	std::cout << "Starting..." << std::endl;

	StopWatchInterface *timer = 0;
	  //	unsigned int timer = 0;
	
	
	heartmain(argc-1, argv+1);
	
	// CUT_SAFE_CALL( cutCreateTimer( &timer));
	// CUT_SAFE_CALL( cutStartTimer( timer));
	sdkCreateTimer(&timer); 
	sdkStartTimer(&timer);
	run_kernels(nelr, elements_surrounding_elements, normals, variables, fluxes);
	cudaThreadSynchronize();
	//	CUT_SAFE_CALL( cutStopTimer(timer) );  
	sdkStopTimer(&timer); 

	std::cout  << "runtime: " << sdkGetAverageTimerValue(&timer) << std::endl;

	
	std::cout << "Cleaning up..." << std::endl;
	dealloc<float>(areas);
	dealloc<int>(elements_surrounding_elements);
	dealloc<float>(normals);
	
	dealloc<float>(variables);
	dealloc<float>(old_variables);
	dealloc<float>(fluxes);
	dealloc<float>(step_factors);

	std::cout << "Done..." << std::endl;

	return 0;
}
