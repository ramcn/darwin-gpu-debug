/*
MIT License

Copyright (c) 2018 Yatish Turakhia, Gill Bejerano and William Dally

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without ion, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

#include "gact.h"
#include <math.h>


enum states {Z, D, I, M};
enum AlnOperands {ZERO_OP, DELETE_OP, INSERT_OP, MATCH_OP};


int dequeue(int *BT_states, int *front, int *queuesize) {
   int data = BT_states[(*front)++];
	
   if(*front == MAX_TILE_SIZE) {
      *front = 0;
   }
	
   (*queuesize)--;
   return data;  
}



#ifdef GPU
std::queue<int> GpuXL(char* ref_str, long long int ref_tile_length, char* query_str, long long int query_tile_length, int* sub_mat, int gap_open, int gap_extend, int ref_pos, int query_pos, bool reverse, bool first_tile, int early_terminate, int* BT_states, int *queuesize, char strand, int *rear, int *front, int mode); 
#endif

Alignment GACT (char* ref_str, char* query_str, std::string ref_name, std::string query_name, int* sub_mat, int gap_open, int gap_extend, int tile_size, int tile_overlap, int ref_pos, int query_pos, uint32_t ref_length, uint32_t query_length, char strand, int first_tile_score_threshold, int mode, int thread_id) {
    std::queue<int> BT_states_std;

    std::string aligned_ref_str = "";
    std::string aligned_query_str = "";

    Alignment alignment;
    alignment.ref_name = ref_name;
    alignment.query_name = query_name;
    alignment.aligned_ref_str = "";
    alignment.aligned_query_str = "";
    alignment.ref_start = ref_pos;
    alignment.query_start = query_pos;
    alignment.aligned_ref_len = 0;
    alignment.aligned_query_len = 0;
    alignment.ref_len = ref_length;
    alignment.query_len = query_length;
    alignment.strand = strand;
    alignment.score = 0;
    alignment.flag = 0;
    
    int ref_tile_length = tile_size;
    int query_tile_length = tile_size;
    int rev_ref_pos = ref_pos;
    int rev_query_pos = query_pos;
   
    int max_ref_pos = 0;
    int max_query_pos = 0;
    int i = 0;
    int j = 0;
    
    int first_tile_score = 0;
    bool first_tile = true;

    
    while ((ref_pos > 0) && (query_pos > 0) && (((i > 0) && (j > 0)) || first_tile)) {
        //change the tile length if elements less than that of the tile size
        ref_tile_length = (ref_pos > tile_size) ? tile_size : ref_pos;
        query_tile_length = (query_pos > tile_size) ? tile_size : query_pos;

    	int *BT_states = (int *) malloc(sizeof(int)*MAX_TILE_SIZE);
    	int *front = (int *) malloc(sizeof(int));
    	int *rear = (int *) malloc(sizeof(int));
    	int *queuesize = (int *) malloc(sizeof(int));
    	*front = 0; *rear=-1; *queuesize=0; 

        if(mode==1) { 
        	BT_states_std = AlignWithBT (ref_str+(ref_pos-ref_tile_length), ref_tile_length, 
				query_str+(query_pos-query_tile_length), query_tile_length, 
				sub_mat, gap_open, gap_extend, query_tile_length, ref_tile_length, false, 
				first_tile, (tile_size - tile_overlap));
        }
#ifdef GPU
        else {
       	       BT_states_std = GpuXL (ref_str+(ref_pos-ref_tile_length), ref_tile_length, 
				query_str+(query_pos-query_tile_length), query_tile_length, 
				sub_mat, gap_open, gap_extend, ref_tile_length, query_tile_length, false, 
				first_tile, (tile_size - tile_overlap), BT_states, queuesize, strand, rear, front, mode);
	}
#endif

        i = 0;
        j = 0;
        int tile_score = BT_states_std.front();
        BT_states_std.pop();

        
        if (first_tile) {
            ref_pos = ref_pos - ref_tile_length + BT_states_std.front();
            max_ref_pos = BT_states_std.front(); BT_states_std.pop();
            query_pos = query_pos - query_tile_length + BT_states_std.front();
            max_query_pos = BT_states_std.front(); BT_states_std.pop();
            rev_ref_pos = ref_pos;
            rev_query_pos = query_pos;
            first_tile_score = tile_score;
        }

        int num_tb = BT_states_std.size();
        char* ref_buf = (char*) malloc(num_tb);
        char* query_buf = (char*) malloc(num_tb);
        int ref_buf_curr = num_tb-1;
        int query_buf_curr = num_tb-1;

        while (!BT_states_std.empty()) {
            first_tile = false;
            int state = BT_states_std.front(); BT_states_std.pop();
            if (state == M) {
                ref_buf[ref_buf_curr--] = ref_str[ref_pos - j - 1];
                query_buf[query_buf_curr--] = query_str[query_pos - i - 1];
                i += 1;
                j += 1;
            }
            if (state == I) {
                ref_buf[ref_buf_curr--] = ref_str[ref_pos - j - 1];
                query_buf[query_buf_curr--] = '-';
                j += 1;
            }
            if (state == D) {
                ref_buf[ref_buf_curr--] = '-';
                query_buf[query_buf_curr--] = query_str[query_pos - i - 1];
                i += 1;
            }
        }
        if (num_tb > 0) {
            aligned_ref_str = std::string(ref_buf, num_tb) + aligned_ref_str;
            aligned_query_str = std::string(query_buf, num_tb) + aligned_query_str;
        }

    	free(BT_states);
    	free(front);
    	free(rear);
    	free(queuesize);
        free(ref_buf);
        free(query_buf);
        ref_pos -= (j);
        query_pos -= (i);
        alignment.aligned_ref_len += j;
        alignment.aligned_query_len += i;
    }
    
    alignment.ref_start = ref_pos;
    alignment.query_start = query_pos;

    ref_pos = rev_ref_pos;
    query_pos = rev_query_pos;
    
    i =  tile_size;
    j = tile_size;
    
    //starts with the first tile
    while ((ref_pos < ref_length) && (query_pos < query_length) && (((i > 0) && (j > 0)) || first_tile)) {
        ref_tile_length = (ref_pos + tile_size < ref_length) ? tile_size : ref_length - ref_pos;
        query_tile_length = (query_pos + tile_size < query_length) ? tile_size : query_length - query_pos;
    	int *BT_states = (int *) malloc(sizeof(int)*MAX_TILE_SIZE);
    	int *front = (int *) malloc(sizeof(int));
    	int *rear = (int *) malloc(sizeof(int));
    	int *queuesize = (int *) malloc(sizeof(int));
    	*front = 0; *rear=-1; *queuesize=0; 

        if(mode == 1) { 
        	BT_states_std = AlignWithBT (ref_str+ref_pos, ref_tile_length, query_str+query_pos, query_tile_length, 
				sub_mat, gap_open, gap_extend, query_tile_length, ref_tile_length, true, 
				first_tile, (tile_size - tile_overlap));
	}
#ifdef GPU
        else {
                BT_states_std = GpuXL (ref_str+ref_pos, ref_tile_length, query_str+query_pos, query_tile_length, 
				sub_mat, gap_open, gap_extend, ref_tile_length, query_tile_length, true, 
				first_tile, (tile_size - tile_overlap), BT_states, queuesize, strand, rear, front, mode);
	}
#endif
        i = 0;
        j = 0;
        int tile_score = BT_states_std.front(); BT_states_std.pop(); 

        if (first_tile) {
            ref_pos = ref_pos + ref_tile_length - BT_states_std.front();
            max_ref_pos = BT_states_std.front(); BT_states_std.pop();
            query_pos = query_pos + query_tile_length - BT_states_std.front();
            max_query_pos = BT_states_std.front(); BT_states_std.pop();
        }

        int num_tb = BT_states_std.size();;
        char* ref_buf = (char*) malloc(num_tb);
        char* query_buf = (char*) malloc(num_tb);
        int ref_buf_curr = 0;
        int query_buf_curr = 0;
        
        while (!BT_states_std.empty()) {
            first_tile = false;
            int state = BT_states_std.front(); BT_states_std.pop();
            if (state == M) {
                ref_buf[ref_buf_curr++] = ref_str[ref_pos + j];
                query_buf[query_buf_curr++] = query_str[query_pos + i];
                i += 1;
                j += 1;
            }
            if (state == I) {
                ref_buf[ref_buf_curr++] = ref_str[ref_pos + j];
                query_buf[query_buf_curr++] = '-';
                j += 1;
            }
            if (state == D) {
                ref_buf[ref_buf_curr++] = '-';
                query_buf[query_buf_curr++] = query_str[query_pos + i];
                i += 1;
            }
        }
        if (num_tb > 0) {
            aligned_ref_str = aligned_ref_str + std::string(ref_buf, num_tb);
            aligned_query_str = aligned_query_str + std::string(query_buf, num_tb);
        }

    	free(BT_states);
    	free(front);
    	free(rear);
    	free(queuesize);
        free(ref_buf);
        free(query_buf);
        ref_pos += (j);
        query_pos += (i);
        alignment.aligned_ref_len += j;
        alignment.aligned_query_len += i;
    }

    int total_score = 0;
    bool open = true;
    for (uint32_t j = 0; j < aligned_ref_str.length(); j++) {
        char ref_nt = aligned_ref_str[j];
        char query_nt = aligned_query_str[j];
        if (ref_nt == '-' || query_nt == '-') {
            total_score += (open) ? gap_open : gap_extend;
            open = false;
        }
        else {
            total_score += sub_mat[5*NtChar2Int(query_nt) + NtChar2Int(ref_nt)];
            open = true;
        }
    }
    alignment.aligned_ref_str = aligned_ref_str;
    alignment.aligned_query_str = aligned_query_str;
    alignment.score = total_score;
    return alignment;
}


#ifdef GPU

#define cudaErrorCheck(err) \
        if(err != cudaSuccess) {\
                printf("error in cuda call at line %d\n", __LINE__);\
        }

void xl_cpu (const char *  a, const  int m,
                                        const char *  b, const int n,
                                        const int gap_open, const int gap_extend, const int jj,
                                        const int query_pos, const int ref_pos, const int reverse, const int first, const int early_terminate,
                                        int *  sub_mat, int *  dir_matrix_arg,   int *  max_score_arg,
                                        int *  max_i_arg,  int *  max_j_arg,  int *  pos_score_arg);
                                        //int *h_matrix_rd, int *m_matrix_rd, int *i_matrix_rd, int *d_matrix_rd,
                                        //int *h_matrix_wr, int *m_matrix_wr, int *i_matrix_wr, int *d_matrix_wr);

__global__ void xl (const char *  a, const  int m,
                                        const char *  b, const int n,
                                        const int gap_open, const int gap_extend, const int jj,
                                        const int query_pos, const int ref_pos, const int reverse, const int first, const int early_terminate,
                                        int *  sub_mat, int *  dir_matrix_arg,   int *  max_score_arg,
                                        int *  max_i_arg,  int *  max_j_arg,  int *  pos_score_arg);
                                        //int *h_matrix_rd, int *m_matrix_rd, int *i_matrix_rd, int *d_matrix_rd,
                                        //int *h_matrix_wr, int *m_matrix_wr, int *i_matrix_wr, int *d_matrix_wr);



std::queue<int> GpuXL(char* a, long long int m, char* b, long long int n, int *sub_mat, int gap_open, int gap_extend, int ref_pos, int query_pos, bool reverse, bool first, int early_terminate, int* BT_states_arg, int *queuesize, char strand, int *rear, int *front, int mode) {

    int score;
    int nbb=1, jj=0;
    int INF = 2147483647;


    int dir_matrix[MAX_TILE_SIZE+1][MAX_TILE_SIZE+1];
    int *max_score = (int *)malloc(sizeof(int));
    int *max_i = (int *)malloc(sizeof(int));
    int *max_j = (int *) malloc(sizeof(int));
    int *pos_score = (int *)malloc(sizeof(int));

    for (int i = 0; i < m+1; i++)
           dir_matrix[i][0] = ZERO_OP;
    for (int j = 0; j < n + 1; j++)
           dir_matrix[0][j] = ZERO_OP;


        int h_matrix_wr[MAX_TILE_SIZE + 1], m_matrix_wr[MAX_TILE_SIZE + 1], i_matrix_wr[MAX_TILE_SIZE + 1], d_matrix_wr[MAX_TILE_SIZE + 1];
        int h_matrix_rd[MAX_TILE_SIZE + 1], m_matrix_rd[MAX_TILE_SIZE + 1], i_matrix_rd[MAX_TILE_SIZE + 1], d_matrix_rd[MAX_TILE_SIZE + 1];

        // MICRO
        for (int i = 0; i < BLOCK_WIDTH + 1; i++) {
                h_matrix_rd[i] = 0; m_matrix_rd[i] = 0; i_matrix_rd[i] = -INF; d_matrix_rd[i] = -INF;
                h_matrix_wr[i] = 0; m_matrix_wr[i] = 0; i_matrix_wr[i] = -INF; d_matrix_wr[i] = -INF;
        }

    char *cl_a, *cl_b;
    int *cl_sub_mat, *cl_dir_matrix;
    int *cl_max_score, *cl_max_i, *cl_max_j, *cl_pos_score;   
    int *cl_h_wr, *cl_m_wr, *cl_i_wr, *cl_d_wr;
    int *cl_h_rd, *cl_m_rd, *cl_i_rd, *cl_d_rd;
    cudaError_t err;
    cudaErrorCheck(cudaMalloc((void **) &cl_a,  m*sizeof(char)));
    cudaErrorCheck(cudaMalloc((void **) &cl_b,  n*sizeof(char)));
    cudaErrorCheck(cudaMalloc((void **) &cl_sub_mat,   sizeof(int)*25));
    cudaErrorCheck(cudaMalloc((void **) &cl_dir_matrix,  sizeof(int)*(MAX_TILE_SIZE+1)*(MAX_TILE_SIZE+1)));
    cudaErrorCheck(cudaMalloc((void **) &cl_max_score,   sizeof(int)));
    cudaErrorCheck(cudaMalloc((void **) &cl_max_i,   sizeof(int)));
    cudaErrorCheck(cudaMalloc((void **) &cl_max_j,   sizeof(int)));
    cudaErrorCheck(cudaMalloc((void **) &cl_pos_score,   sizeof(int)));


    cudaErrorCheck(cudaMalloc((void **) &cl_h_rd, sizeof(int) * (MAX_TILE_SIZE + 1)));
    cudaErrorCheck(cudaMalloc((void **) &cl_m_rd, sizeof(int) * (MAX_TILE_SIZE + 1)));
    cudaErrorCheck(cudaMalloc((void **) &cl_i_rd, sizeof(int) * (MAX_TILE_SIZE + 1)));
    cudaErrorCheck(cudaMalloc((void **) &cl_d_rd, sizeof(int) * (MAX_TILE_SIZE + 1)));
    cudaErrorCheck(cudaMalloc((void **) &cl_h_wr, sizeof(int) * (MAX_TILE_SIZE + 1)));
    cudaErrorCheck(cudaMalloc((void **) &cl_m_wr, sizeof(int) * (MAX_TILE_SIZE + 1)));
    cudaErrorCheck(cudaMalloc((void **) &cl_i_wr, sizeof(int) * (MAX_TILE_SIZE + 1)));
    cudaErrorCheck(cudaMalloc((void **) &cl_d_wr, sizeof(int) * (MAX_TILE_SIZE + 1)));


    cudaErrorCheck(cudaMemcpy(cl_a, a, m*sizeof(char), cudaMemcpyHostToDevice));
    cudaErrorCheck(cudaMemcpy(cl_b, b, n*sizeof(char), cudaMemcpyHostToDevice));
    cudaErrorCheck(cudaMemcpy(cl_sub_mat, sub_mat, sizeof(int)*25, cudaMemcpyHostToDevice));
    cudaErrorCheck(cudaMemcpy(cl_dir_matrix, dir_matrix, sizeof(int)*(MAX_TILE_SIZE+1)*(MAX_TILE_SIZE+1), cudaMemcpyHostToDevice));
    cudaErrorCheck(cudaMemcpy(cl_max_score, max_score, sizeof(int), cudaMemcpyHostToDevice));
    cudaErrorCheck(cudaMemcpy(cl_max_i, max_i, sizeof(int), cudaMemcpyHostToDevice));
    cudaErrorCheck(cudaMemcpy(cl_max_j, max_j, sizeof(int), cudaMemcpyHostToDevice));
    cudaErrorCheck(cudaMemcpy(cl_pos_score, pos_score, sizeof(int), cudaMemcpyHostToDevice));

    cudaErrorCheck(cudaMemcpy(cl_h_rd, h_matrix_rd, sizeof(int) * (MAX_TILE_SIZE+1), cudaMemcpyHostToDevice));
    cudaErrorCheck(cudaMemcpy(cl_m_rd, m_matrix_rd, sizeof(int) * (MAX_TILE_SIZE+1), cudaMemcpyHostToDevice));
    cudaErrorCheck(cudaMemcpy(cl_i_rd, i_matrix_rd, sizeof(int) * (MAX_TILE_SIZE+1), cudaMemcpyHostToDevice));
    cudaErrorCheck(cudaMemcpy(cl_d_rd, d_matrix_rd, sizeof(int) * (MAX_TILE_SIZE+1), cudaMemcpyHostToDevice));
    cudaErrorCheck(cudaMemcpy(cl_h_wr, h_matrix_wr, sizeof(int) * (MAX_TILE_SIZE+1), cudaMemcpyHostToDevice));
    cudaErrorCheck(cudaMemcpy(cl_m_wr, m_matrix_wr, sizeof(int) * (MAX_TILE_SIZE+1), cudaMemcpyHostToDevice));
    cudaErrorCheck(cudaMemcpy(cl_i_wr, i_matrix_wr, sizeof(int) * (MAX_TILE_SIZE+1), cudaMemcpyHostToDevice));
    cudaErrorCheck(cudaMemcpy(cl_d_wr, d_matrix_wr, sizeof(int) * (MAX_TILE_SIZE+1), cudaMemcpyHostToDevice));


    if(mode == 2) {
    	xl<<<1,512>>>(cl_a, m, cl_b, n, gap_open, gap_extend, jj, query_pos, ref_pos, reverse, first, early_terminate,
					cl_sub_mat, cl_dir_matrix, cl_max_score, cl_max_i, cl_max_j, cl_pos_score);
					//cl_h_rd, cl_m_rd, cl_i_rd, cl_d_rd,
					//cl_h_wr, cl_m_wr, cl_i_wr, cl_d_wr);

    	cudaDeviceSynchronize();

    	cudaMemcpy(dir_matrix, cl_dir_matrix, sizeof(int)*(MAX_TILE_SIZE+1)*(MAX_TILE_SIZE+1), cudaMemcpyDeviceToHost);
    	cudaMemcpy(max_score, cl_max_score, sizeof(int), cudaMemcpyDeviceToHost);
    	cudaMemcpy(max_i, cl_max_i, sizeof(int), cudaMemcpyDeviceToHost);
    	cudaMemcpy(max_j, cl_max_j, sizeof(int), cudaMemcpyDeviceToHost);
    	cudaMemcpy(pos_score, cl_pos_score, sizeof(int), cudaMemcpyDeviceToHost);

    } else if(mode ==3) {
    	xl_cpu(a, m, b, n, gap_open, gap_extend, jj, query_pos, ref_pos, reverse, first, early_terminate,
                                        sub_mat, (int *)dir_matrix, max_score, max_i, max_j, pos_score);
                                        //h_matrix_rd, m_matrix_rd, i_matrix_rd,d_matrix_rd,
                                        //h_matrix_wr, m_matrix_wr, i_matrix_wr, d_matrix_wr);
    }

    err = cudaFree(cl_a);
    err = cudaFree(cl_b);
    err = cudaFree(cl_sub_mat);
    err = cudaFree(cl_dir_matrix);
    err = cudaFree(cl_max_score);
    err = cudaFree(cl_max_i);
    err = cudaFree(cl_max_j);
    err = cudaFree(cl_pos_score);


  std::queue<int> BT_states;

  int i_curr=ref_pos, j_curr=query_pos;
  int i_steps = 0, j_steps = 0;

  int open = 0;
  if (first) {
      i_curr = *max_i;
      j_curr = *max_j;
      BT_states.push(*max_score);
      BT_states.push(i_curr);
      BT_states.push(j_curr);
  }
  else {
      BT_states.push(*pos_score);
  }

  int state = dir_matrix[i_curr][j_curr] % 4;

  while (state != Z) {
      if ((i_steps >= early_terminate) || (j_steps >= early_terminate)) {
          break;
      }
      BT_states.push(state);
      if (state == M) {
          state = (dir_matrix[i_curr-1][j_curr-1] % 4);
          i_curr--;
          j_curr--;
          i_steps++;
          j_steps++;
      }
      else if (state == I) {
          state = (dir_matrix[i_curr][j_curr] & (2 << INSERT_OP)) ? M : I;
          i_curr--;
          i_steps++;
      }
      else if (state == D) {
          state = (dir_matrix[i_curr][j_curr] & (2 << DELETE_OP)) ? M : D;
          j_curr--;
          j_steps++;
      }
  };

  free(max_score);
  free(max_i);
  free(max_j);
  free(pos_score);
  return BT_states;

}


#endif
