/*
MIT License

Copyright (c) 2018 Yatish Turakhia, Gill Bejerano and William Dally

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
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

#include <iostream>
#include <queue>
#include <limits>
#include "ntcoding.h"

#define STRING_BUFFER_LEN 1024
#define PRECOMPILED_BINARY "align"
#define NUM_DEVICES 1
#define MAX_NUM_DEVICES 16
#define MAX_TILE_SIZE 512 
#define BLOCK_WIDTH 512
#define GPU


extern int forward_kernel_counter;
extern int reverse_kernel_counter;

std::queue<int> AlignWithBT(char* ref_seq, long long int ref_len, char* query_seq, long long int query_len, int* sub_mat, int gap_open, int gap_extend, int query_pos, int ref_pos, bool reverse, bool first, int early_terminate);

