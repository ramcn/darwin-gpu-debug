#include "dsoft_hash_table.hpp"
#include <iostream>
#include <stdlib.h>
using namespace std;

int DsoftHashTable::GetKmerSize() {
    return kmer_size_;
}

bool DsoftHashTable::IsPresent(uint32_t index) {
    uint32_t start_index = (index > 0) ? index_table_[index-1] : 0;
    uint32_t end_index = index_table_[index];
    return (end_index - start_index <= kmer_max_occurence_);
}

DsoftHashTable::DsoftHashTable(char* ref_str, uint32_t ref_length, std::string shape, uint32_t seed_occurence_multiple, uint32_t bin_size) {
    L2Hashtable l2hashtable;

    hashtable.read_from_file("dsoft_ht");
    l2hashtable.l2_read_from_file("dsoft_ht");
    hashtable.set_l2_hashtable(&l2hashtable);

    std::cerr <<"dsoft hashtable opened successfully" << endl;

    shape_size_ = shape.length(); 
    int kmer_size = 0;
    for (int i = 0; i < shape_size_; i++) {
        kmer_size += (shape[i] == '1');
    }

    assert(kmer_size <= 15);
    assert(kmer_size > 3); 

    kmer_size_ = kmer_size;
    bin_size_  = bin_size;
    log_bin_size_ = (uint32_t) (log2(bin_size_));
    ref_size_ = ref_length;
    kmer_max_occurence_ = seed_occurence_multiple * (1+((ref_length) >> (2*kmer_size)));

    GenerateShapePos(shape);

    num_bins_ = 1 + (ref_size_ >> log_bin_size_);
}

DsoftHashTable::~DsoftHashTable() {
}


int DsoftHashTable::DSOFT(char* query, uint32_t query_length, int N, int threshold, uint64_t* candidate_hit_offset, uint64_t* bin_count_offset_array, uint32_t* nz_bins_array, int max_candidates) {
    uint32_t index, offset;
    uint32_t hit, bin, curr_count, last_offset, new_count, adjacent_bin, adjacent_count;
    uint64_t hit_offset;
    uint64_t num_nz_bins = 0;
    int num_seeds = 0;
    int num_candidates = 0;

    for (int i = 0; i < query_length - shape_size_; i++) {
        //index = GetKmerIndexAtPos(query, i);
        index = hashtable.get_hash(&query[i]);
        if (index != (1 << 31)) { // is this still valid [?]
            offset = i;
            //uint32_t start_index = (index > 0) ? index_table_[index-1] : 0;
            //uint32_t end_index = index_table_[index];
	    uint32_t frequency = hashtable.get_frequency(index);
	    uint32_t start_index = hashtable.get_offset(index);
	    uint32_t end_index = start_index + frequency;

            if (end_index - start_index <= kmer_max_occurence_) {
                num_seeds++;
                //for (uint32_t j = start_index; j < end_index; j++) {
                for (uint32_t j = 0; j < frequency; j++) {
                    //hit = pos_table_[j];
		    hit = hashtable.get_location(start_index+j);
                    if (hit >= offset) {
                        bin = ((hit - offset) / bin_size_);
                        curr_count = (bin_count_offset_array[bin] >> 32);
                        last_offset = ((bin_count_offset_array[bin] << 32) >> 32);
                        if (curr_count < threshold) {
                            new_count = ((offset - last_offset > kmer_size_) || (curr_count == 0)) ? curr_count + kmer_size_ : curr_count + (offset - last_offset);
                            bin_count_offset_array[bin] = ((uint64_t) new_count << 32) + offset; 
                            if (new_count >= threshold) {
                                uint64_t hit_offset = ((uint64_t) hit << 32) + offset;
                                if (num_candidates >= max_candidates) {
                                    break;
                                }
                                candidate_hit_offset[num_candidates++] = hit_offset;
                            }
                            if ((curr_count == 0) && (num_nz_bins < nz_bins)) {
                                nz_bins_array[num_nz_bins] = bin;
                                num_nz_bins++;
                            }
                        }
                    }
                }
                if (num_seeds >= N) {
                    break;
                }
            }
        }
    }
    for (uint32_t i = 0; i < num_nz_bins; i++) {
        bin = nz_bins_array[i];
        bin_count_offset_array[bin] = 0;
    }
    return num_candidates;
}


