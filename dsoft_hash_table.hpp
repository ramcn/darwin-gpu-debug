#include <string>
#include <vector>
#include <map>
#include <assert.h>
#include <stdint.h>
#include <math.h>
#include "ntcoding.h"
#include "l1_hashtable.hpp"
#include "l2_hashtable.hpp"

#define nz_bins 25000000

class DsoftHashTable {
    private:
        uint32_t index_table_size_;
        uint32_t ref_size_;
        uint32_t bin_size_;
        uint32_t log_bin_size_;
        int kmer_size_;
        int shape_size_;
        uint32_t kmer_max_occurence_;
        uint32_t *index_table_;
        uint32_t *pos_table_;

        Hashtable hashtable;

        uint32_t num_bins_;
        uint64_t* bin_count_offset_array_;
        uint32_t* nz_bins_;


    public:
        DsoftHashTable();
        DsoftHashTable(char* ref_str, uint32_t ref_length, std::string shape, uint32_t seed_occurence_multiple, uint32_t bin_size);
        ~DsoftHashTable();

        bool IsPresent(uint32_t index);
        int GetKmerSize();
//        void DSOFT(uint64_t* seed_offset_vector, int N, int threshold, std::vector<uint64_t>& candidate_hit_offset);
//        int DSOFT(char* query, uint32_t query_length, int N, int threshold, uint64_t* candidate_hit_offset, int max_candidates);
        int DSOFT(char* query, uint32_t query_length, int N, int threshold, uint64_t* candidate_hit_offset, uint64_t* bin_count_offset_array, uint32_t* nz_bins_array, int max_candidates);
};



