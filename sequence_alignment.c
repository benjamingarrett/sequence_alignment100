#include<limits.h>
#include<stdint.h>
#include<stdio.h>
#include<string.h>
#include"../memoization/memo_long_int.h"
#include"../lcs_instance_reader/lcs_instance_reader.h"
#include"sequence_alignment.h"

//#define PRINT_REFERENCES  

const char * sequence_align_instance_fname_parameter = "--sequence_align_instance_fname";

int64_t * X, * Y;
int32_t nX, nY;
int64_t memo_calls, dist_calls, cache_misses, prev_cache_misses;
long int cost;
long int pattern_references;
long int cache_reads;
long int cache_writes;
long int cache_misses;
long int pattern_reference_threshold;
long int memo_calls, dist_calls;
long int lru_queue_size;

int64_t alignment_cost(long int len_s, long int len_t);

char (*check_preemptive_halt_sequence_align_)();
char check_preemptive_halt_sequence_align_activated();
char check_preemptive_halt_sequence_align_deactivated();

void align();

/* PUBLIC */

void initialize_seq(int argc, char **argv){
  //printf("initialize edit distance\n");
  char instance_fname[200];
  int64_t g;
  for(g=1; g<argc; g++){
    if(strcmp(argv[g], sequence_align_instance_fname_parameter) == 0){
      if(g+1 < argc){
        strcpy(instance_fname, &argv[++g][0]);
      }
    }
  }
  check_preemptive_halt_sequence_align_ = check_preemptive_halt_sequence_align_deactivated;
  lcs_instance * instance = read_lcs_instance(instance_fname);
  nX=instance->num_first_sequence;
  X=calloc(nX,sizeof(int64_t));
  for(g=0; g<nX; g++){
    X[g]=instance->first_sequence[g];
  }
  nY=instance->num_second_sequence;
  Y=calloc(nY,sizeof(int64_t));
  for(g=0; g<nY; g++){
    Y[g]=instance->second_sequence[g];
  }
  pattern_references = 0;
  cache_reads = 0;
  cache_writes = 0;
  memo_calls = dist_calls = cache_misses = prev_cache_misses = 0;
}

void set_cache_miss_threshold_seq(int64_t t){
  prev_cache_misses = t;
}

void set_preemptive_halt_seq(int p){
  check_preemptive_halt_sequence_align_ = p ? check_preemptive_halt_sequence_align_activated : check_preemptive_halt_sequence_align_deactivated;
}

void reset_seq(long int p){
  memo_calls = dist_calls = cache_misses = 0;
}

long int get_cache_misses_seq(){
  return cache_misses;
}

void solve_seq(){
  //printf("solve sequence alignment\n");
  cost = alignment_cost(nX-1, nY-1);
}

void solve_sequence_alignment_standalone(int argc, char **argv) {
  initialize_long_int_cache(argc, argv);
  initialize_seq(argc, argv);
  cost = alignment_cost(nX-1, nY-1);
  //printf("cost %ld LRU queue size %ld cache_misses %ld\n", cost, lru_queue_size, cache_misses);
  align();
  view_scoring_grid();
}

/* PRIVATE */

void view_scoring_grid(){
  int64_t g,h;
  printf("Scoring Grid:\n");
  for(g = -1; g < nX; g++){
    for(h = -1; h < nY; h++){
      printf("%ld ", alignment_cost(g, h));
    }
    printf("\n");
  }
  printf("===================================\n");
}



char check_preemptive_halt_sequence_align_deactivated(){
  return 0;
}

char check_preemptive_halt_sequence_align_activated(){
  return (cache_misses > prev_cache_misses) ? 1 : 0;
}

uint64_t generate_key(int32_t i, int32_t j) {
  return (((int64_t) i) << 32) | (((int64_t) j) & 0xffffffffu);
}

int64_t memo_alignment_cost(long int len_s, long int len_t);

int64_t alignment_cost(long int len_s, long int len_t){
  int64_t * val, len;
  uint64_t * key;
  if(check_preemptive_halt_sequence_align_()){
    return 0;
  }
  dist_calls++;
  key = calloc(1, sizeof (int64_t));
  *key = generate_key(len_s, len_t);
  val = cache_read_long_int(key);
  cache_reads++;
  if(val == NULL){
    len = memo_alignment_cost(len_s, len_t);
    val = &len;
    cache_write_long_int(key, val);
    cache_writes++;
    cache_misses++;
  }
  len = *val;
  return len;
}

int X_seq(long int k){
  pattern_references++;
  //  if(pattern_reference_threshold < pattern_references){
  //    printf("pattern references exceeded threshold\n");
  //    exit(0);
  //  }
  return X[k];
}

int Y_seq(long int k){
  pattern_references++;
  //  if(pattern_reference_threshold < pattern_references){
  //    printf("pattern references exceeded threshold\n");
  //    exit(0);
  //  }
  return Y[k];
}

int64_t similarity(int64_t i, int64_t j){
  return X_seq(i) == Y_seq(j) ? 1 : -1;
}

int64_t memo_alignment_cost(long int i, long int j){
  memo_calls++;
  //  printf("%ld %ld\n", i,j);
  if(i <  0 && j >= 0) return -1 * j - 1;
  if(i >= 0 && j <  0) return -1 * i - 1;
  if(i <  0 && j <  0) return 0;
  int64_t c1 = alignment_cost(i - 1, j - 1) + similarity(i, j);
  int64_t c2 = alignment_cost(i - 1, j) - 1;
  int64_t c3 = alignment_cost(i, j - 1) - 1;
  //printf("three costs %ld %ld %ld\n", c1, c2, c3);
  int64_t max = (c1 > c2) ? c1 : c2;
  return (c3 > max) ? c3 : max;
}

void align(){
  //printf("align\n");
  int64_t longer_size = nX > nY ? nX : nY;
  int64_t * alignmentA = calloc(longer_size, sizeof(int64_t));
  int64_t * alignmentB = calloc(longer_size, sizeof(int64_t));
  int64_t * buffer = calloc(longer_size, sizeof(int64_t));
  int64_t g = 0;
  int64_t h = 0;
  int64_t i = nX-1;
  int64_t j = nY-1;
  while( i >= 0 || j >= 0){
    if( i >= 0 && j >= 0 && alignment_cost(i, j) == alignment_cost(i-1, j-1) + similarity(i, j) ){
      alignmentA[g++] = X[i--];
      alignmentB[h++] = Y[j--];
    } else 
    if( i >= 0 && alignment_cost(i, j) == alignment_cost(i-1, j) - 1 ){
      alignmentA[g++] = X[i--];
      alignmentB[h++] = LONG_MAX;
    } else {
      alignmentA[g++] = LONG_MAX;
      alignmentB[h++] = Y[j--];
    }
  }
  int64_t q = g - 1;
  int64_t p;
  for(p = 0; p < g; p++){
    buffer[p] = alignmentA[q--];
  }
  for(p = 0; p < g; p++){
    alignmentA[p] = buffer[p];
  }
  q = h - 1;
  for(p = 0; p < h; p++){
    buffer[p] = alignmentB[q--];
  }
  for(p = 0; p < h; p++){
    alignmentB[p] = buffer[p];
  }
  /*
  for(p = 0; p < g; p++){
    if(alignmentA[p] != LONG_MAX){
      printf("%ld ", alignmentA[p]);
    } else {
      printf("  ");
    }
  }
  printf("\n");
  for(p = 0; p < h; p++){
    if(alignmentB[p] != LONG_MAX){
      printf("%ld ", alignmentB[p]);
    } else {
      printf("  ");
    }
  }
  printf("\n");
  */
}


