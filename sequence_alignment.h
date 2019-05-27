/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   edit_distance.h
 * Author: benjamin
 *
 * Created on February 19, 2019, 4:19 PM
 */

#ifndef SEQUENCE_ALIGNMENT_H
#define SEQUENCE_ALIGNMENT_H

#ifdef __cplusplus
extern "C" {
#endif

  void initialize_sequence_alignment(int argc, char **argv);
  void set_cache_miss_threshold_sequence_alignment(int64_t t);
  void set_preemptive_halt_sequence_alignment(int p);
  void reset_sequence_alignment(long int p);
  long int get_cache_misses_sequence_alignment();
  void solve_sequence_alignment();
  void solve_sequence_alignment_standalone(int argc, char **argv);

#ifdef __cplusplus
}
#endif

#endif /* SEQUENCE_ALIGNMENT_H */

