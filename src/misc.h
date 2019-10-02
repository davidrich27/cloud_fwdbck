/*******************************************************************************
 *  @file misc.h
 *  @brief Miscellaneous helper functions.
 *
 *  @author Dave Rich (devrek)
 *  @bug Lots.
 *******************************************************************************/

#ifndef _MISC_H
#define _MISC_H

// constants
#define LOGSUM_SCALE 1000.f
#define LOGSUM_TBL 16000

/* Min/Max Fns for Viterbi */
float calc_Max (float x, float y);
float calc_Min (float x, float y);

/* Logsum fns for Forward/Backward */
void init_Logsum ();
float calc_Logsum (float x, float y);
float calc_Logsum_exact (float x, float y);
void print_Logsum ();

/* DP Matrix fns */
void dp_matrix_Print (const int Q, const int T, 
                      const float st_MX[ NUM_NORMAL_STATES * (Q+1) * (T+1) ], 
                      const float sp_MX[ NUM_SPECIAL_STATES * (Q+1) ]);
void dp_matrix_Clear (const int Q, const int T, 
                      float st_MX[ NUM_NORMAL_STATES * (Q+1) * (T+1) ], 
                      float sp_MX[ NUM_SPECIAL_STATES * (Q+1) ]);
void dp_matrix_Save (const int Q, const int T, 
                     const float st_MX[ NUM_NORMAL_STATES * (Q+1) * (T+1) ], 
                     const float sp_MX[ NUM_SPECIAL_STATES * (Q+1) ],
                     const char *_filename_);

#endif /* _MISC_H */