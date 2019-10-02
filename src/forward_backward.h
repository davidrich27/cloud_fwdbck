/*******************************************************************************
 *  @file forward_backward.h
 *  @brief Function prototypes for the Forward-Backward algorithm.
 *
 *  @author Dave Rich (devrek)
 *  @bug Lots.
 *******************************************************************************/

#ifndef _FORWARDBACKWARD_H
#define _FORWARDBACKWARD_H

void forward_backward_Run (const SEQ* query, 
                           const HMM_PROFILE* target,
                           int Q, int T, 
                           float st_MX[ NUM_NORMAL_STATES * (Q+1) * (T+1) ], 
                           float sp_MX[ NUM_SPECIAL_STATES * (Q+1) ], 
                           RESULTS* res);
float forward_Run (const SEQ* query, 
                   const HMM_PROFILE* target, 
                   int Q, int T, 
                   float st_MX[ NUM_NORMAL_STATES * (Q+1) * (T+1) ], 
                   float sp_MX[ NUM_SPECIAL_STATES * (Q+1) ], 
                   RESULTS* res);
float backward_Run (const SEQ* query, 
                    const HMM_PROFILE* target, 
                    int Q, int T, 
                    float st_MX[ NUM_NORMAL_STATES * (Q+1) * (T+1) ], 
                    float sp_MX[ NUM_SPECIAL_STATES * (Q+1) ], 
                    RESULTS* res);

#endif /* _FORWARDBACKWARD_H */
