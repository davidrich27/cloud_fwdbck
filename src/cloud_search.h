/*******************************************************************************
 *  @file cloud_search.h
 *  @brief The "Cloud Search" Algorithm for the heuristic Forward-Backward.
 *
 *  @author Dave Rich (devrek)
 *  @bug Lots.
 *******************************************************************************/


#ifndef _CLOUD_SEARCH_H
#define _CLOUD_SEARCH_H

void cloud_search_forward_Run(const SEQ* query, 
                              const HMM_PROFILE* target,
                              int Q, int T, 
                              float st_MX[ NUM_NORMAL_STATES * (Q+1) * (T+1) ], 
                              float sp_MX[ NUM_SPECIAL_STATES * (Q+1) ], 
                              RESULTS* res,
                              TRACEBACK* tr);

#endif /* _CLOUD_SEARCH_H */