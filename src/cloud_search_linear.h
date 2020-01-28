/*******************************************************************************
 *  @file cloud_search3.h
 *  @brief The "Cloud Search" Algorithm for the heuristic Forward-Backward.
 *
 *  @author Dave Rich (devrek)
 *  @bug Lots.
 *******************************************************************************/


#ifndef _CLOUD_SEARCH_LINEAR_H
#define _CLOUD_SEARCH_LINEAR_H

void cloud_forward_Run3 (const SEQ* query,
                         const HMM_PROFILE* target,
                         int Q, int T,
                         float* st_MX3,
                         float* sp_MX,
                         TRACEBACK* tr,
                         EDGEBOUNDS* edg,
                         float alpha, int beta,
                         int *sc_final );

void cloud_backward_Run3 (const SEQ* query,
                          const HMM_PROFILE* target,
                          int Q, int T,
                          float* st_MX3,
                          float* sp_MX,
                          TRACEBACK* tr,
                          EDGEBOUNDS* edg,
                          float alpha, int beta,
                          int *sc_final );


#endif /* _CLOUD_SEARCH_LINEAR_H */