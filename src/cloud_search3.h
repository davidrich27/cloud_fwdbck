/*******************************************************************************
 *  @file cloud_search3.h
 *  @brief The "Cloud Search" Algorithm for the heuristic Forward-Backward.
 *
 *  @author Dave Rich (devrek)
 *  @bug Lots.
 *******************************************************************************/


#ifndef _CLOUD_SEARCH_H
#define _CLOUD_SEARCH_H

/* st_MX [ (NUM_NORMAL_STATES) * (Q + 1) * (3) ] */
/* sp_MX [ (NUM_SPECIAL_STATES) * (Q + 1) * (T + 1) ] */

void cloud_forward_Run3 (const SEQ* query,
                         const HMM_PROFILE* target,
                         int Q, int T,
                         float* st_MX3,
                         float* sp_MX,
                         TRACEBACK* tr,
                         EDGEBOUNDS* edg,
                         float alpha, int beta,
                         int *sc_final )

void cloud_backward_Run3 (const SEQ* query,
                          const HMM_PROFILE* target,
                          int Q, int T,
                          float* st_MX3,
                          float* sp_MX,
                          TRACEBACK* tr,
                          EDGEBOUNDS* edg,
                          float alpha, int beta,
                          int *sc_final )

float forward_bounded_Run3(const SEQ* query,
                           const HMM_PROFILE* target,
                           int Q, int T,
                           float* st_MX3,
                           float* sp_MX,
                           EDGEBOUNDS* edg, 
                           int *sc_final )

float backward_bounded_Run3(const SEQ* query,
                            const HMM_PROFILE* target,
                            int Q, int T,
                            float* st_MX3,
                            float* sp_MX,
                            EDGEBOUNDS* edg, 
                            int *sc_final )


#endif /* _CLOUD_SEARCH_H */