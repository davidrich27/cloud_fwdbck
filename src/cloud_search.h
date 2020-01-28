/*******************************************************************************
 *  @file cloud_search.h
 *  @brief The "Cloud Search" Algorithm for the heuristic Forward-Backward.
 *
 *  @author Dave Rich (devrek)
 *  @bug Lots.
 *******************************************************************************/


#ifndef _CLOUD_SEARCH_H
#define _CLOUD_SEARCH_H

void cloud_forward_Run(const SEQ* query, 
                     const HMM_PROFILE* target,
                     int Q, int T, 
                     float* st_MX, 
                     float* sp_MX, 
                     TRACEBACK* tr,
                     EDGEBOUNDS* edg,
                     float alpha, int beta );

void cloud_backward_Run(const SEQ* query, 
                     const HMM_PROFILE* target,
                     int Q, int T, 
                     float* st_MX, 
                     float* sp_MX, 
                     TRACEBACK* tr,
                     EDGEBOUNDS* edg,
                     float alpha, int beta );

float forward_bounded_Run(const SEQ* query, 
                        const HMM_PROFILE* target,
                        int Q, int T, 
                        float* st_MX, 
                        float* sp_MX, 
                        EDGEBOUNDS* edg, 
                        float *sc_final );

float backward_bounded_Run(const SEQ* query, 
                        const HMM_PROFILE* target,
                        int Q, int T, 
                        float* st_MX, 
                        float* sp_MX,
                        EDGEBOUNDS* edg,
                        float *sc_final );

#endif /* _CLOUD_SEARCH_H */