/*******************************************************************************
 *  @file cloud_search.h
 *  @brief Testing for navigating through the matrices.
 *
 *  @author Dave Rich (devrek)
 *  @bug Lots.
 *******************************************************************************/


#ifndef _TESTING_H
#define _TESTING_H

void test_cycle(int Q, int T,
                float st_MX[ NUM_NORMAL_STATES * (Q + 1) * (T + 1) ],
                float sp_MX[ NUM_NORMAL_STATES * (Q + 1) ],
                TRACEBACK *tr);

void rev_test_cycle(int Q, int T,
                    float st_MX[ NUM_NORMAL_STATES * (Q + 1) * (T + 1) ],
                    float sp_MX[ NUM_NORMAL_STATES * (Q + 1) ],
                    TRACEBACK *tr);

void test_cloud(int Q, int T,
                float st_MX[ NUM_NORMAL_STATES * (Q + 1) * (T + 1) ],
                float sp_MX[ NUM_NORMAL_STATES * (Q + 1) ],
                EDGEBOUNDS* edg,
                float val);

#endif /* _TESTING_H */