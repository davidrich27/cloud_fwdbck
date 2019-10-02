/*******************************************************************************
 *  @file forward_backward.c
 *  @brief The "Cloud Search" Algorithm for the heuristic Forward-Backward.
 *
 *  @author Dave Rich (devrek)
 *  @bug Lots.
 *******************************************************************************/

// imports
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>

// local imports (after struct declarations)
#include "structs.h"
#include "misc.h"
#include "hmm_parser.h"
#include "viterbi.h"
#include "forward_backward.h"
#include "cloud_search.h"

// macros
#define getName(var) #var
#define SCALE_FACTOR 1000

// macro functions
#define max(x,y) ((x > y) ? y : x)
#define min(x,y) ((x > y) ? x : y)

/*  
 *  FUNCTION: cloud_search_forward_Run()
 *  SYNOPSIS: Perform Forward part of Cloud Search Algorithm.
 *
 *  PURPOSE:
 *
 *  ARGS:      <query>     query sequence, 
 *             <target>    HMM model,
 *             <Q>         query length, 
 *             <T>         target length,
 *             <st_MX>     Normal State (Match, Insert, Delete) Matrix,
 *             <sp_MX>     Special State (J,N,B,C,E) Matrix,
 *             <res>       Results Data
 *             <tr>        Traceback Data
 *
 *  RETURN: 
 */
void cloud_search_forward_Run(const SEQ* query, 
                              const HMM_PROFILE* target,
                              int Q, int T, 
                              float st_MX[ NUM_NORMAL_STATES * (Q+1) * (T+1) ], 
                              float sp_MX[ NUM_SPECIAL_STATES * (Q+1) ], 
                              RESULTS* res,
                              TRACEBACK* tr)
{
   printf("cloud forward search...\n");


   char   a;                  /* store current character in sequence */
   int    A;                  /* store int value of character */
   int    d,i,j,k = 0;        /* row, column indices */
   char   *seq = query->seq;  /* alias for getting seq */
   int diag_max;
   int d_cnt = 0;             /* number of anti-diags from starting position */
   int left_bound, right_bound;

   /* X-drop ratio */
   float alpha = 1e-4;
   /* number of anti-diag passes before pruning */
   float beta = 3;

   /* format: [number of bounds on anti-diag], [i-idx start (1)], [i-idx end (1)]... */
   /* bounds of current anti-diag */
   int curr_bounds[1024];
   /* bounds of previous anti-diag */
   int prev_bounds[1024];
   /* bounds of previous-previous anti-diag */
   int prev2_bounds[1024];

   float cnt = 1.0;

   float  prev_mat, prev_del, prev_ins, prev_beg, prev_best;
   float  sc, sc_1, sc_2, sc_best, sc_max;

   dp_matrix_Clear (Q, T, st_MX, sp_MX);

   /* get starting anti-diag */
   d = ABS( tr->start.i - tr->start.j );
   /* temp overwrite */
   d = 0;
   left_bound = d;
   right_bound = d;

   /* increment by anti-diagonal */
   for (; d <= Q+T; d++)
   {
      left_bound = d - d_cnt;
      right_bound = d + d_cnt;
      diag_max = -INF;

      /* left to right across anti-diagonal */
      for (i = d - left_bound; i >= 0 + right_bound; i--)
      {
         j = d - i; 
         // printf("d=%d, i=%d, j=%d\n", d, i, j);

         /* edge-bound check */
         if (i < 0 || i > Q || j < 0 || j > T)
         {
            continue;
         }

         MMX(i,j) = cnt;
         cnt++;

         /* TODO: Forward (push values to next state cells) */

         /* Push from Match State to Next */
         a = seq[i+1];
         A = AA_REV[a];
         MMX(i+1,j+1) = calc_Logsum( 
                           MMX(i+1,j+1), 
                           MMX(i,j) + TSC(j,M2M) + MSC(j+1,A) );
         IMX(i+1,j+1) = calc_Logsum(
                           IMX(i+1,j+1),
                           MMX(i,j) + TSC(j,M2I) + ISC(j+1,A) );
         DMX(i+1,j+1) = calc_Logsum(
                           DMX(i+1,j+1),
                           MMX(i,j) + TSC(j,M2D) );

         /* Push from Insert State to Next */
         MMX(i+1,j) = calc_Logsum( 
                           MMX(i+1,j), 
                           IMX(i,j) + TSC(j,I2M) + MSC(j,A) );
         IMX(i+1,j) = calc_Logsum(
                           IMX(i+1,j),
                           IMX(i,j) + TSC(j,I2I) + ISC(j,A) );

         /* Push from Delete State to Next */
         a = seq[i];
         A = AA_REV[a];
         MMX(i,j+1) = calc_Logsum( 
                           MMX(i,j+1), 
                           DMX(i,j) + TSC(j,D2M) + MSC(j+1,A) );
         DMX(i,j+1) = calc_Logsum(
                           DMX(i,j+1),
                           DMX(i,j) + TSC(j,D2D) );

         /* Update E State */

      }

      /* TODO: Determine Pruning Bounds */
      ++d_cnt;
      if (d >= beta)
      {

      }
      else 
      {
         --left_bound;
         ++right_bound;
      }
   }

   dp_matrix_Print(Q, T, st_MX, sp_MX);
   return;

   // TODO: efficiently check bounds
   // MIN = min( Q, T );
   // MAX = max( Q, T );
   // /* increment by anti-diagonal */
   // for (d = 0; d <= MIN; d++)
   // {
   //    /* left to right across anti-diagonal */
   //    for (i = d; i >= 0; i--)
   //    {
   //       j = d - i; 
   //       printf("d=%d, i=%d, j=%d\n", d, i, j);

   //       MAT[i][j] = cnt;
   //       cnt++;

   //       /* TODO: Forward (push values to next state cells) */


   //       /* TODO: Determine Pruning Bounds */
   //    }
   // }
   
   // /* begin dp matrix taper on one side */
   // for (d = MIN+1; d < MAX; d++)
   // {
   //    /* if dp matrix is tall */
   //    if (Q > T) 
   //    {
   //       /* left to right across anti-diagonal */
   //       for (i = d; i > 0; i--)
   //       {

   //       }
   //    }
   //    /* if dp matrix is wide */
   //    if (T > Q)
   //    {
   //       /* left to right across anti-diagonal */
   //       for (i = d; i > 0; i--)
   //       {

   //       }
   //    }

   // }

   // /* begin dp matrix taper on both sides */
   // for (d = MAX+1; d < Q+T; d++)
   // {

   // }

   }


/*  
 *  FUNCTION: cloud_search_backward_Run()
 *  SYNOPSIS: Perform Backward part of Cloud Search Algorithm.
 *
 *  PURPOSE:
 *
 *  ARGS:      <query>     query sequence, 
 *             <target>    HMM model,
 *             <Q>         query length, 
 *             <T>         target length,
 *             <st_MX>     Normal State (Match, Insert, Delete) Matrix,
 *             <sp_MX>     Special State (J,N,B,C,E) Matrix,
 *             <res>       Results Data
 *             <tr>        Traceback Data
 *
 *  RETURN: 
 */
void cloud_search_backward_Run(const SEQ* query, 
                              const HMM_PROFILE* target,
                              int Q, int T, 
                              float st_MX[ NUM_NORMAL_STATES * (Q+1) * (T+1) ], 
                              float sp_MX[ NUM_SPECIAL_STATES * (Q+1) ], 
                              RESULTS* res,
                              TRACEBACK* tr)
{
   printf("cloud backward search...\n");

   char   a;                  /* store current character in sequence */
   int    A;                  /* store int value of character */
   int    d,i,j,k = 0;        /* row, column indices */
   char   *seq = query->seq;  /* alias for getting seq */
   int    diag_max;
   int    d_cnt = 0;             /* number of anti-diags from starting position */
   int    left_bound, right_bound, left_prev, right_prev;

   /* X-drop ratio */
   float alpha = 1e-4;
   /* number of anti-diag passes before pruning */
   float beta = 3;

   /* format: [number of bounds on anti-diag], [i-idx start (1)], [i-idx end (1)]... */
   /* bounds of current anti-diag */
   int curr_bounds[1024];
   /* bounds of previous anti-diag */
   int prev_bounds[1024];
   /* bounds of previous-previous anti-diag */
   int prev2_bounds[1024];

   float cnt = 1.0;

   float  prev_mat, prev_del, prev_ins, prev_beg, prev_best;
   float  sc, sc_1, sc_2, sc_best, sc_max;

   dp_matrix_Clear (Q, T, st_MX, sp_MX);

   /* get starting anti-diag */
   d = ABS( tr->start.i - tr->start.j );
   /* temp overwrite */
   d = 0;
   left_bound = d;
   right_bound = d;

   left_prev;
   right_prev;

   /* increment by anti-diagonal */
   for (; d <= Q+T; d++)
   {
      /* prune and set new edgebounds */
      if (d >= beta)
      {
         for (i = d - left_bound; i >= 0 + right_bound; i--)
         {
            
         }
      }
      else /* else edges expand in square pattern */
      {
         left_bound = 0;
         right_bound = d;
      }
      
      diag_max = -INF;

      /* left to right across anti-diagonal */
      for (i = d - left_bound; i >= 0 + right_bound; i--)
      {
         j = d - i; 
         // printf("d=%d, i=%d, j=%d\n", d, i, j);

         /* edge-bound check */
         if (i < 0 || i > Q || j < 0 || j > T)
         {
            continue;
         }

         MMX(i,j) = cnt;
         cnt++;

         /* TODO: Forward (push values to next state cells) */

         /* Push from Match State to Next */
         a = seq[i+1];
         A = AA_REV[a];
         MMX(i+1,j+1) = calc_Logsum( 
                           MMX(i+1,j+1), 
                           MMX(i,j) + TSC(j,M2M) + MSC(j+1,A) );
         IMX(i+1,j+1) = calc_Logsum(
                           IMX(i+1,j+1),
                           MMX(i,j) + TSC(j,M2I) + ISC(j+1,A) );
         DMX(i+1,j+1) = calc_Logsum(
                           DMX(i+1,j+1),
                           MMX(i,j) + TSC(j,M2D) );

         /* Push from Insert State to Next */
         MMX(i+1,j) = calc_Logsum( 
                           MMX(i+1,j), 
                           IMX(i,j) + TSC(j,I2M) + MSC(j,A) );
         IMX(i+1,j) = calc_Logsum(
                           IMX(i+1,j),
                           IMX(i,j) + TSC(j,I2I) + ISC(j,A) );

         /* Push from Delete State to Next */
         a = seq[i];
         A = AA_REV[a];
         MMX(i,j+1) = calc_Logsum( 
                           MMX(i,j+1), 
                           DMX(i,j) + TSC(j,D2M) + MSC(j+1,A) );
         DMX(i,j+1) = calc_Logsum(
                           DMX(i,j+1),
                           DMX(i,j) + TSC(j,D2D) );

         /* Update E State */

      }

      /* TODO: Determine Pruning Bounds */
      ++d_cnt;
   }

   dp_matrix_Print(Q, T, st_MX, sp_MX);
   return;

   // TODO: efficiently check bounds
   // MIN = min( Q, T );
   // MAX = max( Q, T );
   // /* increment by anti-diagonal */
   // for (d = 0; d <= MIN; d++)
   // {
   //    /* left to right across anti-diagonal */
   //    for (i = d; i >= 0; i--)
   //    {
   //       j = d - i; 
   //       printf("d=%d, i=%d, j=%d\n", d, i, j);

   //       MAT[i][j] = cnt;
   //       cnt++;

   //       /* TODO: Forward (push values to next state cells) */


   //       /* TODO: Determine Pruning Bounds */
   //    }
   // }
   
   // /* begin dp matrix taper on one side */
   // for (d = MIN+1; d < MAX; d++)
   // {
   //    /* if dp matrix is tall */
   //    if (Q > T) 
   //    {
   //       /* left to right across anti-diagonal */
   //       for (i = d; i > 0; i--)
   //       {

   //       }
   //    }
   //    /* if dp matrix is wide */
   //    if (T > Q)
   //    {
   //       /* left to right across anti-diagonal */
   //       for (i = d; i > 0; i--)
   //       {

   //       }
   //    }

   // }

   // /* begin dp matrix taper on both sides */
   // for (d = MAX+1; d < Q+T; d++)
   // {

   // }

   }
