/*******************************************************************************
 *  @file forward_backward.c
 *  @brief The "Cloud Search" Algorithm for the heuristic Forward-Backward.
 *
 *  @synopsis
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
// NOTE: wrap all macro vars in parens!!
#define max(x,y) (((x) > (y)) ? (y) : (x))
#define min(x,y) (((x) > (y)) ? (x) : (y))

/*  
 *  FUNCTION: cloud_search_backward_Run()
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
   int    d,i,j,k,x = 0;        /* row, column indices */
   char   *seq = query->seq;  /* alias for getting seq */
   float  diag_max, diag_limit;  /* max prob score in the diag, and the pruning floor */
   int    d_cnt = 0;             /* number of anti-diags from starting position */
   int    left_bound, right_bound, left_new, right_new, left_edge, right_edge;
   float  cell_max;
   int    num_cells;


}


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

   char   a;                     /* store current character in sequence */
   int    A;                     /* store int value of character */
   int    d,i,j,k;               /* diagonal, row, column indices */
   char   *seq = query->seq;     /* alias for getting seq */
   float  diag_max, diag_limit;  /* max prob score in the diag, and the pruning floor */
   int    d_cnt = 0;             /* number of anti-diags from starting position */
   int    left_bound, right_bound, left_new, right_new, left_edge, right_edge;
   float  cell_max;
   int    num_cells;

   /* local or global (multiple alignments) */
   bool   is_local = false;
   float  sc_E = (is_local) ? 0 : -INF;

   /* X-drop ratio */
   float alpha = 0.01;
   /* number of anti-diag passes before pruning */
   float beta = 10;

   /* these are aligned arrays => diag_num: (lb, rb) */
   /* bounds of current anti-diag */
   int diag_num[1024];
   /* bounds of previous anti-diag */
   int diag_lb[1024];
   /* right bound of anti-diag */
   int diag_rb[1024];

   float cnt = 1.0;

   float  prev_mat, prev_del, prev_ins, prev_beg, prev_sum;
   // float  prev_mat, prev_del, prev_ins, prev_beg, prev_best;
   float  sc, sc_1, sc_2, sc_best, sc_max;

   /* for testing */
   dp_matrix_Clear (Q, T, st_MX, sp_MX);
   /* testing matrix */
   float test_MX[ (Q+1) * (T+1) ];
   for (i = 0; i <= Q; i++)
      for (j = 0; j <= T; j++)
         TMX(i,j) = 0;

   /* initialize special states (?) */
   XMX(SP_N,0) = 0;                                         /* S->N, p=1             */
   XMX(SP_B,0) = XSC(SP_N,SP_MOVE);                         /* S->N->B, no N-tail    */
   XMX(SP_E,0) = XMX(SP_C,0) = XMX(SP_J,0) = -INF;          /* need seq to get here (?)  */

   /* initialize 0 row (top-edge) and ) column (left-edge) */
   for (j = 0; j <= T; j++)
      MMX(0,j) = IMX(0,j) = DMX(0,j) = -INF;
   for (i = 0; i <= Q; i++)
      MMX(i,0) = IMX(i,0) = DMX(i,0) = -INF;

   /* get starting anti-diag */
   d = ABS( tr->start.i - tr->start.j );
   /* temp overwrite */
   d = 0;

   left_bound = right_bound = d;

   diag_max = -INF;
   num_cells = 0;

   /* increment by anti-diagonal */
   for (; d <= Q+T; d++)
   {
      if (d < Q && d < T)
         ++num_cells;
      if (d > Q && d > T)
         --num_cells;

      // printf("d= %d/%d, lb: %d, rb: %d\n", d, Q+T, left_bound, right_bound);

      /* TODO: Update Prune Bounds */
      /* if free passes are over (beta < d), prune and set new edgebounds */
      if (beta < d)
      {
         /* Traverse current bounds to find max score on diag */
         for (k = left_bound; k <= right_bound; k++)
         {
            i = k;
            j = d - i - 1; /* back one diag */
            diag_max = calc_Max( 
                           calc_Max( diag_max, MMX(i,j) ),
                           calc_Max( IMX(i,j), DMX(i,j) ) );
         }
         diag_limit = diag_max * alpha;
         printf("diag_max: %.2f diag_limit: %.2f\n", diag_max, diag_limit);

         /* Find the first cell from the left which passes above threshold */
         for (k = left_bound; k <= right_bound; k++)
         {
            i = k;
            j = d - i - 1;
            cell_max = calc_Max( MMX(i,j),
                           calc_Max( IMX(i,j), DMX(i,j) ) );
            if( cell_max >= diag_limit )
            {
               left_new = k - 2;
               break;
            }
         }

         /* Find the first cell from the right which passes above threshold */
         for (k = right_bound; k >= left_bound; k--)
         {
            i = k;
            j = d - i - 1;
            cell_max = calc_Max( MMX(i,j),
                           calc_Max( IMX(i,j),   DMX(i,j) ) );
            if( cell_max >= diag_limit )
            {
               right_new = k + 2;
               break;
            }
         }
      }
      else /* else edges expand in square pattern */
      {
         left_new = 0;
         right_new = d;
      }
      
      /* update pruned bounds */
      right_bound = right_new;
      left_bound = left_new;
      diag_max = -INF;

      /* TODO: edge-check pruned bounds
      /* bounds-check of the matrix and check that left and right are inside */
      left_edge = calc_Max(0, d - T);
      right_edge = left_edge + num_cells;
      left_bound = calc_Max(left_edge, left_bound);
      right_bound = calc_Min(right_edge, right_bound);

      printf("d= %d/%d, lb: %d, rb: %d, ", d, Q+T, left_bound, right_bound);
      printf("le: %d, re: %d\n", left_edge, right_edge);

      /* left to right across anti-diagonal */
      for (k = left_bound; k <= right_bound; k++)
      {
         i = k;
         j = d - i; 

         printf("d=%d, i=%d, j=%d\n", d, i, j);

         TMX(i,j) = 1;
         cnt++;

         /* TODO: Forward (pull values from previous state cells) */

         /* Push from Match State to Next */
         a = seq[i];
         A = AA_REV[a];
         
         /* FIND BEST PATH TO MATCH STATE (FROM MATCH, INSERT, DELETE, OR BEGIN) */
         /* best previous state transition (match takes the diag element of each prev state) */
         prev_mat = MMX(i-1,j-1)  + TSC(j-1,M2M);
         prev_ins = IMX(i-1,j-1)  + TSC(j-1,I2M);
         prev_del = DMX(i-1,j-1)  + TSC(j-1,D2M);
         prev_beg = XMX(SP_B,i-1) + TSC(j-1,B2M); /* from begin match state (new alignment) */
         /* best-to-match */
         prev_sum = calc_Logsum( 
                        calc_Logsum( prev_mat, prev_ins ),
                        calc_Logsum( prev_del, prev_beg )
                     );
         MMX(i,j) = prev_sum + MSC(j,A);

         /* FIND BEST PATH TO INSERT STATE (FROM MATCH OR INSERT) */
         /* previous states (match takes the left element of each state) */
         prev_mat = MMX(i-1,j) + TSC(j,M2I);
         prev_ins = IMX(i-1,j) + TSC(j,I2I);
         /* best-to-insert */
         prev_sum = calc_Logsum( prev_mat, prev_ins );
         IMX(i,j) = prev_sum + ISC(j,A);

         /* FIND BEST PATH TO DELETE STATE (FOMR MATCH OR DELETE) */
         /* previous states (match takes the left element of each state) */
         prev_mat = MMX(i,j-1) + TSC(j-1,M2D);
         prev_del = DMX(i,j-1) + TSC(j-1,D2D);
         /* best-to-delete */
         prev_sum = calc_Logsum(prev_mat, prev_del);
         DMX(i,j) = prev_sum;

         /* UPDATE E STATE */
         XMX(SP_E,i) = calc_Logsum( XMX(SP_E,i), 
                                    MMX(i,j) + sc_E );
         XMX(SP_E,i) = calc_Logsum( XMX(SP_E,i), 
                                    DMX(i,j) + sc_E );
      }

      ++d_cnt;
   }

   dp_matrix_Print(Q, T, st_MX, sp_MX);
   test_matrix_Print(Q, T, test_MX);
   return;

   }


/*  
 *  FUNCTION: forward_bounded_Run()
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
 *             <bnd>       Bounds Data
 *
 *  RETURN: 
 */
void forward_bounded_Run(const SEQ* query, 
                        const HMM_PROFILE* target,
                        int Q, int T, 
                        float st_MX[ NUM_NORMAL_STATES * (Q+1) * (T+1) ], 
                        float sp_MX[ NUM_SPECIAL_STATES * (Q+1) ], 
                        RESULTS* res)
{
}
