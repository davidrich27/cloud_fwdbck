/*******************************************************************************
 *  @file forward_backward.c
 *  @brief Testing for navigating through the matrices.
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
#define max(x,y) (((x) > (y)) ? (x) : (y))
#define min(x,y) (((x) < (y)) ? (x) : (y))


/*  
 *  FUNCTION: cloud_backward_Run()
 *  SYNOPSIS: Perform Backward part of Cloud Search Algorithm (Quadratic Space Implementation).
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
 *             <edg>       Edge Bounds Tracker Data
 *            tuning variables:
 *             <alpha>     Pruning Drop 
 *             <beta>      Number of Passes before
 *            meta-data:
 *             <comp_num>  Number of Cells Computed
 *
 *  RETURN: 
 */
void cloud_backward_Run(const SEQ* query, 
                        const HMM_PROFILE* target,
                        int Q, int T, 
                        float st_MX[ NUM_NORMAL_STATES * (Q+1) * (T+1) ], 
                        float sp_MX[ NUM_SPECIAL_STATES * (Q+1) ], 
                        RESULTS* res,
                        TRACEBACK* tr,
                        EDGEBOUNDS* edg,
                        float alpha, int beta )
{
   /* vars for navigating matrix */
   int d,i,j,k,b;                /* diagonal, row, column, ... indices */
   int lb, rb, le, re;           /* right/left bounds and edges */
   int lb_new, rb_new;           /* tmp vars for new right/left bounds */
   int num_cells;                /* number of cells in diagonal */
   int d_st, d_end, d_cnt;       /* starting and ending diagonal indices */
   int dim_min, dim_max;         /* diagonal index where num cells reaches highest point and diminishing point */ 
   int dim_T, dim_Q;             /* dimensions of submatrix being searched */
   int total_cnt;                /* number of cells computed */
   COORDS start, end;            /* start and end point of alignment */

   /* vars for computing cells */
   char   a;                     /* store current character in sequence */
   int    A;                     /* store int value of character */
   char   *seq = query->seq;     /* alias for getting seq */

   /* vars for recurrance */
   float  prev_mat, prev_del, prev_ins, prev_beg, prev_end, prev_sum;
   float  sc, sc_1, sc_2, sc_best, sc_max;
   float  sc_M, sc_I, sc_D;

   /* vars for pruning */
   float  cell_max, total_max, diag_max;  /* max scores for current cell, diagonal, and entire matrix */
   float  diag_limit;                     /* pruning score */

   /* local or global? (multiple alignments) */
   bool   is_local = false;
   float  sc_E = (is_local) ? 0 : -INF;

   /* INIT EDGEBOUND DATA STRUCT */
   int min_size = 256;
   edg->size = min_size;
   edg->bounds = (BOUND *)malloc( min_size * sizeof(BOUND) );

   /* track number of cells computed */
   int cpu_num = 0;

   dp_matrix_Clear(Q, T, st_MX, sp_MX);

   /* INIT ANTI-DIAGS */
   start = tr->first_m;
   end = tr->last_m;

   /* TESTING */
   // alpha = 6.0;
   // beta = 5;
   // start = (COORDS) {1, 1};
   // end = (COORDS) {T, Q}; 
   printf("T=%d, Q=%d, alpha=%.1f, beta=%d\n", T, Q, alpha, beta);

   /* diag index at corners of dp matrix */
   d_st = 0;
   d_end = Q + T;
   d_cnt = 0;

   /* diag index of different start points, creating submatrix */
   // d_st = start.i + start.j;
   d_end = end.i + end.j;

   /* dimension of submatrix */
   dim_Q = end.i;
   dim_T = end.j;

   /* diag index where num cells reaches highest point and begins diminishing */
   dim_min = min(end.i, end.j);
   dim_max = max(end.i, end.j);

   /* set bounds using starting cell */
   lb = end.i;
   rb = end.i + 1;
   num_cells = 0;

   /* resize if needed */
   if (edg->N >= edg->size) {
      edg->size *= 2;
      edg->bounds = realloc(edg->bounds, edg->size * sizeof(BOUND) );
   }

   /* INITIALIZE SPECIAL STATES */
   /* Initialize the Q row. */
   XMX(SP_J,Q) = XMX(SP_B,Q) = XMX(SP_N,Q) = -INF;
   XMX(SP_C,Q) = XSC(SP_C,SP_MOVE);
   XMX(SP_E,Q) = XMX(SP_C,Q) + XSC(SP_E,SP_MOVE);

   MMX(Q,T) = DMX(Q,T) = XMX(SP_E,Q);
   IMX(Q,T) = -INF;

   for (j = T-1; j >= 1; j--)
   {
      MMX(Q,j) = calc_Logsum( XMX(SP_E,Q) + sc_E, 
                              DMX(Q,j+1)  + TSC(j,M2D) );
      DMX(Q,j) = calc_Logsum( XMX(SP_E,Q) + sc_E,
                              DMX(Q,j+1)  + TSC(j,D2D) );
      IMX(Q,j) = -INF;
   }

   /* compute special states probs up to first i-index */
   for (i = Q-1; i >= end.i; i--) {

      /* Push from Match State to Next */
      a = seq[i];
      A = AA_REV[a];

      /* SPECIAL STATES */
      XMX(SP_B,i) = MMX(i+1,1) + TSC(0,B2M) + MSC(1,A);

      /* B -> MATCH */
      for (j = 2; j <= T; j++)
      {
         XMX(SP_B,i) = calc_Logsum( XMX(SP_B,i),
                                    MMX(i+1,j) + TSC(j-1,B2M) + MSC(j,A) );
      }

      XMX(SP_J,i) = calc_Logsum( XMX(SP_J,i+1) + XSC(SP_J,SP_LOOP),
                                 XMX(SP_B,i)   + XSC(SP_J,SP_MOVE) );

      XMX(SP_C,i) = XMX(SP_C,i+1) + XSC(SP_C,SP_LOOP);

      XMX(SP_E,i) = calc_Logsum( XMX(SP_J,i) + XSC(SP_E,SP_LOOP),
                                 XMX(SP_C,i) + XSC(SP_E,SP_MOVE) );

      XMX(SP_N,i) = calc_Logsum( XMX(SP_N,i+1) + XSC(SP_N,SP_LOOP),
                                 XMX(SP_B,i)   + XSC(SP_N,SP_MOVE) );
   }

   /* ITERATE THROUGHT ANTI-DIAGONALS */
   for (d = d_end; d >= d_st; d--, d_cnt++)
   {
      /* Is dp matrix diagonal growing or shrinking? */
      if (d >= dim_max)
         num_cells++;
      if (d < dim_min)
         num_cells--;

      /* UPDATE BOUNDS */
      /* if free passes are complete (beta < d), prune and set new edgebounds */
      if (beta < d_cnt)
      {
         lb_new = -100;
         rb_new = -100;

         /* Traverse current bounds to find max score on diag */
         diag_max = -INF;
         for (k = lb; k < rb; k++)
         {
            i = k;
            j = (d + 1) - i;    /* back one diag */
            diag_max = calc_Max( 
                           calc_Max( diag_max, MMX(i,j) ),
                           calc_Max( IMX(i,j), DMX(i,j) ) );
            // printf("diag_max(%d,%d) = %.2f \n", i,j,diag_max);
         }

         /* total max records largest cell score see so far */
         if (diag_max > total_max)
            total_max = diag_max;

         /* set score threshold for pruning */
         // diag_limit = diag_max - alpha;
         diag_limit = total_max - alpha;

         // printf("total_max: %.2f diag_max: %.2f diag_limit: %.2f, alph=%.1f\n", total_max, diag_max, diag_limit, alpha);

         /* Find the first cell from the left which passes above threshold */
         for (k = lb; k < rb; k++)
         {
            i = k;
            j = (d + 1) - i;
            cell_max = calc_Max( MMX(i,j),
                           calc_Max( IMX(i,j), DMX(i,j) ) );
            // printf("> TEST(%d,%d)=%.2f\n", i, j, cell_max);

            if( cell_max >= diag_limit )
            {
               lb_new = i;
               break;
            }
         }

         /* If no boundary edges are found on diag, then branch is pruned entirely and we are done */
         if (lb_new == -100)
            break;

         /* Find the first cell from the right which passes above threshold */
         for (k = rb - 1; k >= lb; k--)
         {
            i = k;
            j = (d + 1) - i;
            cell_max = calc_Max( MMX(i,j),
                           calc_Max( IMX(i,j),   DMX(i,j) ) );
            // printf("< TEST(%d,%d)=%.2f\n", i, j, cell_max);

            if( cell_max >= diag_limit )
            {
               rb_new = (i + 1);
               break;
            }
         }
      }
      else /* else edges expand in square pattern */
      {
         lb_new = lb;
         rb_new = rb;
      }

      /* Update Bounds */
      lb = lb_new - 1;
      rb = rb_new;

      /* Edge-check: find diag cells that are inside matrix bounds */
      le = max(end.i - (d_end - d) + 1, 0);
      re = le + num_cells;

      if (lb < le)
         lb = le;
      if (rb > re)
         rb = re;

      /* ADD NEW ANTI-DIAG TO EDGEBOUNDS */
      edg->bounds[edg->N].lb = lb;
      edg->bounds[edg->N].rb = rb;
      edg->bounds[edg->N].diag = d;
      edg->N += 1;
      /* resize if needed */
      if (edg->N >= edg->size) {
         edg->size *= 2;
         edg->bounds = (BOUND *)realloc(edg->bounds, edg->size * sizeof(BOUND) );
      }

      // printf("%d: new[%d,%d] -> edg[%d,%d] -> [%d,%d]\n", d, lb_new, rb_new, le, re, lb, rb);

      /* ITERATE THROUGH CELLS OF ANTI-DIAGONAL */
      for (k = lb; k < rb; k++, total_cnt++)
      {
         i = k;
         j = d - i;

         /* Get next sequence character */
         a = seq[i];
         A = AA_REV[a];

         /* SPECIAL STATES */
         XMX(SP_B,i) = MMX(i+1,1) + TSC(0,B2M) + MSC(1,A);

         /* B -> MATCH */
         for (b = 2; b <= T; b++)
         {
            XMX(SP_B,i) = calc_Logsum( XMX(SP_B,i),
                                       MMX(i+1,j) + TSC(j-1,B2M) + MSC(j,A) );
         }

         XMX(SP_J,i) = calc_Logsum( XMX(SP_J,i+1) + XSC(SP_J,SP_LOOP),
                                    XMX(SP_B,i)   + XSC(SP_J,SP_MOVE) );

         XMX(SP_C,i) = XMX(SP_C,i+1) + XSC(SP_C,SP_LOOP);

         XMX(SP_E,i) = calc_Logsum( XMX(SP_J,i) + XSC(SP_E,SP_LOOP),
                                    XMX(SP_C,i) + XSC(SP_E,SP_MOVE) );

         XMX(SP_N,i) = calc_Logsum( XMX(SP_N,i+1) + XSC(SP_N,SP_LOOP),
                                 XMX(SP_B,i)   + XSC(SP_N,SP_MOVE) );
         MMX(i,T) = DMX(i,T) = XMX(SP_E,i);
         IMX(i,T) = -INF;
         
         /* MAIN RECURSION */
         sc_M = MSC(j+1,A);
         sc_I = ISC(j+1,A);

         /* FIND SUM OF PATHS FROM MATCH, INSERT, DELETE, OR END STATE (TO PREVIOUS MATCH) */
         prev_mat = MMX(i+1,j+1) + TSC(j,M2M) + sc_M;
         prev_ins = IMX(i+1,j)   + TSC(j,M2I) + sc_I;
         prev_del = DMX(i,j+1)   + TSC(j,M2D);
         prev_end = XMX(SP_E,i)  + sc_E;     /* from end match state (new alignment) */
         /* best-to-match */
         prev_sum = calc_Logsum( 
                           calc_Logsum( prev_mat, prev_ins ),
                           calc_Logsum( prev_del, prev_end ) );
         MMX(i,j) = prev_sum;

         /* FIND SUM OF PATHS FROM MATCH OR INSERT STATE (TO PREVIOUS INSERT) */
         prev_mat = MMX(i+1,j+1) + TSC(j,I2M) + sc_M;
         prev_ins = IMX(i+1,j)   + TSC(j,I2I) + sc_I;
         /* best-to-insert */
         prev_sum = calc_Logsum( prev_mat, prev_ins );
         IMX(i,j) = prev_sum;

         /* FIND SUM OF PATHS FROM MATCH OR DELETE STATE (FROM PREVIOUS DELETE) */
         prev_mat = MMX(i+1,j+1) + TSC(j,D2M) + sc_M;
         prev_del = DMX(i,j+1)   + TSC(j,D2D);
         prev_end = XMX(SP_E,i)  + sc_E;
         /* best-to-delete */
         prev_sum = calc_Logsum( prev_mat, prev_del );
         prev_sum = calc_Logsum( prev_sum, prev_end );
         DMX(i,j) = prev_sum;

         // printf("CALC %d=>(%d,%d): %.1f %.1f %.1f \n", d, i, j, MMX(i,j),IMX(i,j),DMX(i,j));
      }
   }

   // printf("BACKWARD:\n");
   // printf("MAX: %d\n", edg->N);
   // edgebounds_Print(edg);

   /* reverse order of diagonals */
   BOUND tmp;
   for (i = 0; i < (edg->N / 2); ++i)
   {
      tmp = edg->bounds[i];
      edg->bounds[i] = edg->bounds[edg->N-i];
      edg->bounds[edg->N-i] = tmp;
   }
}  


/*  
 *  FUNCTION: cloud_forward_Run()
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
 *             <edg>       Edge Bounds Tracker Data
 *             <alpha>     Pruning Drop 
 *             <beta>      Number of Passes before
 *
 *  RETURN: 
 */
void cloud_forward_Run(const SEQ* query, 
                     const HMM_PROFILE* target,
                     int Q, int T, 
                     float st_MX[ NUM_NORMAL_STATES * (Q+1) * (T+1) ], 
                     float sp_MX[ NUM_SPECIAL_STATES * (Q+1) ], 
                     RESULTS* res,
                     TRACEBACK* tr,
                     EDGEBOUNDS* edg,
                     float alpha, int beta )
{
   /* vars for navigating matrix */
   int d,i,j,k;                  /* diagonal, row, column indices */
   int lb, rb, le, re;           /* right/left bounds and edges */
   int lb_new, rb_new;           /* tmp vars for new right/left bounds */
   int num_cells;                /* number of cells in diagonal */
   int d_st, d_end, d_cnt;       /* starting and ending diagonal indices */
   int dim_min, dim_max;         /* diagonal index where num cells reaches highest point and diminishing point */ 
   int dim_T, dim_Q;             /* dimensions of submatrix being searched */
   int total_cnt;                /* number of cells computed */
   COORDS start, end;            /* start and end point of alignment */

   /* vars for computing cells */
   char   a;                     /* store current character in sequence */
   int    A;                     /* store int value of character */
   char   *seq = query->seq;     /* alias for getting seq */

   /* vars for recurrance */
   float  prev_mat, prev_del, prev_ins, prev_beg, prev_end, prev_sum;
   float  sc, sc_1, sc_2, sc_best, sc_max;
   float  sc_M, sc_I, sc_D;
   /* vars for pruning */
   float  cell_max, total_max, diag_max, diag_limit;

   /* local or global? (multiple alignments) */
   bool   is_local = false;
   float  sc_E = (is_local) ? 0 : -INF;

   /* INIT EDGEBOUND DATA STRUCT */
   int min_size = 128;
   edg->size = min_size;
   edg->bounds = (BOUND *)malloc( min_size * sizeof(BOUND) );

   /* track number of cells computed */
   int cpu_num = 0;

   dp_matrix_Clear(Q, T, st_MX, sp_MX);

   /* INIT ANTI-DIAGS */
   /* test coords */
   start = tr->first_m;
   end = tr->last_m;

   /* TESTING */
   // alpha = 6.0;
   // beta = 5;
   // start = (COORDS) {1, 1};
   // end = (COORDS) {T, Q};
   printf("T=%d, Q=%d, alpha=%.1f, beta=%d\n", T, Q, alpha, beta);
   printf("start: (%d,%d), end: (%d,%d)\n", start.i, start.j, end.i, end.j);

   /* diag index at corners of dp matrix */
   d_st = 0;
   d_end = Q + T;
   d_cnt = 0;

   /* diag index of different start points, creating submatrix */
   d_st = start.i + start.j;
   // d_end = end.i + end.j;

   /* dimension of submatrix */
   dim_Q = Q - start.i;
   dim_T = T - start.j;

   /* diag index where num cells reaches highest point and begins diminishing */
   dim_min = min(Q + start.i, T + start.j);
   dim_max = max(Q + start.i, T + start.j);

   /* set bounds using starting cell */
   lb = start.i;
   rb = start.i + 1;
   num_cells = 0;

   /* INITIALIZE VALUES */
   /* initialize special states (?) */
   XMX(SP_N,0) = 0;                                         /* S->N, p=1             */
   XMX(SP_B,0) = XSC(SP_N,SP_MOVE);                         /* S->N->B, no N-tail    */
   XMX(SP_E,0) = XMX(SP_C,0) = XMX(SP_J,0) = -INF;          /* need seq to get here (?)  */

   /* initialize 0 row (top-edge) and ) column (left-edge) */
   for (j = 0; j <= T; j++)
      for (i = 0; i <= Q; i++)
         MMX(0,j) = IMX(0,j) = DMX(0,j) = -INF;

   /* keeps largest number seen on current diagonal */
   diag_max = -INF;
   total_max = -INF;

   /* ITERATE THROUGH ANTI-DIAGONALS */
   for (d = d_st; d <= d_end; d++, d_cnt++)
   {
      /* is dp matrix diagonal growing or shrinking? */
      if (d <= dim_min)
         num_cells++;
      if (d > dim_max)
         num_cells--;

      // printf("d=%d/%d, num_cells=%d\n", d, d_end, num_cells);

      /* UPDATE BOUNDS */
      /* if free passes are complete (beta < d), prune and set new edgebounds */
      if (beta < d_cnt)
      {
         lb_new = -100;
         rb_new = -100;

         /* Traverse current bounds to find max score on diag */
         diag_max = -INF;
         for (k = lb; k < rb; k++)
         {
            i = k;
            j = (d - 1) - i;    /* back one diag */
            diag_max = calc_Max( 
                           calc_Max( diag_max, MMX(i,j) ),
                           calc_Max( IMX(i,j), DMX(i,j) ) );
            // printf("diag_max(%d,%d) = %.2f \n", i,j,diag_max);
         }

         /* total max records largest cell score see so far */
         if (diag_max > total_max)
            total_max = diag_max;

         /* set score threshold for pruning */
         // diag_limit = diag_max - alpha;
         diag_limit = total_max - alpha;

         // printf("total_max: %.2f diag_max: %.2f diag_limit: %.2f\n", total_max, diag_max, diag_limit);

         /* Find the first cell from the left which passes above threshold */
         for (k = lb; k < rb; k++)
         {
            i = k;
            j = (d - 1) - i;
            cell_max = calc_Max( MMX(i,j),
                           calc_Max( IMX(i,j), DMX(i,j) ) );
            // printf("> TEST(%d,%d)=%.2f\n", i, j, cell_max);

            if( cell_max >= diag_limit )
            {
               lb_new = i;
               break;
            }
         }

         /* If no boundary edges are found on diag, then branch is pruned entirely and we are done */
         if (lb_new == -100)
            break;

         /* Find the first cell from the right which passes above threshold */
         for (k = rb - 1; k >= lb; k--)
         {
            i = k;
            j = (d - 1) - i;
            cell_max = calc_Max( MMX(i,j),
                           calc_Max( IMX(i,j),   DMX(i,j) ) );
            // printf("< TEST(%d,%d)=%.2f\n", i, j, cell_max);

            if( cell_max >= diag_limit )
            {
               rb_new = (i + 1);
               break;
            }
         }
      }
      else /* else edges expand in square pattern */
      {
         lb_new = lb;
         rb_new = rb;
      }

      /* Update bounds */
      lb = lb_new;
      rb = rb_new + 1;

      /* Edge-checks: find if diag cells that are inside matrix bounds */
      le = max(start.i, d - T);
      re = le + num_cells;

      /* Update bounds */
      if (lb < le)
         lb = le;
      if (rb > re)
         rb = re;

      // printf("%d: new[%d, %d] -> edg[%d, %d] -> [%d, %d]\n", d_cnt, lb_new, rb_new, le, re, lb, rb );

      /* ADD NEW ANTI-DIAG TO EDGEBOUNDS */
      edg->bounds[edg->N].lb = lb;
      edg->bounds[edg->N].rb = rb;
      edg->bounds[edg->N].diag = d;
      edg->N += 1;

      /* resize if needed */
      if (edg->N >= edg->size) {
         edg->size *= 2;
         edg->bounds = realloc(edg->bounds, edg->size * sizeof(BOUND) );
      }

      /* ITERATE THROUGH CELLS OF ANTI-DIAGONAL */
      for (k = lb; k < rb; k++, total_cnt++)
      {
         i = k;
         j = d - i;

         a = seq[i];
         A = AA_REV[a];

         /* MAIN RECURSION */

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

         /* SPECIAL STATES */
         /* J state */
         sc_1 = XMX(SP_J,i-1) + XSC(SP_J,SP_LOOP);       /* J->J */
         sc_2 = XMX(SP_E,i)   + XSC(SP_E,SP_LOOP);       /* E->J is E's "loop" */
         XMX(SP_J,i) = calc_Logsum( sc_1, sc_2 );         

         /* C state */
         sc_1 = XMX(SP_C,i-1) + XSC(SP_C,SP_LOOP);
         sc_2 = XMX(SP_E,i)   + XSC(SP_E,SP_MOVE);
         XMX(SP_C,i) = calc_Logsum( sc_1, sc_2 );

         /* N state */
         XMX(SP_N,i) = XMX(SP_N,i-1) + XSC(SP_N,SP_LOOP);

         /* B state */
         sc_1 = XMX(SP_N,i) + XSC(SP_N,SP_MOVE);         /* N->B is N's move */
         sc_2 = XMX(SP_J,i) + XSC(SP_J,SP_MOVE);         /* J->B is J's move */
         XMX(SP_B,i) = calc_Logsum( sc_1, sc_2 );   

         // printf("CALC (%d,%d): %.1f %.1f %.1f \n", i, j, MMX(i,j),IMX(i,j),DMX(i,j));
      }
   }
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
                        RESULTS* res,
                        EDGEBOUNDS* edg)
{
   /* vars for navigating matrix */
   int d,i,j,k;                  /* diagonal, row, column indices */
   int lb, rb, le, re;           /* right/left bounds and edges */
   int lb_new, rb_new;           /* tmp vars for new right/left bounds */
   int num_cells;                /* number of cells in diagonal */
   int d_st, d_end, d_cnt;       /* starting and ending diagonal indices */
   int dim_min, dim_max;         /* diagonal index where num cells reaches highest point and diminishing point */ 
   int dim_T, dim_Q;             /* dimensions of submatrix being searched */
   int total_cnt;                /* number of cells computed */
   COORDS start, end;            /* start and end point of alignment */

   /* vars for computing cells */
   char   a;                     /* store current character in sequence */
   int    A;                     /* store int value of character */
   char   *seq = query->seq;     /* alias for getting seq */

   /* vars for recurrance */
   float  prev_mat, prev_del, prev_ins, prev_beg, prev_end, prev_sum;
   float  sc, sc_1, sc_2, sc_best, sc_max;
   float  sc_M, sc_I, sc_D;
   /* vars for pruning */
   float  cell_max, total_max, diag_max, diag_limit;

   /* local or global? (multiple alignments) */
   bool   is_local = false;
   float  sc_E = (is_local) ? 0 : -INF;

   dp_matrix_Clear(Q, T, st_MX, sp_MX);
   d_cnt = 0;

   /* INITIALIZE VALUES */
   /* initialize special states (?) */
   XMX(SP_N,0) = 0;                                         /* S->N, p=1             */
   XMX(SP_B,0) = XSC(SP_N,SP_MOVE);                         /* S->N->B, no N-tail    */
   XMX(SP_E,0) = XMX(SP_C,0) = XMX(SP_J,0) = -INF;          /* need seq to get here (?)  */

   /* initialize 0 row (top-edge) and ) column (left-edge) */
   for (j = 0; j <= T; j++)
      for (i = 0; i <= Q; i++)
         MMX(0,j) = IMX(0,j) = DMX(0,j) = -INF;

   d_st = edg->bounds[0].diag;
   d_end = edg->bounds[edg->N-1].diag;

   /* ITERATE THROUGH ANTI-DIAGONALS */
   for (d = d_st; d <= d_end; d++, d_cnt++)
   {
      if (d < 2) {
         continue;
      }

      lb = max(1, edg->bounds[d_cnt].lb);
      rb = edg->bounds[d_cnt].rb;

      printf("%d: (%d,%d)\n", d,lb,rb);

      /* ITERATE THROUGH CELLS OF ANTI-DIAGONAL */
      for (k = lb; k < rb; k++, total_cnt++)
      {
         i = k;
         j = d - i;

         a = seq[i];
         A = AA_REV[a];

         /* MAIN RECURSION */

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

         /* SPECIAL STATES */
         /* J state */
         sc_1 = XMX(SP_J,i-1) + XSC(SP_J,SP_LOOP);       /* J->J */
         sc_2 = XMX(SP_E,i)   + XSC(SP_E,SP_LOOP);       /* E->J is E's "loop" */
         XMX(SP_J,i) = calc_Logsum( sc_1, sc_2 );         

         /* C state */
         sc_1 = XMX(SP_C,i-1) + XSC(SP_C,SP_LOOP);
         sc_2 = XMX(SP_E,i)   + XSC(SP_E,SP_MOVE);
         XMX(SP_C,i) = calc_Logsum( sc_1, sc_2 );

         /* N state */
         XMX(SP_N,i) = XMX(SP_N,i-1) + XSC(SP_N,SP_LOOP);

         /* B state */
         sc_1 = XMX(SP_N,i) + XSC(SP_N,SP_MOVE);         /* N->B is N's move */
         sc_2 = XMX(SP_J,i) + XSC(SP_J,SP_MOVE);         /* J->B is J's move */
         XMX(SP_B,i) = calc_Logsum( sc_1, sc_2 );   

         // printf("CALC (%d,%d): %.1f %.1f %.1f \n", i, j, MMX(i,j),IMX(i,j),DMX(i,j));
      }
   }

   /* T state */
   sc_best = XMX(SP_C,Q) + TSC(SP_C,SP_MOVE);

   // return sc_best;
}