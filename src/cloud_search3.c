/*******************************************************************************
 *  @file cloud_search3.c
 *  @brief Cloud Forward-Backward Algorithm (Linear Space Implementation)
 *
 *  @synopsis
 *       NOTE: HOW TO CONVERT row-coords to diag-coords
 *       MMX(i-1,j-1) => MMX3(, d_2)
 *       MMX(i,  j-1) => MMX3(, d_1)
 *       MMX(i,  j  ) => MMX3(, d_1)
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
#include "cloud_search3.h"

// macros
// #define getName(var) #var
// #define SCALE_FACTOR 1000

// macro functions
// NOTE: wrap all macro vars in parens!!
#define max(x,y) (((x) > (y)) ? (x) : (y))
#define min(x,y) (((x) < (y)) ? (x) : (y))


/*
 *  FUNCTION: cloud_forward3_Run()
 *  SYNOPSIS: Perform Forward part of Cloud Search Algorithm (Linear Space Implementation).
 *
 *  PURPOSE:
 *
 *  ARGS:      <query>     query sequence,
 *             <target>    HMM model,
 *             <Q>         query length,
 *             <T>         target length,
 *             <st_MX3>    Normal State (Match, Insert, Delete) 3x (linear space),
 *             <sp_MX>     Special State (J,N,B,C,E) Matrix,
 *             <res>       Results Data
 *             <tr>        Traceback Data
 *             <edg>       Edge Bounds Tracker Data
 *          free variables:
 *             <alpha>     Pruning Drop
 *             <beta>      Number of Passes before
 *
 *  RETURN:
 */
void cloud_forward_Run3 (const SEQ* query,
                         const HMM_PROFILE* target,
                         int Q, int T,
                         float* st_MX3,
                         float* sp_MX,
                         TRACEBACK* tr,
                         EDGEBOUNDS* edg,
                         float alpha, int beta,
                         int *sc_final )
{
   /* vars for navigating matrix */
   int d, i, j, k;               /* diagonal, row, column indices */
   int lb, rb, le, re;           /* right/left bounds and edges */
   int lb_new, rb_new;           /* tmp vars for new right/left bounds */
   int num_cells;                /* number of cells in diagonal */
   int d_st, d_end, d_cnt;       /* starting and ending diagonal indices */
   int dim_min, dim_max;         /* diagonal index where num cells reaches highest point and diminishing point */
   int dim_T, dim_Q;             /* dimensions of submatrix being searched */
   int total_cnt;                /* number of cells computed */
   COORDS start, end;            /* start and end point of alignment */
   BOUND old_b, new_b;           /* bounds for clearing old data */

   /* vars for computing cells */
   char   a;                     /* store current character in sequence */
   int    A;                     /* store int value of character */
   char   *seq = query->seq;     /* alias for getting seq */

   /* vars for recurrance */
   int    d_mod3, d_0, d_1, d_2; /* d (mod 3) for assigning prev array ptrs */
   float  prev_mat, prev_del, prev_ins, prev_beg, prev_end, prev_sum;
   float  sc, sc_1, sc_2, sc_best, sc_max;
   float  sc_M, sc_I, sc_D;
   /* vars for pruning */
   float  cell_max, total_max, diag_max, diag_limit;

   /* local or global? (multiple alignments) */
   bool   is_local = target->isLocal;
   float  sc_E = (is_local) ? 0 : -INF;

   /* INIT EDGEBOUND DATA STRUCT */
   int min_size = 128;
   edg->size = min_size;
   edg->bounds = (BOUND *)malloc( min_size * sizeof(BOUND) );

   /* --------------------------------------------------------------------------------- */

   // dp_matrix_Clear(Q, T, st_MX, sp_MX);

   /* INIT ANTI-DIAGS */
   /* test coords */
   start = tr->first_m;
   end = tr->last_m;

   /* TESTING */
   // printf("T=%d, Q=%d, alpha=%.1f, beta=%d\n", T, Q, alpha, beta);
   // alpha = 6.0;
   // beta = 5;
   // start = (COORDS) {1, 1};
   // end = (COORDS) {T, Q};

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

   // /* INITIALIZE VALUES */
   // /* initialize special states (?) */
   // XMX(SP_N,0) = 0;                                         /* S->N, p=1             */
   // XMX(SP_B,0) = XSC(SP_N,SP_MOVE);                         /* S->N->B, no N-tail    */
   // XMX(SP_E,0) = XMX(SP_C,0) = XMX(SP_J,0) = -INF;          /* need seq to get here (?)  */

   /* initialize zero row */
   for (j = 0; j <= 3; j++)
      for (i = 0; i <= Q; i++)
         MMX3(i, j) = IMX3(i, j) = DMX3(i, j) = -INF;

   /* keeps largest number seen on current diagonal */
   diag_max = -INF;
   total_max = -INF;

   /* begin state probability begins at 0 */
   prev_beg = 0;

   /* ITERATE THROUGH ANTI-DIAGONALS */
   for (d = d_st; d <= d_end; d++, d_cnt++)
   {
      d_0 = d % 3;      /* current diagonal */
      d_1 = (d - 1) % 3; /* look back 1 diagonal */
      d_2 = (d - 2) % 3; /* look back 2 diagonals */

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

         /* Traverse current bounds to find max score on prev diag */
         diag_max = -INF;
         for (k = lb; k < rb; k++)
         {
            i = k;
            /* back one diag */
            diag_max = calc_Max(
                          calc_Max( diag_max,    MMX3(k, d_1) ),
                          calc_Max( IMX3(k, d_1), DMX3(k, d_1) ) );
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
            /* looking back one diag */
            cell_max = calc_Max( MMX3(k, d_1),
                                 calc_Max( IMX3(k, d_1), DMX3(k, d_1) ) );
            // printf("> TEST(%d,%d)=%.2f\n", i, j, cell_max);

            if ( cell_max >= diag_limit )
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
            cell_max = calc_Max( MMX3(k, d_1),
                                 calc_Max( IMX3(k, d_1),   DMX3(k, d_1) ) );
            // printf("< TEST(%d,%d)=%.2f\n", i, j, cell_max);

            if ( cell_max >= diag_limit )
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

      printf("%d: new[%d, %d] -> edg[%d, %d] -> [%d, %d]\n", d_cnt, lb_new, rb_new, le, re, lb, rb );

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

      /* SCRUB OLD DATA OF NEW CURRENT ANTI-DIAG */
      if (d_cnt >= 2) {
         old_b = edg->bounds[(edg->N - 1) - 2]; /* look back-3 diags to find bounds of data to be scrubbed */
         for (k = old_b.lb; k < old_b.rb; k++)
         {
            MMX3(k, d_0) = IMX3(k, d_0) = DMX3(k, d_0) = -INF;
         }
      }


      /* ITERATE THROUGH CELLS OF ANTI-DIAGONAL */
      for (k = lb; k < rb; k++, total_cnt++)
      {
         i = k;
         j = d - i;

         a = seq[i];
         A = AA_REV[a];

         /* MAIN RECURSION */

         /* FIND SUM OF PATHS TO MATCH STATE (FROM MATCH, INSERT, DELETE, OR BEGIN) */
         /* best previous state transition (match takes the diag element of each prev state) */
         prev_mat = MMX3(k - 1, d_2)  + TSC(j - 1, M2M); /* MMX(i-1,j-1) */
         prev_ins = IMX3(k - 1, d_2)  + TSC(j - 1, I2M); /* MMX(i-1,j-1) */
         prev_del = DMX3(k - 1, d_2)  + TSC(j - 1, D2M); /* MMX(i-1,j-1) */
         /* from begin match state (new alignment) */
         // prev_beg = XMX(SP_B,i-1) + TSC(j-1,B2M);
         prev_beg = TSC(j - 1, B2M);
         /* best-to-match */
         prev_sum = calc_Logsum(
                       calc_Logsum( prev_mat, prev_ins ),
                       calc_Logsum( prev_del, prev_beg ) );
         MMX3(k, d_0) = prev_sum + MSC(j, A);      /* MMX(i,j) */

         /* FIND SUM OF PATHS TO INSERT STATE (FROM MATCH OR INSERT) */
         /* previous states (match takes the left element of each state) */
         prev_mat = MMX3(k, d_1) + TSC(j, M2I);  /* MMX(i-1,j) */
         prev_ins = IMX3(k, d_1) + TSC(j, I2I);  /* IMX(i-1,j) */
         /* best-to-insert */
         prev_sum = calc_Logsum( prev_mat, prev_ins );
         IMX3(k, d_0) = prev_sum + ISC(j, A);    /* IMX(i,j) */

         /* FIND SUM OF PATHS TO DELETE STATE (FOMR MATCH OR DELETE) */
         /* previous states (match takes the left element of each state) */
         prev_mat = MMX3(k - 1, d_1) + TSC(j - 1, M2D); /* MMX(i,j-1) */
         prev_del = DMX3(k - 1, d_1) + TSC(j - 1, D2D); /* DMX(i,j-1) */
         /* best-to-delete */
         prev_sum = calc_Logsum(prev_mat, prev_del);
         DMX3(k, d_0) = prev_sum;                  /* DMX(i,j) */

         /* UPDATE E STATE */
         // XMX(SP_E,i) = calc_Logsum( XMX(SP_E,i),
         //                            MMX(k,d_0) + sc_E );
         // XMX(SP_E,i) = calc_Logsum( XMX(SP_E,i),
         //                            DMX(k,d_0) + sc_E );

         // /* SPECIAL STATES */
         // /* J state */
         // sc_1 = XMX(SP_J,i-1) + XSC(SP_J,SP_LOOP);       /* J->J */
         // sc_2 = XMX(SP_E,i)   + XSC(SP_E,SP_LOOP);       /* E->J is E's "loop" */
         // XMX(SP_J,i) = calc_Logsum( sc_1, sc_2 );

         // /* C state */
         // sc_1 = XMX(SP_C,i-1) + XSC(SP_C,SP_LOOP);
         // sc_2 = XMX(SP_E,i)   + XSC(SP_E,SP_MOVE);
         // XMX(SP_C,i) = calc_Logsum( sc_1, sc_2 );

         // /* N state */
         // XMX(SP_N,i) = XMX(SP_N,i-1) + XSC(SP_N,SP_LOOP);

         // /* B state */
         // sc_1 = XMX(SP_N,i) + XSC(SP_N,SP_MOVE);         /* N->B is N's move */
         // sc_2 = XMX(SP_J,i) + XSC(SP_J,SP_MOVE);         /* J->B is J's move */
         // XMX(SP_B,i) = calc_Logsum( sc_1, sc_2 );

         // printf("CALC (%d,%d): %.1f %.1f %.1f \n", i, j, MMX(i,j),IMX(i,j),DMX(i,j));
      }
   }
}


/*
 *  FUNCTION: cloud_backward_Run3()
 *  SYNOPSIS: Perform Backward part of Cloud Search Algorithm (Linear Space Implementation).
 *
 *  PURPOSE:
 *
 *  ARGS:      <query>     query sequence,
 *             <target>    HMM model,
 *             <Q>         query length,
 *             <T>         target length,
 *             <st_MX3>    Normal State (Match, Insert, Delete) 3x (linear space),
 *             <sp_MX>     Special State (J,N,B,C,E) Matrix,
 *             <res>       Results Data
 *             <tr>        Traceback Data
 *             <edg>       Edge Bounds Tracker Data
 *             <alpha>     Pruning Drop
 *             <beta>      Number of Passes before
 *
 *  RETURN:
 */
void cloud_backward_Run3 (const SEQ* query,
                          const HMM_PROFILE* target,
                          int Q, int T,
                          float* st_MX3,
                          float* sp_MX,
                          TRACEBACK* tr,
                          EDGEBOUNDS* edg,
                          float alpha, int beta,
                          int *sc_final )
{
   /* vars for navigating matrix */
   int d, i, j, k, b;            /* diagonal, row, column, ... indices */
   int lb, rb, le, re;           /* right/left bounds and edges */
   int lb_new, rb_new;           /* tmp vars for new right/left bounds */
   int num_cells;                /* number of cells in diagonal */
   int d_st, d_end, d_cnt;       /* starting and ending diagonal indices */
   int dim_min, dim_max;         /* diagonal index where num cells reaches highest point and diminishing point */
   int dim_T, dim_Q;             /* dimensions of submatrix being searched */
   int total_cnt;                /* number of cells computed */
   COORDS start, end;            /* start and end point of alignment */
   BOUND old_b, new_b;           /* bounds for clearing old data */

   /* vars for computing cells */
   char   a;                     /* store current character in sequence */
   int    A;                     /* store int value of character */
   char   *seq = query->seq;     /* alias for getting seq */

   /* vars for recurrance */
   int    d_mod3, d_0, d_1, d_2; /* d (mod 3) for assigning prev array ptrs */
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
   int min_size = 128;
   edg->size = min_size;
   edg->bounds = (BOUND *)malloc( min_size * sizeof(BOUND) );

   /* --------------------------------------------------------------------------------- */

   // dp_matrix_Clear(Q, T, st_MX, sp_MX);

   /* INIT ANTI-DIAGS */
   start = tr->first_m;
   end = tr->last_m;

   /* TESTING */
   // alpha = 6.0;
   // beta = 5;
   // start = (COORDS) {1, 1};
   // end = (COORDS) {T, Q};
   // printf("T=%d, Q=%d, alpha=%.1f, beta=%d\n", T, Q, alpha, beta);

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

   /* ADD INITIAL ANTI-DIAG TO EDGEBOUNDS */
   edg->bounds[edg->N].lb = lb;
   edg->bounds[edg->N].rb = rb;
   edg->bounds[edg->N].diag = d;
   edg->N += 1;
   /* resize if needed */
   if (edg->N >= edg->size) {
      edg->size *= 2;
      edg->bounds = realloc(edg->bounds, edg->size * sizeof(BOUND) );
   }

   // /* INITIALIZE SPECIAL STATES */
   // /* Initialize the Q row. */
   // XMX(SP_J,Q) = XMX(SP_B,Q) = XMX(SP_N,Q) = -INF;
   // XMX(SP_C,Q) = XSC(SP_C,SP_MOVE);
   // XMX(SP_E,Q) = XMX(SP_C,Q) + XSC(SP_E,SP_MOVE);

   // MMX(Q,T) = DMX(Q,T) = XMX(SP_E,Q);
   // IMX(Q,T) = -INF;

   // for (j = T-1; j >= 1; j--)
   // {
   //    MMX(Q,j) = calc_Logsum( XMX(SP_E,Q) + sc_E,
   //                            DMX(Q,j+1)  + TSC(j,M2D) );
   //    DMX(Q,j) = calc_Logsum( XMX(SP_E,Q) + sc_E,
   //                            DMX(Q,j+1)  + TSC(j,D2D) );
   //    IMX(Q,j) = -INF;
   // }

   /* initialize zero row */
   for (j = 0; j <= 3; j++)
      for (i = 0; i <= Q; i++)
         MMX3(i, j) = IMX3(i, j) = DMX3(i, j) = -INF;

   // /* compute special states probs up to first i-index */
   // for (i = Q-1; i >= end.i; i--) {

   //    /* Push from Match State to Next */
   //    a = seq[i];
   //    A = AA_REV[a];

   //    /* SPECIAL STATES */
   //    XMX(SP_B,i) = MMX(i+1,1) + TSC(0,B2M) + MSC(1,A);

   //    /* B -> MATCH */
   //    for (j = 2; j <= T; j++)
   //    {
   //       XMX(SP_B,i) = calc_Logsum( XMX(SP_B,i),
   //                                  MMX(i+1,j) + TSC(j-1,B2M) + MSC(j,A) );
   //    }

   //    XMX(SP_J,i) = calc_Logsum( XMX(SP_J,i+1) + XSC(SP_J,SP_LOOP),
   //                               XMX(SP_B,i)   + XSC(SP_J,SP_MOVE) );

   //    XMX(SP_C,i) = XMX(SP_C,i+1) + XSC(SP_C,SP_LOOP);

   //    XMX(SP_E,i) = calc_Logsum( XMX(SP_J,i) + XSC(SP_E,SP_LOOP),
   //                               XMX(SP_C,i) + XSC(SP_E,SP_MOVE) );

   //    XMX(SP_N,i) = calc_Logsum( XMX(SP_N,i+1) + XSC(SP_N,SP_LOOP),
   //                               XMX(SP_B,i)   + XSC(SP_N,SP_MOVE) );
   // }

   /* end state starts at 0 */
   prev_end = 0;

   /* ITERATE THROUGHT ANTI-DIAGONALS */
   for (d = d_end; d >= d_st; d--, d_cnt++)
   {
      d_0 = d % 3;      /* current diagonal */
      d_1 = (d + 1) % 3; /* look back 1 diagonal */
      d_2 = (d + 2) % 3; /* look back 2 diagonals */

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

         /* Traverse current bounds to find max score on prev diag */
         diag_max = -INF;
         for (k = lb; k < rb; k++)
         {
            i = k;
            /* back one diag */
            diag_max = calc_Max(
                          calc_Max( diag_max,    MMX3(k, d_1) ),
                          calc_Max( IMX3(k, d_1), DMX3(k, d_1) ) );
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
            /* looking back one diag */
            cell_max = calc_Max( MMX3(k, d_1),
                                 calc_Max( IMX3(k, d_1), DMX3(k, d_1) ) );
            // printf("> TEST(%d,%d)=%.2f\n", i, j, cell_max);

            if ( cell_max >= diag_limit )
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
            cell_max = calc_Max( MMX3(k, d_1),
                                 calc_Max( IMX3(k, d_1),   DMX3(k, d_1) ) );
            // printf("< TEST(%d,%d)=%.2f\n", i, j, cell_max);

            if ( cell_max >= diag_limit )
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

      /* SCRUB OLD DATA OF NEW CURRENT ANTI-DIAG */
      if (d_cnt >= 2) {
         old_b = edg->bounds[(edg->N - 1) - 2]; /* look back 3 diags to find bounds of data to be scrubbed */
         for (k = old_b.lb; k < old_b.rb; k++)
         {
            MMX3(k, d_0) = IMX3(k, d_0) = DMX3(k, d_0) = -INF;
         }
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

         // /* SPECIAL STATES */
         // XMX(SP_B,i) = MMX(i+1,1) + TSC(0,B2M) + MSC(1,A);

         // /* B -> MATCH */
         // for (b = 2; b <= T; b++)
         // {
         //    XMX(SP_B,i) = calc_Logsum( XMX(SP_B,i),
         //                               MMX(i+1,j) + TSC(j-1,B2M) + MSC(j,A) );
         // }

         // XMX(SP_J,i) = calc_Logsum( XMX(SP_J,i+1) + XSC(SP_J,SP_LOOP),
         //                            XMX(SP_B,i)   + XSC(SP_J,SP_MOVE) );

         // XMX(SP_C,i) = XMX(SP_C,i+1) + XSC(SP_C,SP_LOOP);

         // XMX(SP_E,i) = calc_Logsum( XMX(SP_J,i) + XSC(SP_E,SP_LOOP),
         //                            XMX(SP_C,i) + XSC(SP_E,SP_MOVE) );

         // XMX(SP_N,i) = calc_Logsum( XMX(SP_N,i+1) + XSC(SP_N,SP_LOOP),
         //                         XMX(SP_B,i)   + XSC(SP_N,SP_MOVE) );
         // MMX(i,T) = DMX(i,T) = XMX(SP_E,i);
         // IMX(i,T) = -INF;

         /* MAIN RECURSION */
         sc_M = MSC(j + 1, A);
         sc_I = ISC(j + 1, A);

         /* FIND SUM OF PATHS FROM MATCH, INSERT, DELETE, OR END STATE (TO PREVIOUS MATCH) */
         prev_mat = MMX3(i + 1, j + 1) + TSC(j, M2M) + sc_M;
         prev_ins = IMX3(i + 1, j)   + TSC(j, M2I) + sc_I;
         prev_del = DMX3(i, j + 1)   + TSC(j, M2D);
         // prev_end = XMX(SP_E,i)  + sc_E;     /* from end match state (new alignment) */
         // prev_end = sc_E;
         /* best-to-match */
         prev_sum = calc_Logsum(
                       calc_Logsum( prev_mat, prev_ins ),
                       calc_Logsum( prev_del, prev_end ) );
         MMX3(i, j) = prev_sum;

         /* FIND SUM OF PATHS FROM MATCH OR INSERT STATE (TO PREVIOUS INSERT) */
         prev_mat = MMX3(i + 1, j + 1) + TSC(j, I2M) + sc_M;
         prev_ins = IMX3(i + 1, j)   + TSC(j, I2I) + sc_I;
         /* best-to-insert */
         prev_sum = calc_Logsum( prev_mat, prev_ins );
         IMX3(i, j) = prev_sum;

         /* FIND SUM OF PATHS FROM MATCH OR DELETE STATE (FROM PREVIOUS DELETE) */
         prev_mat = MMX3(i + 1, j + 1) + TSC(j, D2M) + sc_M;
         prev_del = DMX3(i, j + 1)   + TSC(j, D2D);
         prev_end = XMX(SP_E, i)  + sc_E;
         /* best-to-delete */
         prev_sum = calc_Logsum( prev_mat, prev_del );
         prev_sum = calc_Logsum( prev_sum, prev_end );
         DMX3(i, j) = prev_sum;

         // printf("CALC %d=>(%d,%d): %.1f %.1f %.1f \n", d, i, j, MMX(i,j),IMX(i,j),DMX(i,j));
      }
   }

   /* reverse order of diagonals */
   BOUND tmp;
   for (i = 0; i < (edg->N / 2); ++i)
   {
      tmp = edg->bounds[i];
      edg->bounds[i] = edg->bounds[edg->N - i];
      edg->bounds[edg->N - i] = tmp;
   }
}



/*
 *  FUNCTION: forward_bounded_Run()
 *  SYNOPSIS: Perform Edge-Bounded Forward part of Cloud Search Algorithm.
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
float forward_bounded_Run3(const SEQ* query,
                           const HMM_PROFILE* target,
                           int Q, int T,
                           float* st_MX3,
                           float* sp_MX,
                           EDGEBOUNDS* edg, 
                           int *sc_final )
{
   char   a;                              /* store current character in sequence */
   int    A;                              /* store int value of character */
   int    i,j,k = 0;                      /* row, column indices */
   char   *seq = query->seq;              /* alias for getting seq */
   int    N = edg->N;                     /* length of edgebound list */

   int    x, y1, y2;                      /* row, leftcol and rightcol bounds in row */
   int    x_0, row_cur, x_1, row_prv;     /* real index of current and previous rows */
   int    r_0, r_1;                       /* row offset -> r_0: row_cur % 2, r_1: row_prv % 2 */
   int    r_0b, r_0e, r_1b, r_1e;         /* begin and end indices for row in edgebound list */

   float  prev_mat, prev_del, prev_ins;   /* temp placeholder sums */
   float  prev_beg, prev_end, prev_sum;   /* temp placeholder sums */
   float  sc, sc1, sc2, sc3, sc4;         /* temp placeholder sums (testing) */
   float  sc_best;                        /* alignment score (return value) */
   float  sc_M, sc_I, sc_D;

   /* local or global? (multiple alignments) */
   bool   is_local = target->isLocal;
   float  sc_E = (is_local) ? 0 : -INF;

   /* --------------------------------------------------------------------------------- */

   /* initialize zero row / clear all pre-existing data from matrix */
   for (j = 0; j <= 3; j++)
      for (i = 0; i <= Q; i++)
         MMX3(i, j) = IMX3(i, j) = DMX3(i, j) = -INF;

   /* initialize special states (?) */
   XMX(SP_N,0) = 0;                                         /* S->N, p=1             */
   XMX(SP_B,0) = XSC(SP_N,SP_MOVE);                         /* S->N->B, no N-tail    */
   XMX(SP_E,0) = XMX(SP_C,0) = XMX(SP_J,0) = -INF;          /* need seq to get here (?)  */

   /* initialize zero row (top-edge) */
   for (j = 0; j < T; j++)
      MMX3(0,j) = IMX3(0,j) = DMX3(0,j) = -INF;              /* need seq to get here (?)  */

   /* pass over top-row (i=0) edgebounds from list */
   row_cur = r_0 = 0;   /* current row in matrix */
   k = 0;               /* current index in edgebounds */
   r_0b = 0;
   while ( k < N && edg->bounds[k].diag == row_cur ) {
      k++;
   }
   r_0e = k;

   /* init look back 1 (r_1) */
   x_1 = 0;
   r_1 = r_0;
   r_1b = r_0b;
   r_1e = r_0e;

   /* MAIN RECURSION */
   /* FOR every position in QUERY sequence (row in matrix) */
   for (x_0 = 1; x_0 <= Q; x_0++)
   {
      /* convert quadratic space row index to linear space row index (ex % 2) */
      row_cur = x_0;
      r_0 = x_0;        /* for use in linear space alg (mod-mapping) */
      r_1 = x_0 - 1;    /* for use in linear space alg (mod-mapping) */

      /* add every edgebound from current row */
      r_0b = k;
      // printf("k: %d => row_cur: %d, k.diag: %d\n", k, row_cur, edg->bounds[k].diag);
      while ( k < N && edg->bounds[k].diag == row_cur ) {
         k++;
      }
      r_0e = k;

      /* Get next sequence character */
      a = seq[x_1];  /* off-by-one */
      A = AA_REV[a];

      /* Initialize zero column (left-edge) */
      MMX3(r_0, 0) = IMX3(r_0, 0) = DMX3(r_0, 0) = -INF;
      XMX(SP_E, x_0) = -INF;

      /* FOR every EDGEBOUND in current ROW */
      for (i = r_0b; i < r_0e; i++)
      {
         /* in this context, "diag" represents the "row" */
         x = edg->bounds[i].diag;            /* NOTE: this is always the same as cur_row, x_0 */
         y1 = max(1, edg->bounds[i].lb);     /* can't overflow the left edge */
         y2 = min(edg->bounds[i].rb, T);   /* can't overflow the right edge */
         // printf("%d: (%d->%d,%d->%d)\n", x, edg->bounds[i].lb, y1, edg->bounds[i].rb, y2);

         /* MAIN RECURSION */
         /* FOR every position in TARGET profile */
         for (j = y1; j < y2; j++)
         {
            /* FIND SUM OF PATHS TO MATCH STATE (FROM MATCH, INSERT, DELETE, OR BEGIN) */
            /* best previous state transition (match takes the diag element of each prev state) */
            prev_mat = MMX3(r_1,j-1)  + TSC(j-1,M2M);
            prev_ins = IMX3(r_1,j-1)  + TSC(j-1,I2M);
            prev_del = DMX3(r_1,j-1)  + TSC(j-1,D2M);
            prev_beg = XMX(SP_B,i-1) + TSC(j-1,B2M); /* from begin match state (new alignment) */

            /* best-to-match */
            prev_sum = calc_Logsum( 
                           calc_Logsum( prev_mat, prev_ins ),
                           calc_Logsum( prev_del, prev_beg )
                        );
            MMX3(r_0,j) = prev_sum + MSC(j,A);
            // printf("MMX(%d,%d): %f\t MSC(%d,%d): %f\n", i, j, MMX(r_1,j), j, A, MSC(j,A) );

            /* FIND SUM OF PATHS TO INSERT STATE (FROM MATCH OR INSERT) */
            /* previous states (match takes the left element of each state) */
            prev_mat = MMX3(r_1,j) + TSC(j,M2I);
            prev_ins = IMX3(r_1,j) + TSC(j,I2I);
            /* best-to-insert */
            prev_sum = calc_Logsum( prev_mat, prev_ins );
            IMX3(r_0,j) = prev_sum + ISC(j,A);

            /* FIND SUM OF PATHS TO DELETE STATE (FROM MATCH OR DELETE) */
            /* previous states (match takes the left element of each state) */
            prev_mat = MMX3(r_0,j-1) + TSC(j-1,M2D);
            prev_del = DMX3(r_0,j-1) + TSC(j-1,D2D);
            /* best-to-delete */
            prev_sum = calc_Logsum(prev_mat, prev_del);
            DMX3(r_0,j) = prev_sum;

            /* UPDATE E STATE */
            prev_mat = MMX3(r_0, j) + sc_E;
            prev_del = DMX3(r_0, j) + sc_E;
            XMX(SP_E, x_0) = calc_Logsum( 
                                    calc_Logsum( prev_mat, prev_del ),
                                    XMX(SP_E, x_0) );

            // printf("(%d,%d) -> E: %f, Esc: %f, MMX: %f, DMX: %f, MSC: %f \n", x_0, j, XMX(SP_E, x_0), sc_E, MMX(r_0, j), DMX(r_0, j), MSC(j,A) );
         }

         /* UNROLLED FINAL LOOP ITERATION */
         j = y2; 

         /* FIND SUM OF PATHS TO MATCH STATE (FROM MATCH, INSERT, DELETE, OR BEGIN) */
         /* best previous state transition (match takes the diag element of each prev state) */
         prev_mat = MMX3(r_1,j-1)  + TSC(j-1,M2M);
         prev_ins = IMX3(r_1,j-1)  + TSC(j-1,I2M);
         prev_del = DMX3(r_1,j-1)  + TSC(j-1,D2M);
         prev_beg = XMX(SP_B,r_1) + TSC(j-1,B2M);    /* from begin match state (new alignment) */
         /* best-to-match */
         prev_sum = calc_Logsum( 
                           calc_Logsum( prev_mat, prev_ins ),
                           calc_Logsum( prev_del, prev_beg )
                     );
         MMX3(r_0,j) = prev_sum + MSC(j,A);

         /* FIND SUM OF PATHS TO INSERT STATE */
         IMX3(r_0,j) = -INF;

         /* FIND SUM OF PATHS TO DELETE STATE (FROM MATCH OR DELETE) */
         /* previous states (match takes the left element of each state) */
         prev_mat = MMX3(i,j-1) + TSC(j-1,M2D);
         prev_del = DMX3(i,j-1) + TSC(j-1,D2D);
         /* best-to-delete */
         prev_sum = calc_Logsum( prev_mat, prev_del );
         DMX3(r_0,j) = prev_sum;

         /* UPDATE E STATE */
         prev_mat = MMX3(r_0, j);
         prev_del = DMX3(r_0, j);
         XMX(SP_E, x_0) = calc_Logsum( 
                              calc_Logsum( prev_mat, prev_del ),
                              XMX(SP_E, x_0)
                           );
      }

      /* ONCE ROW IS COMPLETED, UPDATE SPECIAL STATES */
      {
         /* SPECIAL STATES */
         /* J state */
         sc1 = XMX(SP_J, x_1) + XSC(SP_J, SP_LOOP);   /* J->J */
         sc2 = XMX(SP_E, x_0) + XSC(SP_E, SP_LOOP);   /* E->J is E's "loop" */
         XMX(SP_J, x_0) = calc_Logsum( sc1, sc2 );

         /* C state */
         sc1 = XMX(SP_C, x_1) + XSC(SP_C, SP_LOOP);
         sc2 = XMX(SP_E, x_0) + XSC(SP_E, SP_MOVE);
         XMX(SP_C, x_0) = calc_Logsum( sc1, sc2 );
         // printf("x_0: %d -> C(-1): %f, Cloop: %f, E(0): %f, Emove: %f\n", x_0, XMX(SP_C, x_1), XSC(SP_C, SP_LOOP), XMX(SP_E, x_0), XSC(SP_E, SP_MOVE) );

         /* N state */
         XMX(SP_N, x_0) = XMX(SP_N, x_1) + XSC(SP_N, SP_LOOP);

         /* B state */
         sc1 = XMX(SP_N, x_0) + XSC(SP_N, SP_MOVE);       /* N->B is N's move */
         sc2 = XMX(SP_J, x_0) + XSC(SP_J, SP_MOVE);       /* J->B is J's move */
         XMX(SP_B, x_0) = calc_Logsum( sc1, sc2 );

         // printf("x_0: %d -> J: %f, C: %f, N: %f, B: %f\n", x_0, XMX(SP_J, x_0), XMX(SP_C, x_0), XMX(SP_N, x_0), XMX(SP_N, x_0));
      }

      /* SCRUB PREVIOUS ROW VALUES */
      for (i = r_1b; i < r_1e; i++) 
      {
         /* in this context, "diag" represents the "row" */
         // x  = edg->bounds[i].diag;
         y1 = max(0, edg->bounds[i].lb);     /* can't overflow the left edge */
         y2 = min(edg->bounds[i].rb, T);     /* can't overflow the right edge */

         for (j = y1; j <= y2; j++) {
            MMX3(r_1, j) = IMX3(r_1, j) = DMX3(r_1, j) = -INF;
         }
      }

      /* SET CURRENT ROW TO PREVIOUS ROW */
      row_prv = x_1 = row_cur;
      r_1 = r_0;
      r_1b = r_0b;
      r_1e = r_0e;
   }

   /* T state */
   sc_best = XMX(SP_C,Q) + XSC(SP_C,SP_MOVE);
   printf("sc_best: %f\n", sc_best);

   sc_final[0] = sc_best;
   return sc_best;
}


/*
 *  FUNCTION: backward_bounded_Run()
 *  SYNOPSIS: Perform Edge-Bounded Backward part of Cloud Search Algorithm.
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
float backward_bounded_Run3(const SEQ* query,
                            const HMM_PROFILE* target,
                            int Q, int T,
                            float* st_MX3,
                            float* sp_MX,
                            EDGEBOUNDS* edg, 
                            int *sc_final )
{
   char   a;                              /* store current character in sequence */
   int    A;                              /* store int value of character */
   int    i,j,k = 0;                      /* row, column indices */
   char   *seq = query->seq;              /* alias for getting seq */
   int    N = edg->N;                     /* length of edgebound list */

   int    x, y1, y2;                      /* row, leftcol and rightcol bounds in row */
   int    x_0, row_cur, x_1, row_prv;     /* real index of current and previous rows */
   int    r_0, r_1;                       /* row offset -> r_0: row_cur % 2, r_1: row_prv % 2 */
   int    r_0b, r_0e, r_1b, r_1e;         /* begin and end indices for row in edgebound list */

   float  prev_mat, prev_del, prev_ins;   /* temp placeholder sums */
   float  prev_beg, prev_end, prev_sum;   /* temp placeholder sums */
   float  sc, sc1, sc2, sc3, sc4;         /* temp placeholder sums (testing) */
   float  sc_best;                        /* alignment score (return value) */
   float  sc_M, sc_I, sc_D;

   /* local or global? (multiple alignments) */
   bool   is_local = target->isLocal;
   float  sc_E = (is_local) ? 0 : -INF;

   /* --------------------------------------------------------------------------------- */

   /* initialize zero row / clear all pre-existing data from matrix */
   for (j = 0; j <= 3; j++)
      for (i = 0; i <= Q; i++)
         MMX3(i, j) = IMX3(i, j) = DMX3(i, j) = -INF;

   /* Initialize the Q row. */
   XMX(SP_J,Q) = XMX(SP_B,Q) = XMX(SP_N,Q) = -INF;
   XMX(SP_C,Q) = XSC(SP_C,SP_MOVE);
   XMX(SP_E,Q) = XMX(SP_C,Q) + XSC(SP_E,SP_MOVE);

   MMX3(r_0,T) = DMX3(r_0,T) = XMX(SP_E,Q);
   IMX3(r_0,T) = -INF;

   for (j = T-1; j >= 1; j--)
   {
      sc1 = XMX(SP_E,Q) + sc_E;
      sc2 = DMX3(r_0,j+1)  + TSC(j,M2D);
      MMX3(r_0,j) = calc_Logsum( XMX(SP_E,Q) + sc_E, 
                              DMX3(r_0,j+1)  + TSC(j,M2D) );

      sc1 = XMX(SP_E,Q) + sc_E;
      sc2 = DMX3(r_0,j+1)  + TSC(j,D2D);
      DMX3(r_0,j) = calc_Logsum( XMX(SP_E,Q) + sc_E,
                              DMX3(r_0,j+1)  + TSC(j,D2D) );
      
      IMX3(r_0,j) = -INF;
   }

   /* pass over (Q) bottom-row edgebounds from list */
   row_cur = r_0 = Q;
   k = N-1;
   r_0b = N-1;
   while ( k >= 0 && edg->bounds[k].diag == row_cur ) {
      k--;
   }
   r_0e = k;

   x_1 = Q;    /* initial one-row-back */
   r_1 = x_1;  /* for use in linear space alg (mod-mapping) */

   /* MAIN RECURSION */
   /* FOR every bound in EDGEBOUND */
   for (x_0 = Q-1; x_0 > 0; --x_0)
   {
      /* convert quadratic space row index to linear space row index (ex % 2) */
      row_cur = x_0;
      r_0 = x_0 % 2;          /* for use in linear space alg (mod-mapping) */
      r_1 = (x_0 + 1) % 2;    /* for use in linear space alg (mod-mapping) */

      /* add every edgebound from current row */
      r_0b = k;
      // printf("k: %d => row_cur: %d, k.diag: %d\n", k, row_cur, edg->bounds[k].diag);
      while ( k >= 0 && edg->bounds[k].diag >= row_cur ) {
         k--;
      }
      r_0e = k;

      /* Get next sequence character */
      a = seq[x_0];
      A = AA_REV[a];

      /* UPDATE SPECIAL STATES at the start of EACH ROW */
      {
         /* SPECIAL STATES */
         j = 1;
         XMX(SP_B, x_0) = MMX3(r_1, j) + TSC(j-1, B2M) + MSC(j, A);
         // printf("B(%d)\t (%d)%f\t", x_0, j, XMX(SP_B, x_0) );

         /* B -> MATCH */
         for (j = 2; j <= T; j++) {
            XMX(SP_B, x_0) = calc_Logsum( XMX(SP_B, x_0),
                                       MMX3(r_1, j) + TSC(j-1, B2M) + MSC(j, A) );
            // printf("B(%d,%d): B=%f, M=%f, TSC=%f, MSC=%f\n", x_0, j, XMX(SP_B, x_0), MMX(r_1, j), TSC(j-1, B2M), MSC(j, A) );
         }
         // printf("\n");

         XMX(SP_J, x_0) = calc_Logsum( XMX(SP_J, x_1) + XSC(SP_J, SP_LOOP),
                                       XMX(SP_B, x_0) + XSC(SP_J, SP_MOVE) );

         XMX(SP_C, x_0) = XMX(SP_C, x_1) + XSC(SP_C, SP_LOOP);

         XMX(SP_E, x_0) = calc_Logsum( XMX(SP_J, x_0) + XSC(SP_E,SP_LOOP),
                                       XMX(SP_C, x_0) + XSC(SP_E,SP_MOVE) );

         XMX(SP_N, x_0) = calc_Logsum( XMX(SP_N, x_1) + XSC(SP_N,SP_LOOP),
                                       XMX(SP_B, x_0) + XSC(SP_N,SP_MOVE) );

         MMX3(r_0, T) = DMX3(r_0, T) = XMX(SP_E, x_0);
         IMX3(r_0, T) = -INF;

         // printf("x_0: %d -> B: %f, J: %f, C: %f, E: %f, N: %f, MSC(%d): %f \n", x_0, XMX(SP_B, x_0), XMX(SP_J, x_0), XMX(SP_C, x_0),  XMX(SP_E, x_0), XMX(SP_N, x_0), A, MSC(j, A) );
      }
      // printf("x_0: %d (%d, %d) -> (%d, %d)\n", x_0, r_0b, r_0e, edg->bounds[r_0b], edg->bounds[r_0e]);

      /* FOR every EDGEBOUND in current ROW */
      for (i = r_0b; i > r_0e; i--)
      {
         /* in this context, "diag" represents the "row" */
         // x_0  = edg->bounds[i].diag;
         y1 = max(1, edg->bounds[i].lb);     /* can't overflow the left edge */
         y2 = min(edg->bounds[i].rb, T);     /* can't overflow the right edge */

         /* FOR every position in TARGET profile */
         for (j = y2-1; j >= y1; --j)
         {
            sc_M = MSC(j+1,A);
            sc_I = ISC(j+1,A);
            // printf("(%d,%d): A=%d, MSC=%f, ISC=%f\n", x_0, j, A, sc_M, sc_I);

            /* FIND SUM OF PATHS FROM MATCH, INSERT, DELETE, OR END STATE (TO PREVIOUS MATCH) */
            prev_mat = MMX3(r_1, j+1)  + TSC(j, M2M) + sc_M;
            prev_ins = IMX3(r_1, j)    + TSC(j, M2I) + sc_I;
            prev_del = DMX3(r_0, j+1)  + TSC(j, M2D);
            prev_end = XMX(SP_E, r_0) + sc_E;     /* from end match state (new alignment) */
            /* best-to-match */
            prev_sum = calc_Logsum( 
                              calc_Logsum( prev_mat, prev_ins ),
                              calc_Logsum( prev_end, prev_del ) );
            MMX3(r_0, j) = prev_sum;

            /* FIND SUM OF PATHS FROM MATCH OR INSERT STATE (TO PREVIOUS INSERT) */
            prev_mat = MMX3(r_1, j+1) + TSC(j, I2M) + sc_M;
            prev_ins = IMX3(r_1, j)   + TSC(j, I2I) + sc_I;
            /* best-to-insert */
            prev_sum = calc_Logsum( prev_mat, prev_ins );
            IMX3(r_0,j) = prev_sum;

            /* FIND SUM OF PATHS FROM MATCH OR DELETE STATE (FROM PREVIOUS DELETE) */
            prev_mat = MMX3(r_1, j+1)  + TSC(j,D2M) + sc_M;
            prev_del = DMX3(r_0, j+1)  + TSC(j,D2D);
            prev_end = XMX(SP_E, x_0) + sc_E;
            /* best-to-delete */
            prev_sum = calc_Logsum( 
                              prev_mat,
                              calc_Logsum( prev_del, prev_end ) );
            DMX3(r_0,j) = prev_sum;

            // printf("(%d,%d): M=%f, I=%f, D=%f\n", r_0, j, MMX(r_0, j), IMX(r_0,j), DMX(r_0,j) );
         }
      }

      /* SCRUB PREVIOUS ROW VALUES */
      for (i = r_1b; i > r_1e; i--) 
      {
         /* in this context, "diag" represents the "row" */
         // x  = edg->bounds[i].diag;
         y1 = max(0, edg->bounds[i].lb);     /* can't overflow the left edge */
         y2 = min(edg->bounds[i].rb, T);     /* can't overflow the right edge */

         for (j = y1-1; j >= y2; --j) {
            MMX3(r_1, j) = IMX3(r_1, j) = DMX3(r_1, j) = -INF;
         }
      }

      /* SET CURRENT ROW TO PREVIOUS ROW */
      row_prv = x_1 = row_cur;
      r_1 = r_0;
      r_1b = r_0b;
      r_1e = r_0e;
   }

   /* FINAL ROW */
   row_cur = x_0 = 1;
   r_0 = row_cur % 2;
   /* FINAL i = 0 row */
   a = seq[x_0];
   A = AA_REV[a];

   XMX(SP_B,0) = MMX3(1,1) + TSC(0,B2M) + MSC(1,A);

   for (j = 2; j >= T; j++) {
      XMX(SP_B,0) = calc_Logsum( XMX(SP_B,0),
                                 MMX3(1,j) + TSC(j-1,B2M) + MSC(j,A) );
   }

   XMX(SP_J,i) = -INF;
   XMX(SP_C,i) = -INF;
   XMX(SP_E,i) = -INF;

   XMX(SP_N,i) = calc_Logsum( XMX(SP_N,1) + XSC(SP_N,SP_LOOP),
                              XMX(SP_B,0) + XSC(SP_N,SP_MOVE) );

   for (j = T; j >= 1; j--) {
      MMX3(i,j) = IMX3(i,j) = DMX3(i,j) = -INF;
   }

   sc_best = XMX(SP_N,0);

   sc_final[0] = sc_best;
   return sc_best;
}
