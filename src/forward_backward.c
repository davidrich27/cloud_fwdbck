/*******************************************************************************
 *  @file forward_backward.c
 *  @brief The Forward-Backward Algorithm for Sequence Alignment Search.
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

/*  
 *  FUNCTION: forward_Run()
 *  SYNOPSIS: Perform Forward-Backward Algorithm (general, unoptimized)
 *
 *  PURPOSE:
 *
 *  ARGS:      <query>     query sequence, 
 *             <target>    HMM model,
 *             <res>       Results Data
 *
 *  RETURN: 
 */
void forward_backward_Run (const SEQ* query, 
                           const HMM_PROFILE* target,
                           int Q, int T, 
                           float st_MX[ NUM_NORMAL_STATES * (Q+1) * (T+1) ], 
                           float sp_MX[ NUM_SPECIAL_STATES * (Q+1) ], 
                           RESULTS* res)
{
   init_Logsum();
   // print_Logsum();

   /* allocate DP matrix */

   /* NORMAL STATE (MATCH, INSERT, DELETE) MATRIX */
   float st_MX2[ NUM_NORMAL_STATES * (Q+1) * (T+1) ]; 
   /* SPECIAL STATES MATRIX */
   float sp_MX2[ NUM_SPECIAL_STATES * (Q+1) ];

   /* NORMAL STATE (MATCH, INSERT, DELETE) MATRIX */
   float st_MX3[ NUM_NORMAL_STATES * (Q+1) * (T+1) ]; 
   /* SPECIAL STATES MATRIX */
   float sp_MX3[ NUM_SPECIAL_STATES * (Q+1) ];

   for (int i = 0; i <= Q; i++)
   {
      for (int j = 0; j <= T; j++)
      {
         ST_MX(st_MX,MAT_ST,i,j)  = ST_MX(st_MX,INS_ST,i,j)  = ST_MX(st_MX,DEL_ST,i,j)  = -INF;
         ST_MX(st_MX2,MAT_ST,i,j) = ST_MX(st_MX2,INS_ST,i,j) = ST_MX(st_MX2,DEL_ST,i,j) = -INF;
         ST_MX(st_MX3,MAT_ST,i,j) = ST_MX(st_MX3,INS_ST,i,j) = ST_MX(st_MX3,DEL_ST,i,j) = -INF;
      }

      for (int j = 0; j < NUM_SPECIAL_STATES; j++)
      {
         SP_MX(sp_MX,j,i) = -INF;
         SP_MX(sp_MX2,j,i) = -INF;
         SP_MX(sp_MX3,j,i) = -INF;
      }
   }

	forward_Run(query, target, Q, T, st_MX, sp_MX, res);
   backward_Run(query, target, Q, T, st_MX2, sp_MX2, res);

   int i,j;
   for (i = 0; i < Q; i++)
   {
      for (j = 0; j < T; j++)
      {
         ST_MX(st_MX3,MAT_ST,i,j) = ST_MX(st_MX,MAT_ST,i,j) + ST_MX(st_MX2,MAT_ST,i,j);
         ST_MX(st_MX3,INS_ST,i,j) = ST_MX(st_MX,INS_ST,i,j) + ST_MX(st_MX2,INS_ST,i,j);
         ST_MX(st_MX3,DEL_ST,i,j) = ST_MX(st_MX,DEL_ST,i,j) + ST_MX(st_MX2,DEL_ST,i,j);
      }
   }
}

/*  
 *  FUNCTION: forward_Run()
 *  SYNOPSIS: Perform Forward part of Forward-Backward Algorithm.
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
 *
 *  RETURN: 
 */
float forward_Run (const SEQ* query, 
                   const HMM_PROFILE* target, 
                   int Q, int T, 
                   float st_MX[ NUM_NORMAL_STATES * (Q+1) * (T+1) ], 
                   float sp_MX[ NUM_SPECIAL_STATES * (Q+1) ], 
                   RESULTS* res)
{
   char   a;           /* store current character in sequence */
   int    A;           /* store int value of character */
   int    i,j,k = 0;   /* row, column indices */
   char   *seq = query->seq; /* alias for getting seq */

   float  prev_mat, prev_del, prev_ins, prev_beg, prev_sum;
   float  sc, sc_1, sc_2;
   float  sc_max, sc_best;
   COORDS tr_end; /* ending match state of optimal alignment (for traceback) */

   /* local or global (multiple alignments) */
   bool   is_local = true;
   float  sc_E = (is_local) ? 0 : -INF;

   /* initialize matrix */
   for (i = 0; i <= Q; i++)
   {
      for (j = 0; j <= T; j++)
      {
         MMX(i,j) = IMX(i,j) = DMX(i,j) = -INF;
      }
   }

   /* initialize special states (?) */
   XMX(SP_N,0) = 0;                                         /* S->N, p=1             */
   XMX(SP_B,0) = XSC(SP_N,SP_MOVE);                         /* S->N->B, no N-tail    */
   XMX(SP_E,0) = XMX(SP_C,0) = XMX(SP_J,0) = -INF;          /* need seq to get here (?)  */

   /* initialize 0 row (top-edge) */
   for (j = 0; j < T; j++)
      { MMX(0,j) = IMX(0,j) = DMX(0,j) = -INF; }       /* need seq to get here (?)  */

   /* FOR every position in QUERY seq */
   for (i = 1; i <= Q; i++)
   {  
      /* Get next sequence character */
      a = seq[i-1];
      A = AA_REV[a];

      /* Initialize zero column (left-edge) */
      MMX(i,0) = IMX(i,0) = DMX(i,0) = -INF;
      XMX(SP_E,i) = -INF;

      /* MAIN RECURSION */
      /* FOR every position in TARGET profile */
      for (j = 1; j < T; j++)
      {
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

         /* SAVE FOR DEBUGGING (transition states slightly off?) */
         // printf("(%d, %d): ", i, j);
         // printf("Mp: %.2f, Mtr: %.2f, ", MMX(i-1,j-1),TSC(j-1,M2M));
         // printf("Ip: %.2f, Itr: %.2f, ", IMX(i-1,j-1),TSC(j-1,I2M));
         // printf("Dp: %.2f, Dtr: %.2f, ", DMX(i-1,j-1),TSC(j-1,D2M));
         // printf("Bp: %.2f, Btr: %.2f, ", XMX(SP_B,i-1),TSC(j-1,B2M));
         // printf("sum: %.2f, Mc: %.2f\n", prev_sum, MMX(i,j));

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

      /* UNROLLED FINAL LOOP ITERATION */
      j = T; 

      /* FIND BEST PATH TO MATCH STATE (FROM MATCH, INSERT, DELETE, OR BEGIN) */
      /* best previous state transition (match takes the diag element of each prev state) */
      prev_mat = MMX(i-1,j-1)  + TSC(j-1,M2M);
      prev_ins = IMX(i-1,j-1)  + TSC(j-1,I2M);
      prev_del = DMX(i-1,j-1)  + TSC(j-1,D2M);
      prev_beg = XMX(SP_B,i-1) + TSC(j-1,B2M);    /* from begin match state (new alignment) */
      /* best-to-match */
      prev_sum = calc_Logsum( 
                        calc_Logsum( prev_mat, prev_ins ),
                        calc_Logsum( prev_del, prev_beg )
                  );
      MMX(i,j) = prev_sum + MSC(j,A);

      /* FIND BEST PATH TO INSERT STATE */
      IMX(i,j) = -INF;

      /* FIND BEST PATH TO DELETE STATE (FOMR MATCH OR DELETE) */
      /* previous states (match takes the left element of each state) */
      prev_mat = MMX(i,j-1) + TSC(j-1,M2D);
      prev_del = DMX(i,j-1) + TSC(j-1,D2D);
      /* best-to-delete */
      prev_sum = calc_Logsum( prev_mat, prev_del );
      DMX(i,j) = prev_sum;

      /* UPDATE E STATE */
      XMX(SP_E,i) = calc_Logsum( XMX(SP_E,i), MMX(i,j) );
      XMX(SP_E,i) = calc_Logsum( XMX(SP_E,i), DMX(i,j) );

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
   }

   /* T state */
   sc_best = XMX(SP_C,Q) + TSC(SP_C,SP_MOVE);

   return sc_best;
}

/* FUNCTION: backward_Run()
 * SYNOPSIS: Perform Backward part of Forward-Backward Algorithm.
 *
 * PURPOSE:
 *
 *  ARGS:      <query>     query sequence, 
 *             <target>    HMM model,
 *             <Q>         query length, 
 *             <T>         target length,
 *             <st_MX>     Normal State (Match, Insert, Delete) Matrix,
 *             <sp_MX>     Special State (J,N,B,C,E) Matrix,
 *             <res>       Results Data
 *
 * RETURN: 
*/
float backward_Run (const SEQ* query, 
                    const HMM_PROFILE* target, 
                    int Q, int T, 
                    float st_MX[ NUM_NORMAL_STATES * (Q+1) * (T+1) ], 
                    float sp_MX[ NUM_SPECIAL_STATES * (Q+1) ], 
                    RESULTS* res)
{
   char   a;           /* store current character in sequence */
   int    A;           /* store int value of character */
   int    i,j,k = 0;   /* row, column indices */
   char   *seq = query->seq; /* alias for getting seq */

   float  prev_mat, prev_del, prev_ins, prev_end, prev_sum;
   float  sc, sc_1, sc_2, sc_M, sc_I;
   float  sc_max, sc_best;

   /* local or global (multiple alignments) */
   bool   is_local = true;
   float  sc_E = (is_local) ? 0 : -INF;

   /* Initialize the Q row. */
   XMX(SP_J,Q) = XMX(SP_B,Q) = XMX(SP_N,Q) = -INF;
   XMX(SP_C,Q) = XSC(SP_C,SP_MOVE);
   XMX(SP_E,Q) = XMX(SP_C,Q) + XSC(SP_E,SP_MOVE);

   printf("(%d) C: %.2f E: %.2f Cm: %.2f Em: %.2f\n", Q, XMX(SP_C,Q), XMX(SP_E,Q), XSC(SP_C,SP_MOVE), XSC(SP_E,SP_MOVE));

   MMX(Q,T) = DMX(Q,T) = XMX(SP_E,Q);
   IMX(Q,T) = -INF;

   printf("(%d,%d) M: %.2f D: %.2f I: %2.f\n", Q,T, MMX(Q,T), DMX(Q,T), IMX(Q,T) );

   for (j = T-1; j >= 1; j--)
   {
      MMX(Q,j) = calc_Logsum( XMX(SP_E,Q) + sc_E, 
                              DMX(Q,j+1)  + TSC(j,M2D) );
      DMX(Q,j) = calc_Logsum( XMX(SP_E,Q) + sc_E,
                              DMX(Q,j+1)  + TSC(j,D2D) );
      IMX(Q,j) = -INF;

      printf("(%d,%d) M: %.2f D: %.2f I: %2.f, ", Q,j, MMX(Q,j), DMX(Q,j), IMX(Q,j) );
      printf("E: %.2f Dp: %.2f Tmd: %.2f Tdd: %.2f\n", XMX(SP_E,Q), DMX(Q,j+1), TSC(j,M2D), TSC(j,D2D) );
   }

   /* MAIN RECURSION */
   /* FOR every position in QUERY seq */
   for (i = Q-1; i >= 1; i--)
   {
      /* Get next sequence character */
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

      MMX(i,T) = DMX(i,T) = XMX(SP_E,i);
      IMX(i,T) = -INF;

      /* FOR every position in TARGET profile */
      for (j = T-1; j >= 1; j--)
      {
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
                           calc_Logsum( prev_del, prev_end )
                     );
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

         /* SAVE FOR DEBUGGING (transition states slightly off?) */
         printf("(%d, %d): ", i, j);
         printf("Mp: %.2f, Mtr: %.2f Msc: %2f, ", MMX(i+1,j+1), TSC(j,M2M), sc_M);
         printf("Ip: %.2f, Itr: %.2f Isc: %.2f, ", IMX(i+1,j), TSC(j,M2I), sc_I);
         printf("Dp: %.2f, Dtr: %.2f, ", DMX(i,j+1), TSC(j,M2D));
         printf("Ep: %.2f, Esc: %.2f, ", XMX(SP_E,i), sc_E);
         printf("sum: %.2f, Mc: %.2f\n", prev_sum, MMX(i,j));
      }
   }

   /* FINAL i = 0 row */
   a = seq[1];
   A = AA_REV[a];

   XMX(SP_B,0) = MMX(1,1) + TSC(0,B2M) + MSC(1,A);

   for (j = 2; j >= T; j++)
   {
      XMX(SP_B,0) = calc_Logsum( XMX(SP_B,0),
                                 MMX(1,j) + TSC(j-1,B2M) + MSC(j,A) );
   }

   XMX(SP_J,i) = -INF;
   XMX(SP_C,i) = -INF;
   XMX(SP_E,i) = -INF;

   XMX(SP_N,i) = calc_Logsum( XMX(SP_N,1) + XSC(SP_N,SP_LOOP),
                              XMX(SP_B,0) + XSC(SP_N,SP_MOVE) );

   for (j = T; j >= 1; j--)
   {
      MMX(i,j) = IMX(i,j) = DMX(i,j) = -INF;
   }

   sc_best = XMX(SP_N,0);

   return sc_best;
}








