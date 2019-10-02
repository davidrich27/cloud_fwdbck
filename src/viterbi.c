/*******************************************************************************
 *  @file viterbi.c
 *  @brief The Viterbi Algorithm for Sequence Alignment Search.
 *
 *  @author Dave Rich (devrek)
 *  @bug Lots.
 *******************************************************************************/

/* macros */
#define getName(var) #var
#define SCALE_FACTOR 1000

/* external imports */
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>

/* local imports (after struct declarations) */
#include "structs.h"
#include "misc.h"
#include "hmm_parser.h"
#include "viterbi.h"

/* 
 *  FUNCTION:  viterbi_Run()
 *  SYNOPSIS:  Run Viterbi Algorithm (Seq-to-Profile, general unoptimized)
 *
 *  PURPOSE:
 *
 *  ARGS:      <query>     query sequence, 
 *             <target>    HMM model,
 *             <res>       Results Data
 *             <tr>        Traceback Data
 *
 *  RETURN: 
 */
float viterbi_Run (const SEQ* query, 
                  const HMM_PROFILE* target, 
                  int Q, int T, 
                  float st_MX[ NUM_NORMAL_STATES * (Q+1) * (T+1) ], 
                  float sp_MX[ NUM_SPECIAL_STATES * (Q+1) ], 
                  RESULTS* res,
                  TRACEBACK* tr)
{
   char   a;           /* store current character in sequence */
   int    A;           /* store int value of character */
   int    i,j,k = 0;   /* row, column indices */
   char   *seq  = query->seq; /* alias for getting seq */

   float  prev_mat, prev_del, prev_ins, prev_beg, prev_best;
   float  sc, sc_1, sc_2, sc_best, sc_max;

   /* local or global (multiple alignments) */
   bool   is_local = false;
   float  sc_E = (is_local) ? 0 : -INF;

   /* initialize special states (?) */
   XMX(SP_N,0) = 0;                                         /* S->N, p=1             */
   XMX(SP_B,0) = XSC(SP_N,SP_MOVE);                         /* S->N->B, no N-tail    */
   XMX(SP_E,0) = XMX(SP_C,0) = XMX(SP_J,0) = -INF;          /* need seq to get here  */

   /* initialize zero row (top-edge) */
   for (j = 0; j <= T; j++)
      { MMX(0,j) = IMX(0,j) = DMX(0,j) = -INF; }            /* need seq to get here  */


   /* FOR every position in QUERY seq */
   for (i = 1; i <= Q; i++)
   {
      /* Get next character in Query */
      a = seq[i-1];
      A = AA_REV[a];

      /* Initialize zero column (left-edge) */
      MMX(i,0) = IMX(i,0) = DMX(i,0) = -INF;
      XMX(SP_E,i) = -INF;

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
         prev_best = calc_Max( calc_Max( prev_mat, prev_ins ) , calc_Max( prev_del, prev_beg ) );
         MMX(i,j)  = prev_best + MSC(j,A);

         /* FIND BEST PATH TO INSERT STATE (FROM MATCH OR INSERT) */
         /* previous states (match takes the left element of each state) */
         prev_mat = MMX(i-1,j) + TSC(j,M2I);
         prev_ins = IMX(i-1,j) + TSC(j,I2I);
         /* best-to-insert */
         prev_best = calc_Max(prev_mat, prev_ins);
         IMX(i,j) = prev_best + ISC(j,A);
         
         /* FIND BEST PATH TO DELETE STATE (FOMR MATCH OR DELETE) */
         /* previous states (match takes the left element of each state) */
         prev_mat = MMX(i,j-1) + TSC(j-1,M2D);
         prev_del = DMX(i,j-1) + TSC(j-1,D2D);
         /* best-to-delete */
         prev_best = calc_Max(prev_mat, prev_del);
         DMX(i,j) = prev_best;

         /* UPDATE E STATE */
         XMX(SP_E,i) = calc_Max( XMX(SP_E,i), MMX(i,j) + sc_E );

         /* check for best score seen so far (best score should end in a match state) */
         if (sc_max < MMX(i,j))
         {
            tr->sc_max = MMX(i,j);
            tr->end.i = i;
            tr->end.j = j;
         }
         
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
      prev_best = calc_Max( 
                     calc_Max( prev_mat, prev_ins ), 
                     calc_Max( prev_del, prev_beg ) 
                  );
      MMX(i,j) = prev_best + MSC(j,A);

      /* FIND BEST PATH TO DELETE STATE (FOMR MATCH OR DELETE) */
      /* previous states (match takes the left element of each state) */
      prev_mat = MMX(i,j-1) + TSC(j-1,M2D);
      prev_del = DMX(i,j-1) + TSC(j-1,D2D);
      /* best-to-delete */
      prev_best = calc_Max( prev_mat, prev_del );
      DMX(i,j) = prev_best;

      /* UPDATE E STATE */
      XMX(SP_E,i) = calc_Max( XMX(SP_E,i), MMX(i,j) );
      XMX(SP_E,i) = calc_Max( XMX(SP_E,i), DMX(i,j) );

      /* SPECIAL STATES */
      /* J state */
      sc_1 = XMX(SP_J,i-1) + XSC(SP_J,SP_LOOP);     /* J->J */
      sc_2 = XMX(SP_E,i)   + XSC(SP_E,SP_LOOP);     /* E->J is E's "loop" */
      XMX(SP_J,i) = calc_Max( sc_1, sc_2 );         

      /* C state */
      sc_1 = XMX(SP_C,i-1) + XSC(SP_C,SP_LOOP);  
      sc_2 = XMX(SP_E,i)   + XSC(SP_E,SP_MOVE);            
      XMX(SP_C,i) = calc_Max( sc_1, sc_2 );

      /* N state */
      XMX(SP_N,i) = XMX(SP_N,i-1) + XSC(SP_N,SP_LOOP);

      /* B state */
      sc_1 = XMX(SP_N,i) + XSC(SP_N,SP_MOVE);       /* N->B is N's move */
      sc_2 = XMX(SP_J,i) + XSC(SP_J,SP_MOVE);       /* J->B is J's move */
      XMX(SP_B,i) = calc_Max( sc_1, sc_2 );        
   }

   /* T state */
   sc_best = XMX(SP_C,Q) + XSC(SP_C,SP_MOVE);
   return sc_best;
}

/* 
 *  FUNCTION:  viterbi_Traceback()
 *  SYNOPSIS:  Run Viterbi Traceback to find Optimal Alignment
 *
 *  PURPOSE:
 *
 *  ARGS:      <query>     query sequence, 
 *             <target>    HMM model,
 *             <Q>         query length, 
 *             <T>         target length,
 *             <st_MX>     Normal State (Match, Insert, Delete) Matrix,
 *             <sp_MX>     Special State (J,N,B,C,E) Matrix,
 *             <res>       Results Data,
 *             <tr>        Traceback Alignment
 *
 *  RETURN: 
 */
void viterbi_Traceback (const SEQ* query, 
                        const HMM_PROFILE* target, 
                        int Q, int T, 
                        float st_MX[ NUM_NORMAL_STATES * (Q+1) * (T+1) ], 
                        float sp_MX[ NUM_SPECIAL_STATES * (Q+1) ], 
                        RESULTS* res,
                        TRACEBACK* tr)
{
   char   a;                  /* store current character in sequence */
   int    A;                  /* store int value of character */
   int    i,j,k = 0;          /* row, column indices */
   int    st_cur, st_prev;    /* current and previous state */
   int    len = 0;

   float  tol = 1e-5;         /* acceptable tolerance for equality */
   char   *seq = query->seq;  /* alias for getting seq */
   float  total = tr->sc_max; /* running total of alignment */

   bool is_local = false;

   /* allocate memory for trace */
   tr->i = malloc(1024 * sizeof(int));
   tr->j = malloc(1024 * sizeof(int));
   tr->trace = malloc(1024 * sizeof(int));

   i = tr->end.i;
   j = tr->end.j;

   /* final state of trace is match state */
   st_cur = -1;
   st_prev = C_ST;

   tr->i[len] = i;
   tr->j[len] = j;
   tr->trace[len] = st_prev;

   /* End of trace is S state */
   while(st_prev != S_ST)
   {
      a = seq[i-1];
      A = AA_REV[a];

      /* jump from current state to the prev state */
      switch(st_prev)
      {
         case M_ST:  /* M connects from i-1,k-1, or B */
            if (MMX(i,j) == -INF ) {
               printf("ERROR: Impossible M_ST reached at (%d,%d)\n", i,j);
               exit(1);
            }

            if ( CMP_TOL( MMX(i,j), MMX(i-1,j-1) + TSC(j-1,M2M) + MSC(j,A) ) ) 
               st_cur = M_ST;
            else if ( CMP_TOL( MMX(i,j), IMX(i-1,j-1) + TSC(j-1,I2M) + MSC(j,A) ) ) 
               st_cur = I_ST;
            else if ( CMP_TOL( MMX(i,j), MMX(i-1,j-1) + TSC(j-1,D2M) + MSC(j,A) ) ) 
               st_cur = D_ST;
            else if ( CMP_TOL( MMX(i,j), XMX(SP_B,j-1) + TSC(j-1,B2M) + MSC(j,A) ) ) 
               st_cur = B_ST;
            else {
               printf("ERROR: Failed to trace from M_ST at (%d,%d)\n", i, j);
               exit(1);
            }
            k--; i--;
            break;

         case I_ST:  /* I connects from M,I at i-1,k */
            if (IMX(i,j) == -INF ) {
               printf("ERROR: Impossible I_ST reached at (%d,%d)\n", i,j);
               exit(1);
            }

            if ( CMP_TOL( MMX(i,j), MMX(i-1,j) + TSC(j,M2I) + ISC(j,A) ) ) 
               st_cur = M_ST;
            else if ( CMP_TOL( MMX(i,j), IMX(i-1,j) + TSC(j,I2I) + ISC(j,A) ) ) 
               st_cur = I_ST;
            else {
               printf("ERROR: Failed to trace from I_ST at (%d,%d)\n", i, j);
               exit(1);
            }
            i--;
            break;

         case D_ST:  /* D connects from M,D at i,k-1 */
            if (DMX(i,j) == -INF ) {
               printf("ERROR: Impossible D_ST reached at (%d,%d)\n", i,j);
               exit(1);
            }

            if ( CMP_TOL( MMX(i,j), MMX(i-1,j-1) + TSC(j-1,M2M) ) ) 
               st_cur = M_ST;            
            else if ( CMP_TOL( MMX(i,j), MMX(i-1,j-1) + TSC(j-1,D2M) ) ) 
               st_cur = D_ST;
            else {
               printf("ERROR: Failed to trace from D_ST at (%d,%d)\n", i, j);
               exit(1);
            }
            k--;
            break;

         case E_ST:  /* E connects from any M state. k set here */
            if (XMX(SP_E,i) == -INF ) {
               printf("ERROR: Impossible E_ST reached at (%d,%d)\n", i,j);
               exit(1);
            }

            if ( is_local )  /* local mode: ends in M */
            {
               st_cur = M_ST;    /* can't come from D, in a *local* Viterbi trace. */
               for (j = T; j >= 1; j--) {
                  if ( CMP_TOL( XMX(SP_E, i), MMX(i,j) ) )
                     break;
               }
               if (j == 0) {
                  printf("ERROR: Failed to trace from E_ST at (%d,%d)\n", i, j);
                  exit(1);
               }
            }
            else     /* glocal mode: we either come from D_M or M_M */
            {
               if ( CMP_TOL( XMX(SP_E,i), MMX(i,T) ) ) {
                  st_cur = M_ST;
                  j = T;
               }
               else if ( CMP_TOL( XMX(SP_E,i), DMX(i,T) ) ) {
                  st_cur = D_ST;
                  j = T;
               }
               else {
                  printf("ERROR: Failed to trace from E_ST at (%d,%d)\n", i, j);
                  exit(1);
               }
            }
            break;

         case N_ST:  /* N connects from S, N */
            if (XMX(SP_N,i) == -INF ) {
               printf("ERROR: Impossible N_ST reached at (%d,%d)\n", i,j);
               exit(1);
            }

            st_cur = ( (i == 0) ? S_ST : N_ST );
            break;

         case J_ST:  /* J connects from E(i) or J(i-1) */
            if (XMX(SP_J,i) == -INF ) {
               printf("ERROR: Impossible J_ST reached at (%d,%d)\n", i,j);
               exit(1);
            }

            if ( CMP_TOL( XMX(SP_J,i), XMX(SP_J,i-1) + XSC(SP_J, SP_LOOP) ) )
               st_cur = J_ST;
            else if ( CMP_TOL( XMX(SP_J,i), XMX(SP_E,i) + XSC(SP_E, SP_LOOP) ) )
               st_cur = E_ST;
            else {
                  printf("ERROR: Failed to trace from J_ST at (%d,%d)\n", i, j);
                  exit(1);
            }
            break;

         case C_ST:  /* C(i) comes from C(i-1) or E(i) */
            if (XMX(SP_C,i) == -INF ) {
               printf("ERROR: Impossible C_ST reached at (%d,%d)\n", i,j);
               exit(1);
            }

            if ( CMP_TOL( XMX(SP_C,i), XMX(SP_C,i-1) + XSC(SP_C,SP_LOOP) ) )
               st_cur = N_ST;
            else if ( CMP_TOL( XMX(SP_C,i), XMX(SP_E,i) + XSC(SP_E,SP_MOVE) ) )
               st_cur = J_ST;
            else {
                  printf("ERROR: Failed to trace from B_ST at (%d,%d)\n", i, j);
                  exit(1);
            }
            break;

         case B_ST:  /* B connects from N, J */
            if ( CMP_TOL( XMX(SP_B,i), XMX(SP_N,i) + XSC(SP_N,SP_MOVE) ) )
               st_cur = N_ST;
            else if ( CMP_TOL( XMX(SP_B,i), XMX(SP_J,i) + XSC(SP_J,SP_MOVE) ) )
               st_cur = J_ST;
            else {
                  printf("ERROR: Failed to trace from B_ST at (%d,%d)\n", i, j);
                  exit(1);
            }
            break;

         default:
            printf("ERROR: Hit Bogus State!!!\n");
            tr->start.i = i;
            tr->start.j = j;
            return;
      }

      /* Add new state and (i,j) to trace */
      st_cur = st_prev;
      len++;
      tr->i[len] = i;
      tr->j[len] = j;
      tr->trace[len] = st_prev;

      /* For NCJ, defer i decrement. */
      if ( (st_cur == N_ST || st_cur == J_ST || st_cur == C_ST) && st_cur == st_prev) 
         i--;
   }

   return;
}
