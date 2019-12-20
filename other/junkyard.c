/* code pulled from the project */

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
 *             <edg>       Edge Bounds Tracker Data
 *             <alpha>     Pruning Drop 
 *             <beta>      Number of Passes before
 *
 *  RETURN: 
 */
void cloud_search_forward_Run(const SEQ* query, 
                              const HMM_PROFILE* target,
                              int Q, int T, 
                              float st_MX[ NUM_NORMAL_STATES * (Q+1) * (T+1) ], 
                              float sp_MX[ NUM_SPECIAL_STATES * (Q+1) ], 
                              RESULTS* res,
                              TRACEBACK* tr,
                              EDGEBOUNDS* edg,
                              float alpha, int beta)
{
   printf("cloud forward search...\n");

   char   a;                     /* store current character in sequence */
   int    A;                     /* store int value of character */

   int    d,i,j,k;               /* diagonal, row, column indices */
   char   *seq = query->seq;     /* alias for getting seq */
   float  diag_max, diag_limit;  /* max prob score in the diag, and the pruning floor */
   int    d_cnt = 0;             /* number of anti-diags from starting position */
   int    left_bound, right_bound;
   int    left_new, right_new; 
   int    left_edge, right_edge;
   float  cell_max, total_max;
   int    num_cells;

   /* local or global (multiple alignments) */
   bool   is_local = false;
   float  sc_E = (is_local) ? 0 : -INF;

   /* these are struct arrays => { diag_num, lb, rb } */
   int min_size = 128;
   edg->size = min_size;
   edg->bounds = (BOUND *)malloc( min_size * sizeof(BOUND) );

   /* recurrance vars */
   float  prev_mat, prev_del, prev_ins, prev_beg, prev_end, prev_sum;
   float  sc, sc_1, sc_2, sc_best, sc_max;
   float  sc_M, sc_I, sc_D;

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
      for (i = 0; i <= Q; i++)
         MMX(0,j) = IMX(0,j) = DMX(0,j) = -INF;

   /* TODO: temp testing var */
   tr->first_m.i = 1;
   tr->first_m.j = 1;

   /* get starting anti-diag (d = i + j) */
   d = tr->first_m.i + tr->first_m.j;

   /* position on anti-diag equal to the ith column, unless d is greater than Q (we've past bottom left corner) */
   left_bound = right_bound = tr->first_m.i;
   if (d > Q) {
      left_bound -= (d - Q);
      right_bound -= (d - Q);
   }

   /* keeps largest number seen on current diagonal */
   diag_max = -INF;
   total_max = -INF;

   /* number of cells in current diagonal */
   if (d < min(Q,T)) {
      num_cells = d;
   }
   else if (d < max(Q,T)) {
      num_cells = min(Q,T);
   } 
   else {
      num_cells = Q + T - d;
   }

   /* compute special states probs up to first diag */
   for (i = 0; i < d; i++) {
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
   }

   /* MAIN DIAG RECURSION */

   /* increment by anti-diagonal */
   for (; d <= Q+T; d++)
   {
      d_cnt++;
      /* number of cells expanding or tapering? */
      if (d < Q && d < T)
         ++num_cells;
      if (d > Q && d > T)
         --num_cells;

      // printf("d= %d/%d, lb: %d, rb: %d\n", d, Q+T, left_bound, right_bound);

      /* TODO: Update Prune Bounds */
      /* if free passes are over (beta < d), prune and set new edgebounds */
      if (beta < d_cnt)
      {
         left_new = -100;
         right_new = -100;
         /* Traverse current bounds to find max score on diag */
         for (k = left_bound; k <= right_bound; k++)
         {
            i = k;
            j = d - i - 1; /* back one diag */
            diag_max = calc_Max( 
                           calc_Max( diag_max, MMX(i,j) ),
                           calc_Max( IMX(i,j), DMX(i,j) ) );
         }

         /* total max records largest cell score see so far */
         if (diag_max > total_max)
            total_max = diag_max;

         // diag_limit = diag_max - alpha;
         diag_limit = total_max - alpha;
         printf("total_max: %.2f diag_max: %.2f diag_limit: %.2f\n", total_max, diag_max, diag_limit);

         /* Find the first cell from the left which passes above threshold */
         for (k = left_bound; k <= right_bound; k++)
         {
            i = k;
            j = d - i - 1;
            cell_max = calc_Max( MMX(i,j),
                           calc_Max( IMX(i,j), DMX(i,j) ) );
            if( cell_max >= diag_limit )
            {
               left_new = i - 1;
               break;
            }
         }

         /* If no boundary edges are found on diag, then branch is pruned entirely and we are done */
         if (left_new == -100)
            break;

         /* Find the first cell from the right which passes above threshold */
         for (k = right_bound; k >= left_bound; k--)
         {
            i = k;
            j = d - i - 1;
            cell_max = calc_Max( MMX(i,j),
                           calc_Max( IMX(i,j),   DMX(i,j) ) );
            if( cell_max >= diag_limit )
            {
               right_new = i + 1;
               break;
            }
         }
      }
      else /* else edges expand in square pattern */
      {
         left_new = left_bound;
         right_new = right_bound+1;

         if (d > Q) {
            left_new -= 1;
            right_new -= 1;
         }
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

      /* add new bounds to edgebound tracker */
      edg->N += 1;
      edg->bounds[edg->N].lb = left_bound;
      edg->bounds[edg->N].rb = right_bound;
      edg->bounds[edg->N].diag = d;

      if (edg->N >= edg->size) {
         edg->size *= 2;
         edg->bounds = realloc(edg->bounds, edg->size * sizeof(BOUND) );
      }

      printf("d= %d/%d, lb: %d, rb: %d, ", d, Q+T, left_bound, right_bound);
      printf("le: %d, re: %d\n", left_edge, right_edge);

      /* left to right across anti-diagonal */
      for (k = left_bound; k <= right_bound; k++)
      {
         i = k;
         j = d - i; 

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

         printf("CALC (%d,%d): %.1f %.1f %.1f \n", i, j, MMX(i,j),IMX(i,j),DMX(i,j));
      }
   }

   return;
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
 *             <edg>       Edge Bounds Tracker Data
 *             <alpha>     Pruning Drop 
 *             <beta>      Number of Passes before
 *
 *  RETURN: 
 */
void cloud_search_backward_Run(const SEQ* query, 
                              const HMM_PROFILE* target,
                              int Q, int T, 
                              float st_MX[ NUM_NORMAL_STATES * (Q+1) * (T+1) ], 
                              float sp_MX[ NUM_SPECIAL_STATES * (Q+1) ], 
                              RESULTS* res,
                              TRACEBACK* tr,
                              EDGEBOUNDS* edg,
                              float alpha, 
                              int beta)
{
   printf("cloud backward search...\n");

   char   a;                     /* store current character in sequence */
   int    A;                     /* store int value of character */
   int    d,i,j,k;               /* diagonal, row, column indices */
   char   *seq = query->seq;     /* alias for getting seq */
   float  diag_max, diag_limit;  /* max prob score in the diag, and the pruning floor */
   int    d_cnt = 0;             /* number of anti-diags from starting position */
   int    left_bound, right_bound;
   int    left_new, right_new; 
   int    left_edge, right_edge;
   float  cell_max, total_max;
   int    num_cells;

   /* local or global (multiple alignments) */
   bool   is_local = false;
   float  sc_E = (is_local) ? 0 : -INF;

   /* these are struct arrays => { diag_num, lb, rb } */
   int min_size = 128;
   edg->size = min_size;
   edg->bounds = (BOUND *)malloc( min_size * sizeof(BOUND) );

   /* recurrance vars */
   float  prev_mat, prev_del, prev_ins, prev_beg, prev_end, prev_sum;
   float  sc, sc_1, sc_2, sc_best, sc_max;
   float  sc_M, sc_I, sc_D;

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
      for (i = 0; i <= Q; i++)
         MMX(0,j) = IMX(0,j) = DMX(0,j) = -INF;

   /* TODO: temp testing var */
   tr->last_m.i = Q;
   tr->last_m.j = T;

   /* get starting anti-diag (d = i - j) */
   d = tr->last_m.i + tr->last_m.j;

   /* position on anti-diag equal to the ith column, unless d is greater than Q (we've past bottom left corner) */
   left_bound = right_bound = tr->last_m.i;
   // if (d > T) {
   //    left_bound = (d - T);
   //    right_bound = (d - T);
   // }
   printf("bounds: %d <=> %d, i=%d, j=%d, d=%d, Q=%d, T=%d\n", left_bound, right_bound, tr->last_m.i, tr->last_m.j, d, Q, T);

   /* keeps largest number seen on current diagonal */
   diag_max = -INF;
   total_max = -INF;

   /* number of cells in current diagonal */
   if (d < min(Q,T)) {
      num_cells = d;
   }
   else if (d < max(Q,T)) {
      num_cells = min(Q,T);
   } 
   else {
      num_cells = Q + T - d;
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

   /* compute special states probs up to first diag */
   for (i = Q; i > d; i--) {
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

   /* MAIN DIAG RECURSION */

   /* decrement by anti-diagonal */
   for (; d >= 0; d--)
   {
      d_cnt++;

      /* number of cells expanding or tapering? */
      if (d < Q && d < T)
         --num_cells;
      if (d > Q && d > T)
         ++num_cells;

      // printf("d= %d/%d, lb: %d, rb: %d\n", d, Q+T, left_bound, right_bound);

      /* TODO: Update Prune Bounds */
      /* if free passes are over (beta < d_cnt), prune and set new edgebounds */
      if (beta < d_cnt)
      {
         left_new = -1;
         right_new = -1;
         /* Traverse current bounds to find max score on diag */
         for (k = left_bound; k <= right_bound; k++)
         {
            i = k;
            j = d - i - 1; /* back one diag */
            diag_max = calc_Max( 
                           calc_Max( diag_max, MMX(i,j) ),
                           calc_Max( IMX(i,j), DMX(i,j) ) );
         }

         /* total max records largest cell score see so far */
         if (diag_max > total_max)
            total_max = diag_max;

         // diag_limit = diag_max - alpha;
         diag_limit = total_max - alpha;
         printf("total_max: %.2f diag_max: %.2f diag_limit: %.2f\n", total_max, diag_max, diag_limit);

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

         /* If no boundary edges are found on diag, then branch is pruned entirely and we are done */
         if (left_new == -1)
            break;

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
         left_new = left_bound-1;
         right_new = right_bound;

         /* if matrix is tapering, diag indices shifts right 1 */
         // if (d < T) {
         //    left_new -= 1;
         //    right_new -= 1;
         // }
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
      /* TODO: FIX THIS!!! */
      right_bound = calc_Max(left_edge, right_bound);

      /* add new bounds to edgebound tracker */
      edg->N += 1;
      edg->bounds[edg->N].lb = left_bound;
      edg->bounds[edg->N].rb = right_bound;
      edg->bounds[edg->N].diag = d;

      if (edg->N >= edg->size) {
         edg->size *= 2;
         edg->bounds = realloc(edg->bounds, edg->size * sizeof(BOUND) );
      }

      printf("d= %d/%d, lb: %d, rb: %d, ", d, Q+T, left_bound, right_bound);
      printf("le: %d, re: %d\n", left_edge, right_edge);

      /* left to right across anti-diagonal */
      for (k = left_bound; k <= right_bound; k++)
      {
         i = k;
         j = d - i; 
         printf("CALC (%d,%d) \n", i, j);

         /* Push from Match State to Next */
         a = seq[i];
         A = AA_REV[a];
         
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

         printf("CALC (%d,%d): %.1f %.1f %.1f \n", i, j, MMX(i,j),IMX(i,j),DMX(i,j));
      }
   }

   return;
}
