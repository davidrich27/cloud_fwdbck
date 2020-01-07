/*******************************************************************************
 *  @file forward_backward.c
 *  @brief Functions for EDGEBOUNDS object.
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
 *  FUNCTION: edgebounds_Create()
 *  SYNOPSIS: Create new edgebounds object.
 *
 *  PURPOSE:
 *
 *  ARGS:      <edg>      Edgebounds Object
 *
 *  RETURN:
 */
void edgebounds_Create(EDGEBOUNDS *edg)
{
   static int min_size = 128;
   edg = (EDGEBOUNDS *)malloc(sizeof(EDGEBOUNDS));
   edg->N = 0;
   edg->size = min_size;
   edg->bounds = (BOUND *)malloc(min_size * sizeof(BOUND));
}

/*
 *  FUNCTION: edgebounds_Destroy()
 *  SYNOPSIS: Destroy edgebounds object.
 *
 *  PURPOSE:
 *
 *  ARGS:      <edg>      Edgebounds Object
 *
 *  RETURN:
 */
void edgebounds_Destroy(EDGEBOUNDS *edg)
{
   free(edg->bounds);
   free(edg);
}


/*
 *  FUNCTION: edgebounds_Add()
 *  SYNOPSIS: Add bound to Edgebound
 *
 *  PURPOSE:
 *
 *  ARGS:      <edg>       Edgebounds,
 *             <bnd>       Bound
 *
 *  RETURN:
 */
void edgebounds_Add(EDGEBOUNDS *edg,
                    BOUND *bnd)
{
   edg->N++;

   if (edg->N == edg->size) {
      edg->bounds = (BOUND *)realloc(edg->bounds, edg->size * 2);
   }

   edg->bounds[edg->N] = *bnd;
}


/*
 *  FUNCTION: edgebounds_Resize()
 *  SYNOPSIS: Resize number of bounds in edgebound object.
 *
 *  PURPOSE:
 *
 *  ARGS:      <edg>      Edgebounds Object
 *
 *  RETURN:
 */
void edgebounds_Resize(EDGEBOUNDS *edg)
{
   int size = edg->size * 2;
   edg->bounds = (BOUND *)realloc(edg->bounds, size * sizeof(BOUND));
}


/*
 *  FUNCTION: edgebounds_Print()
 *  SYNOPSIS: Print edgebound object.
 *
 *  PURPOSE:
 *
 *  ARGS:      <bnd>      Bounds Object
 *
 *  RETURN:
 */
void edgebounds_Print(EDGEBOUNDS *edg)
{
   printf("printing edgebounds...\n");
   for (unsigned int i = 0; i < edg->N; ++i)
   {
      printf("[%d] d: %d, lb: %d, rb: %d\n", i, edg->bounds[i].diag, edg->bounds[i].lb, edg->bounds[i].rb);
   }
}


/*
 *  FUNCTION: edgebounds_Merge()
 *  SYNOPSIS: Combine edgebounds from forward and backward searches.
 *
 *  PURPOSE:
 *
 *  ARGS:      <edg_fwd>      Forward Edgebounds,
 *             <edg_bck>      Backward Edgebounds,
 *             <edg_res>      Result Edgebounds
 *
 *  RETURN:
 */
void edgebounds_Merge(EDGEBOUNDS *edg_fwd,
                      EDGEBOUNDS *edg_bck,
                      EDGEBOUNDS *edg_new)
{
   int i, j, k, x;
   bool not_merged;          /* checks whether to merge or add bounds */
   int diag_cur;             /* currently examined diagonal */
   static int tol = 0;       /* if clouds are within tolerance range, clouds are merged */

   /* allocate new edgebounds */
   int min_size = 128;
   edg_new = (EDGEBOUNDS *)malloc( sizeof(EDGEBOUNDS) );
   edg_new->size = min_size;
   edg_new->N = 0;
   edg_new->bounds = (BOUND *)malloc( min_size * sizeof(BOUND) );
   // edgebounds_Create(edg_new);

   /* temp edgebound for merging current diagonal */
   EDGEBOUNDS *edg_tmp;
   edg_tmp = (EDGEBOUNDS *)malloc( sizeof(EDGEBOUNDS) );
   edg_tmp->size = min_size;
   edg_tmp->N = 0;
   edg_tmp->bounds = (BOUND *)malloc( min_size * sizeof(BOUND) );
   // edgebounds_Create(edg_tmp);

   printf("FORWARD:\n");
   edgebounds_Print(edg_fwd);
   printf("BACKWARD:\n");
   edgebounds_Print(edg_bck);

   BOUND bound_cur;

   /* begin with minimum diagonal */
   diag_cur = min(edg_fwd->bounds[0].diag, edg_bck->bounds[0].diag);

   printf("TOTALS: fwd: %d, bck: %d\n", edg_fwd->N, edg_bck->N);

   /* iterate over all bounds */
   i = 0; j = 0; k = 0;
   while (i < edg_fwd->N || j < edg_bck->N)
   {
      printf("diag_cur: %d, diag_fwd: %d, diag_bck: %d\n", diag_cur, edg_fwd->bounds[i].diag, edg_bck->bounds[j].diag);

      printf("edg_tmp(k): %d, edg_fwd(i): %d, edg_bck(j): %d\n", k, i, j);

      /* merge all forward bounds from current diagonal */
      while (i < edg_fwd->N && edg_fwd->bounds[i].diag == diag_cur)
      {
         /* merge bounds if applicable */
         bound_cur = edg_fwd->bounds[i];
         not_merged = true;
         for (x = 0; x < k; ++x)
         {
            /* if bounds overlap, merge it into */
            if ( !(bound_cur.rb - edg_tmp->bounds[x].lb <= tol || bound_cur.lb - edg_tmp->bounds[x].rb >= tol ) )
            {
               printf("merge fwd (%d)...\n", i);
               not_merged = false;
               edg_tmp->bounds[x].lb = min(bound_cur.lb, edg_tmp->bounds[x].lb);
               edg_tmp->bounds[x].rb = max(bound_cur.rb, edg_tmp->bounds[x].rb);
               break;
            }
         }
         /* if bounds don't overlap with any previous, add new bounds  */
         if (not_merged)
         {
            printf("add fwd (%d)...\n", i);
            edg_tmp->bounds[k] = edg_fwd->bounds[i];
            edg_tmp->N++;
            ++k;
         }
         ++i;
      }

      printf("edg_tmp(k): %d, edg_fwd(i): %d, edg_bck(j): %d\n", k, i, j);
      /* merge all backward bounds from current diagonal */
      while (j < edg_bck->N && edg_bck->bounds[i].diag == diag_cur)
      {
         /* get next bound */
         bound_cur = edg_bck->bounds[i];
         not_merged = true;
         /* compare against every other  */
         for (x = 0; x < k; ++x)
         {
            /* if bounds overlap, merge it into */
            if ( !(bound_cur.rb - edg_tmp->bounds[x].lb <= tol || bound_cur.lb - edg_tmp->bounds[x].rb >= tol ) )
            {
               printf("merge bck (%d)...\n", j);
               not_merged = false;
               edg_tmp->bounds[x].lb = min(bound_cur.lb, edg_tmp->bounds[x].lb);
               edg_tmp->bounds[x].rb = max(bound_cur.rb, edg_tmp->bounds[x].rb);
               break;
            }
         }
         /* if bounds don't overlap with any previous, add new bounds  */
         if (not_merged)
         {
            printf("add bck (%d)...\n", j);
            edg_tmp->bounds[k] = edg_bck->bounds[j];
            edg_tmp->N++;
            ++k;
         }
         ++j;
      }

      int curr_diag = min(edg_fwd->bounds[i].diag, edg_bck->bounds[j].diag);
      exit(0);
   }
}


/*
 *  FUNCTION: edgebounds_Reorient()
 *  SYNOPSIS: Change edgebounds from by-diagonal to by-row.
 *
 *  PURPOSE:
 *
 *  ARGS:      <edg_diag>      Edgebounds (by-diag)
 *             <edg_row>       Edgebounds (by-row)
 *
 *  RETURN:
 */
void edgebounds_Reorient(EDGEBOUNDS *edg_diag,
                         EDGEBOUNDS *edg_row)
{
   int i, j, k, d;
   int N = edg_diag->N;

   /* convert edgebounds from (diag, leftbound, rightbound) to {(x1,y1),(x2,y2)} coords */
   for (i = 0; i < N; ++i)
   {

   }
}


/*
 *  FUNCTION: edgebounds_Merge_Reorient()
 *  SYNOPSIS: Merge and Reorient edgebounds from forward and backward searches.
 *
 *  PURPOSE:
 *
 *  ARGS:      <edg_fwd>      Forward Edgebounds,
 *             <edg_bck>      Backward Edgebounds,
 *             <edg_res>      Result Edgebounds,
 *             <st_MX>        State Matrix
 *
 *  RETURN:
 */
int edgebounds_Merge_Reorient(EDGEBOUNDS*edg_fwd,
                               EDGEBOUNDS*edg_bck,
                               EDGEBOUNDS*edg_new,
                               int Q, int T,
                               float st_MX[ NUM_NORMAL_STATES * (Q + 1) * (T + 1) ],
                               float sp_MX[ NUM_NORMAL_STATES * (Q + 1) ])
{
   /* initialize new edgebound */
   static int min_size = 128;
   edg_new->N = 0;
   edg_new->size = min_size;
   edg_new->bounds = (BOUND *)malloc(min_size * sizeof(BOUND));

   /* merge edgebounds */
   dp_matrix_Clear_X(Q, T, st_MX, sp_MX, 0);
   test_cloud(Q, T, st_MX, sp_MX, edg_fwd, 1);
   test_cloud(Q, T, st_MX, sp_MX, edg_bck, 1);
   // dp_matrix_Print(Q, T, st_MX, sp_MX);

   int num_cells = 0;
   int tot_cells = Q * T;
   float perc_cells;

   int i, j;
   int x, y1, y2;
   bool in_cloud = false;

   /* reorient cloud from (diag, bottom-left-offset, top-right-offset) to (row, left-col-offset, right-col-offset) */
   for (i = 0; i < Q; ++i)
   {
      for (j = 0; j < T; ++j)
      {
         if (in_cloud) {
            if (IMX(i, j) == 0) {
               in_cloud = false;
               y2 = j;
               edg_new->bounds[edg_new->N] = (BOUND) {y1, y2, x};
               edg_new->N += 1;
               // printf("Adding: {%d,%d,%d}...\n", x, y1, y2);
            }
            num_cells++;
         }
         else
         {
            if (IMX(i, j) > 0) {
               in_cloud = true;
               x = i;
               y1 = j;
            }
         }
      }

      if (in_cloud) {
         in_cloud = false;
         y2 = T;
         edg_new->bounds[edg_new->N] = (BOUND) {y1, y2, x};
         edg_new->N += 1;
         // printf("Adding: {%d,%d,%d}...\n", x, y1, y2);
      }
   }

   return num_cells;
}

