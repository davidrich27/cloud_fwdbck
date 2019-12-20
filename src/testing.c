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
#include "testing.h"

// macro functions
// NOTE: wrap all macro vars in parens!!
#define max(x,y) (((x) > (y)) ? (x) : (y))
#define min(x,y) (((x) < (y)) ? (x) : (y))

/* test to cycle through all diags */
void test_cycle(int Q, int T,
                float st_MX[ NUM_NORMAL_STATES * (Q*1) * (T+1) ],
                float sp_MX[ NUM_NORMAL_STATES * (Q+1) ])
{
   int d,i,j,k;
   int lb,rb,le,re;              /* right/left bounds and edges */
   int num_cells;                /* number of cells in current diagonal */
   int d_cnt, total_cnt;         /* count of cells in current diag, count of total cells */
   int d_st, d_end;              /* starting and ending diagonal indices */
   int dim_min, dim_max;         /* diagonal index where num cells reaches highest point and begins diminishing */ 
   int dim_T, dim_Q;             /* dimensions of submatrix being searched */
   COORDS start, end;            /* start and end point of alignment */

   dp_matrix_Clear(Q, T, st_MX, sp_MX);

   total_cnt = 0;

   // test coords
   start = (COORDS) {4, 8};
   end = (COORDS) {12, 14}; 

   // diag index at corners of dp matrix
   d_st = 0;
   d_end = Q + T;

   // diag index of different start points, creating submatrix
   d_st = start.i + start.j;
   // d_end = end.i + end.j;

   // dimension of submatrix
   dim_Q = Q - start.i;
   dim_T = T - start.j;

   // diag index where num cells reaches highest point and begins diminishing
   dim_min = min(Q + start.i, T + start.j);
   dim_max = max(Q + start.i, T + start.j);

   // set bounds using starting cell
   lb = start.i;
   rb = start.i + 1;
   num_cells = 0;

   printf("test cycle. (Q=%d,T=%d) (dim_Q=%d, dim_T=%d)...\n", Q, T, dim_Q, dim_T);

   /* iterate through diags */
   for (d = d_st; d <= d_end; d++, rb++)
   {
      printf("d=%d, dim_min=%d, dim_max=%d, num_cells=%d \n", d, dim_min, dim_max, num_cells);

      d_cnt = 0;

      /* is dp matrix diagonal growing or shrinking? */
      if (d <= dim_min)
         num_cells++;
      if (d > dim_max)
         num_cells--;

      printf("[%d,%d]\n", lb, rb);

      /* find diag cells that are inside matrix bounds */
      le = max(start.i, d - T);
      re = le + num_cells;

      printf("[%d,%d]\n", le, re);

      /* for testing, bounds = edges */
      if (lb < le)
         lb = le;
      if (rb > re)
         rb = re;

      printf("[%d,%d]\n", lb, rb);

      /* iterate through cells of diag */
      for (k = lb; k < rb; k++)
      {
         i = k;
         j = d - i;
         printf("(%d,%d)->", i,j);

         MMX(i,j) = i;
         IMX(i,j) = j;
         DMX(i,j) = total_cnt;

         d_cnt++;
         total_cnt++;
      }
      printf("\n");
   }

   dp_matrix_Print(Q, T, st_MX, sp_MX);
}  

/* test to cycle through all diags in reverse */
void rev_test_cycle(int Q, int T,
                float st_MX[ NUM_NORMAL_STATES * (Q*1) * (T+1) ],
                float sp_MX[ NUM_NORMAL_STATES * (Q+1) ])
{
   int d,i,j,k;
   int lb,rb,le,re;                       /* right/left bounds and edges */
   int num_cells;                         /* number of cells in current diagonal */
   int d_cnt, d_total_cnt, total_cnt;     /* count of diags, count of cells in current diag, count of total cells */
   int d_st, d_end;                       /* starting and ending diagonal indices */
   int dim_min, dim_max;                  /* diagonal index where num cells reaches highest point and begins diminishing */ 
   int dim_T, dim_Q;                      /* dimensions of submatrix being searched */
   COORDS start, end;                     /* start and end point of alignment */

   dp_matrix_Clear(Q, T, st_MX, sp_MX);

   total_cnt = 0;

   // test coords
   start = (COORDS) {4, 8};
   end = (COORDS) {14, 12}; 

   // diag index at corners of dp matrix
   d_st = 0;
   d_end = Q + T;

   // diag index of different start points, creating submatrix
   // d_st = start.i + start.j;
   d_end = end.i + end.j;

   // dimension of submatrix
   dim_Q = end.i;
   dim_T = end.j;

   // diag index where num cells reaches highest point and begins diminishing
   dim_min = min(end.i, end.j);
   dim_max = max(end.i, end.j);

   // set bounds using starting cell
   lb = end.i;
   rb = end.i + 1;
   num_cells = 0;

   printf("test cycle. (Q=%d,T=%d) (dim_Q=%d, dim_T=%d)...\n", Q, T, dim_Q, dim_T);

   /* iterate through diags */
   for (d = d_end; d >= d_st; d--, lb--)
   {
      printf("d=%d, dim_min=%d, dim_max=%d, num_cells=%d \n", d, dim_min, dim_max, num_cells);

      d_cnt = 0;

      /* is dp matrix diagonal growing or shrinking? */
      if (d >= dim_max)
         num_cells++;
      if (d < dim_min)
         num_cells--;

      printf("[%d,%d]\n", lb, rb);

      /* find diag cells that are inside matrix bounds */
      le = max(end.i - (d_end - d), 0);
      re = le + num_cells;

      printf("[%d,%d]\n", le, re);

      /* for testing, bounds = edges */
      if (lb < le)
         lb = le;
      if (rb > re)
         rb = re;

      printf("[%d,%d]\n", lb, rb);

      /* iterate through cells of diag */
      for (k = lb; k < rb; k++)
      {
         i = k;
         j = d - i;
         printf("(%d,%d)->", i,j);

         MMX(i,j) = i;
         IMX(i,j) = j;
         DMX(i,j) = total_cnt;

         d_cnt++;
         total_cnt++;
      }
      printf("\n");
   }

   dp_matrix_Print(Q, T, st_MX, sp_MX);
}  

/* test to show the cloud area */
void test_cloud(int Q, int T,
                float st_MX[ NUM_NORMAL_STATES * (Q*1) * (T+1) ],
                float sp_MX[ NUM_NORMAL_STATES * (Q+1) ],
                EDGEBOUNDS* edg,
                float val)
{
   int d_st, d_end, d_cnt;
   int d,i,j,k;
   int rb,lb;

   d_st = edg->bounds[0].diag;
   d_end = edg->bounds[edg->N-1].diag;
   d_cnt = 0;

   for (d = d_st; d < d_end; d++, d_cnt++)
   {
      lb = edg->bounds[d_cnt].lb;
      rb = edg->bounds[d_cnt].rb;

      for (k = lb; k < rb; k++)
      {
         i = k;
         j = d - i;

         MMX(i,j) += 1;
         IMX(i,j) = 1;
         DMX(i,j) += val;
      }
   }
}