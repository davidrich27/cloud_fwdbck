/*******************************************************************************
 *  @file misc.h
 *  @brief Miscellaneous helper functions.
 *
 *  @author Dave Rich (devrek)
 *  @bug Lots.
 *******************************************************************************/

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

/* table of logsum values */
static float logsum_lookup[LOGSUM_TBL];
bool logsum_initialized = false;

/* Get max value of two floats */
float calc_Max (float x, float y)
{
   if (x > y)
   {
      return x;
   }
   return y;
}

float calc_Min (float x, float y)
{
   if (x < y)
   {
      return x;
   }
   return y;
}

/* initialize the logsum table */
void init_Logsum ()
{
   if (logsum_initialized == false)
   {
      logsum_initialized = true;
      int i;
      for (i = 0; i < LOGSUM_TBL; i++)
      {
         logsum_lookup[i] = log(1. + exp((double) -i / LOGSUM_SCALE));
      }
   }
}

/* takes in two log numbers and returns the log of their real sum (approx) */
/* log( exp(x) + exp(y) ) */
float calc_Logsum (float x, float y)
{
   float max, min;

   if (x > y)
   {
      max = x;
      min = y;
   }
   else
   {
      max = y;
      min = x;
   }

   return (min == log(0) || (max - min) >= 15.7f) ? 
      max : max + logsum_lookup[(int)((max - min) * LOGSUM_SCALE)];
}

/* takes in two log numbers and returns the log of their real sum (exact) */
float calc_Logsum_exact(float x, float y)
{
   return log( exp(x) + exp(y) );
}

void print_Logsum()
{
   printf("=== LOGSUM TABLE ===\n");
   for (int i = 0; i < 16000; i+=160)
   {
      printf("[%d] %.2f\n", i, logsum_lookup[i]);
   }
   printf("\n\n");
}

/* 
 *  FUNCTION:  dp_matrix_Print()
 *  SYNOPSIS:  Print out dynamic programming matrix to screen.
 *
 *  PURPOSE:
 *
 *  ARGS:      <Q>         query length, 
 *             <T>         target length,
 *             <st_MX>     Normal State (Match, Insert, Delete) Matrix,
 *             <sp_MX>     Special State (J,N,B,C,E) Matrix
 *
 *  RETURN: 
 */
void dp_matrix_Print (const int Q, const int T, 
                      const float st_MX[ NUM_NORMAL_STATES * (Q+1) * (T+1) ], 
                      const float sp_MX[ NUM_SPECIAL_STATES * (Q+1) ])
{
   /* PRINT resulting dp matrix */
   printf("\n\n==== DP MATRIX BEGIN ==== \n");
   /* Header */
   printf("\t");
   for (int i = 0; i <= T; i++)
   {
      printf("%d\t", i);
   }
   printf("\n");

   /* Row-by-Row */
   for (int i = 0; i < Q+1; i++)
   {
      printf( "%d M\t", i );
      for (int j = 0; j <= T; j++)
      {
         printf( "%.3f\t", MMX(i,j) );
      }
      printf("\n");

      printf( "%d I\t", i );
      for (int j = 0; j <= T; j++)
      {
         printf( "%.3f\t", IMX(i,j) );
      }
      printf("\n");

      printf( "%d D\t", i );
      for (int j = 0; j <= T; j++)
      {
         printf( "%.3f\t", DMX(i,j) );
      }
      printf("\n\n");
   }

     printf("=== SPECIAL STATES ===\n");
      printf("N:\t");
      for (int i = 0; i <= Q; i++)
         { printf( "%.3f\t", XMX(SP_N,i) ); }
      printf("\n");
      printf("J:\t");
      for (int i = 0; i <= Q; i++)
         { printf( "%.3f\t", XMX(SP_J,i) ); }
      printf("\n");
      printf("E:\t");
      for (int i = 0; i <= Q; i++)
         { printf( "%.3f\t", XMX(SP_E,i) ); }
      printf("\n");
      printf("C:\t");
      for (int i = 0; i <= Q; i++)
         { printf( "%.3f\t", XMX(SP_C,i) ); }
      printf("\n");
      printf("B:\t");
      for (int i = 0; i <= Q; i++)
         { printf( "%.3f\t", XMX(SP_B,i) ); }
      printf("\n");

   printf("==== DP MATRIX END ==== \n\n");
}

/* 
 *  FUNCTION:  dp_matrix_Clear()
 *  SYNOPSIS:  Clear all matrix values to -INF. (for testing)
 *
 *  PURPOSE:
 *
 *  ARGS:      <Q>         query length, 
 *             <T>         target length,
 *             <st_MX>     Normal State (Match, Insert, Delete) Matrix,
 *             <sp_MX>     Special State (J,N,B,C,E) Matrix
 *
 *  RETURN: 
 */
void dp_matrix_Clear (const int Q, const int T, 
                      float st_MX[ NUM_NORMAL_STATES * (Q+1) * (T+1) ], 
                      float sp_MX[ NUM_SPECIAL_STATES * (Q+1) ])
{
   for (int i = 0; i < Q+1; i++)
   {
      for (int j = 0; j <= T; j++)
         MMX(i,j) = IMX(i,j) = DMX(i,j) = -INF;

      for (int j = 0; j < NUM_SPECIAL_STATES; j++)
         XMX(j,i) = -INF;
   }
}


/* 
 *  FUNCTION:  dp_matrix_Save()
 *  SYNOPSIS:  Clear all matrix values to -INF. (for testing)
 *
 *  PURPOSE:
 *
 *  ARGS:      <Q>         query length, 
 *             <T>         target length,
 *             <st_MX>     Normal State (Match, Insert, Delete) Matrix,
 *             <sp_MX>     Special State (J,N,B,C,E) Matrix
 *             <f>         Filename
 *
 *  RETURN: 
 */
void dp_matrix_Save (const int Q, const int T, 
                     const float st_MX[ NUM_NORMAL_STATES * (Q+1) * (T+1) ], 
                     const float sp_MX[ NUM_SPECIAL_STATES * (Q+1) ],
                     const char *_filename_)
{
   FILE *fp;
   fp = fopen(_filename_, "w");

   /* PRINT resulting dp matrix */
   fprintf(fp, "==== DP MATRIX ==== \n");
   fprintf(fp, "DIM: %d,%d\n\n", Q, T);
   /* Header */
   fprintf(fp, "\t");
   for (int i = 0; i <= T; i++)
   {
      fprintf(fp, "%d\t", i);
   }
   fprintf(fp, "\n");

   /* Row-by-Row */
   for (int i = 0; i < Q+1; i++)
   {
      fprintf(fp, "%d M\t", i );
      for (int j = 0; j <= T; j++)
      {
         fprintf(fp, "%.3f\t", MMX(i,j) );
      }
      fprintf(fp, "\n");

      fprintf(fp, "%d I\t", i );
      for (int j = 0; j <= T; j++)
      {
         fprintf(fp, "%.3f\t", IMX(i,j) );
      }
      fprintf(fp, "\n");

      fprintf(fp, "%d D\t", i );
      for (int j = 0; j <= T; j++)
      {
         fprintf(fp, "%.3f\t", DMX(i,j) );
      }
      fprintf(fp, "\n\n");
   }

      fprintf(fp, "=== SPECIAL STATES ===\n");
      fprintf(fp, "N:\t");
      for (int i = 0; i <= Q; i++)
         { fprintf(fp, "%.3f\t", XMX(SP_N,i) ); }
      fprintf(fp, "\n");
      fprintf(fp, "J:\t");
      for (int i = 0; i <= Q; i++)
         { fprintf(fp, "%.3f\t", XMX(SP_J,i) ); }
      fprintf(fp, "\n");
      fprintf(fp, "E:\t");
      for (int i = 0; i <= Q; i++)
         { fprintf(fp, "%.3f\t", XMX(SP_E,i) ); }
      fprintf(fp, "\n");
      fprintf(fp, "C:\t");
      for (int i = 0; i <= Q; i++)
         { fprintf(fp, "%.3f\t", XMX(SP_C,i) ); }
      fprintf(fp, "\n");
      fprintf(fp, "B:\t");
      for (int i = 0; i <= Q; i++)
         { fprintf(fp, "%.3f\t", XMX(SP_B,i) ); }
      fprintf(fp, "\n");

   fclose(fp);
}