/*******************************************************************************
 *  @file main.c
 *  @brief Parses files, 
 *
 *  @author Dave Rich (devrek)
 *  @bug Lots.
 *******************************************************************************/

/* imports */
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include <ctype.h>

/* local imports (after struct declarations) */

/* data structures and file parsers */
#include "structs.h"
#include "misc.h"
#include "hmm_parser.h"
#include "edgebounds_obj.h"

/* quadratic space algs */
#include "viterbi.h"
#include "forward_backward.h"
#include "cloud_search.h"

/* linear space algs */
// #include "viterbi3.h"
// #include "forward_backward3.h"
#include "cloud_search3.h"

/* naive algs */
#include "cloud_search_naive.h"

#include "testing.h"


void parse_args(int argc, char *argv[], ARGS *args);
void test(char *hmm_file, char *fasta_file, float alpha, int beta);
void test_barrage(char *hmm_file, char *fasta_file);

/* MAIN */
int main (int argc, char *argv[])
{
   printf("begin main...\n");
   ARGS *args = (ARGS *)malloc( sizeof(ARGS) );

   char *hmm_file, *fasta_file;
   char *arg;
   int  i, j;

   /* DEFAULTS */
   args->target_hmm_file = "data/test1_2.hmm";
   args->query_fasta_file = "data/test1_1.fa";
   args->alpha = 20.0f;
   args->beta = 5;

   printf("parsing args...\n");
   parse_args(argc, argv, args);

   printf("begin test...\n");
   test(args->target_hmm_file, args->query_fasta_file, args->alpha, args->beta);
}

/* Parses Arguments from the command line */
void parse_args (int argc, char *argv[], ARGS *args)
{
   int i, len;
   int num_main_args, max_main_args;
   
   num_main_args = 0;
   max_main_args = 2;

   for (i = 1; i < argc; ++i)
   {
      if ( argv[i][0] == '-' )
      {
         switch (argv[i][1]) {
            case 'a':   /* alpha */
               i++;
               if (i < argc) {
                  args->alpha = atof(argv[i]);
                  break;
               } else {
                  fprintf(stderr, "Error: -a flag requires argument.\n");
               }
               break;
            case 'b':   /* beta */
               i++;
               if (i < argc) {
                  args->beta = atoi(argv[i]);
               } else {
                  fprintf(stderr, "Error: -b flag requires argument.\n");
               }
               break;
            case 'o':   /* outfile */
               i++;
               if (i < argc) {
                  args->outfile = argv[i];
               } else {
                  fprintf(stderr, "Error: -o flag requires argument.\n");
               } 
               break;
            default:
               fprintf(stderr, "Error: -%c is not a valid flag.\n", argv[i][1]);
               exit(1);
               break;
         }
      }
      else
      {
         len = get_str_len(argv[i]);
         if (num_main_args == 0)
         {
            args->target_hmm_file = argv[i];
         }
         else if (num_main_args == 1)
         {
            args->query_fasta_file = argv[i];
         }
         else
         {
            fprintf(stderr, "Error: Too many main arguments.\n");
         }
         num_main_args++;
      }
   }
}

/* unit test */
void test(char *hmm_file, char *fasta_file, float alpha, int beta)
{
   /* PRINT ARGS */
   printf("HMM_FILENAME: %s\n", hmm_file);
   printf("FASTA_FILENAME: %s\n", fasta_file);
   printf("ALPHA: %f\n", alpha);
   printf("BETA: %d\n\n", beta);
   SCORES *scores = (SCORES*)malloc( sizeof(SCORES) );

   printf("building hmm profile...\n");
   /* get target profile */
   HMM_PROFILE *target_prof = (HMM_PROFILE *)malloc( sizeof(HMM_PROFILE) );
   hmmprofile_Create(target_prof, hmm_file);
   // hmmprofile_Display(target_prof);
   hmmprofile_Save(target_prof, "output/myversion.pre-profile.tsv");

   printf("configuring...\n");
   /* mode choices: MODE_NONE, MODE_MULTILOCAL, MODE_MULTIGLOCAL, MODE_UNILOCAL, MODE_UNIGLOCAL */
   char* modes[] = { "None", "Multi-local", "Multi-glocal", "Uni-local", "Uni-glocal" }; 
   int mode; 
   // mode = MODE_MULTILOCAL;    /* HMMER standard mode (allows jumps) */
   mode = MODE_UNILOCAL;   /* Cloud mode (prohibiits jumps) */
   printf("MODE: %s\n", modes[mode]);
   hmmprofile_Config(target_prof, mode);
   // hmmprofile_Display(target_prof);
   hmmprofile_Save(target_prof, "output/myversion.post-profile.tsv");
   int T = target_prof->leng;

   printf("building query sequence...\n");
   /* get query sequence */
   SEQ *query_seq = (SEQ *)malloc( sizeof(SEQ) );
   seq_Create(query_seq, fasta_file);
   seq_Display(query_seq);
   int Q = query_seq->leng;

   printf("I/O was successful!\n\n");

   /* allocate memory to store results */
   float sc, perc_cells;
   int num_cells, window_cells; 
   int tot_cells = (Q + 1) * (T + 1);
   RESULTS *res = (RESULTS *)malloc( sizeof(RESULTS) );
   TRACEBACK *tr = (TRACEBACK *)malloc( sizeof(TRACEBACK) );
   EDGEBOUNDS *edg_fwd = (EDGEBOUNDS *)malloc( sizeof(EDGEBOUNDS) );
   EDGEBOUNDS *edg_bck = (EDGEBOUNDS *)malloc( sizeof(EDGEBOUNDS) );
   EDGEBOUNDS *edg = (EDGEBOUNDS *)malloc( sizeof(EDGEBOUNDS) );

   /* allocate memory for quadratic algs (for testing) */
   float *st_MX = (float *) malloc( sizeof(float) * (NUM_NORMAL_STATES * (Q+1) * (T+1)) );
   float *sp_MX = (float *) malloc( sizeof(float) * (NUM_SPECIAL_STATES * (Q+1)) );
   /* allocate memory for cloud matrices (for testing) */
   float *st_MX_cloud = (float *) malloc( sizeof(float) * (NUM_NORMAL_STATES * (Q+1) * (T+1)) );
   float *sp_MX_cloud = (float *) malloc( sizeof(float) * (NUM_SPECIAL_STATES * (Q+1)) );
   /* allocate memory for linear algs */
   float *st_MX3 = (float *) malloc( sizeof(float) * (NUM_NORMAL_STATES * (Q+1) * 3) );

   /* run viterbi algorithm */
   printf("=== VITERBI -> START ===\n");
   sc = viterbi_Run(query_seq, target_prof, Q, T, st_MX, sp_MX, res, tr);
   printf("Viterbi Score: %f\n", sc);
   scores->viterbi_sc = sc;
   // dp_matrix_Print(Q, T, st_MX, sp_MX);
   dp_matrix_Save(Q, T, st_MX, sp_MX, "output/myversion.viterbi.tsv");
   printf("=== VITERBI -> END ===\n\n");

   /* run traceback of viterbi */
   printf("=== TRACEBACK -> START ===\n");
   traceback_Build(query_seq, target_prof, Q, T, st_MX, sp_MX, tr);
   // traceback_Print(tr);
   traceback_Save(tr, "output/myversion.traceback_list.tsv");
   dp_matrix_Clear(Q, T, st_MX, sp_MX);
   traceback_Show(Q, T, st_MX, sp_MX, tr);
   printf("START: (%d,%d) -> END: (%d,%d)\n", tr->first_m.i, tr->first_m.j, tr->last_m.i, tr->last_m.j);
   window_cells = (tr->last_m.i - tr->first_m.i) * (tr->last_m.j - tr->first_m.j);
   dp_matrix_trace_Save(Q, T, st_MX, sp_MX, tr, "output/myversion.traceback.tsv");
   printf("=== TRACEBACK -> END ===\n\n");

   /* run forward/backward algorithms */
   init_Logsum();

   printf("=== FORWARD -> START ===\n");
   dp_matrix_Clear(Q, T, st_MX, sp_MX);
   sc = forward_Run(query_seq, target_prof, Q, T, st_MX, sp_MX, res);
   printf("Forward Score: %f\n", sc);
   scores->fwd_sc = sc;
   // dp_matrix_Print(Q, T, st_MX, sp_MX);
   dp_matrix_trace_Save(Q, T, st_MX, sp_MX, tr, "output/myversion.forward.tsv");
   printf("=== FORWARD -> END ===\n\n");

   printf("=== BACKWARD -> START ===\n");
   sc = backward_Run(query_seq, target_prof, Q, T, st_MX, sp_MX, res);
   printf("Backward Score: %f\n", sc);
   scores->bck_sc = sc;
   // dp_matrix_Print(Q, T, st_MX, sp_MX);
   dp_matrix_trace_Save(Q, T, st_MX, sp_MX, tr, "output/myversion.backward.tsv");
   printf("=== BACKWARD -> END ===\n\n");

   /* TEST */
   printf("=== TEST -> START ===\n");
   dp_matrix_Clear_X(Q, T, st_MX, sp_MX, 0);
   fwd_test_cycle(Q, T, st_MX, sp_MX, tr);
   bck_test_cycle(Q, T, st_MX, sp_MX, tr);
   dp_matrix_trace_Save(Q, T, st_MX, sp_MX, tr, "output/myversion.test_fwd.tsv");
   printf("=== TEST -> END ===\n\n");

   /* cloud forward */
   printf("=== CLOUD FORWARD -> START ===\n");
   cloud_forward_Run(query_seq, target_prof, Q, T, st_MX, sp_MX, tr, edg_fwd, alpha, beta);
   dp_matrix_trace_Save(Q, T, st_MX, sp_MX, tr, "output/myversion.cloud_fwd_vals.tsv");
   dp_matrix_Clear_X(Q, T, st_MX, sp_MX, 0);
   cloud_Fill(Q, T, st_MX, sp_MX, edg_fwd, 1, MODE_DIAG);
   dp_matrix_trace_Save(Q, T, st_MX, sp_MX, tr, "output/myversion.cloud_fwd_diag.tsv");
   edgebounds_Save(edg_fwd, "output/myversion.edgebounds_fwd_diag.tsv");
   printf("=== CLOUD FORWARD -> END ===\n\n");

   /* cloud backward */
   printf("=== CLOUD BACKWARD -> START ===\n");
   cloud_backward_Run(query_seq, target_prof, Q, T, st_MX, sp_MX, tr, edg_bck, alpha, beta);
   dp_matrix_trace_Save(Q, T, st_MX, sp_MX, tr, "output/myversion.cloud_bck_vals.tsv");
   dp_matrix_Clear_X(Q, T, st_MX, sp_MX, 0);
   cloud_Fill(Q, T, st_MX, sp_MX, edg_bck, 1, MODE_DIAG);
   dp_matrix_trace_Save(Q, T, st_MX, sp_MX, tr, "output/myversion.cloud_bck_diag.tsv");
   edgebounds_Save(edg_bck, "output/myversion.edgebounds_bck_diag.tsv");
   printf("=== CLOUD BACKWARD -> END ===\n\n");

   /* FORWARD/BACKWARD MODE OPTIONS = {MODE_LINEAR, MODE_QUAD, MODE_NAIVE} */
   int fwdbck_mode = MODE_NAIVE;

   if (fwdbck_mode == MODE_LINEAR) 
   {
      // printf("Using linear space algs...\n");

      // /* cloud forward */
      // printf("=== CLOUD FORWARD -> START ===\n");
      // cloud_forward_Run3(query_seq, target_prof, Q, T, st_MX3, sp_MX, res, tr, edg_fwd, alpha, beta);
      // // edgebounds_Print(edg_fwd);
      // // dp_matrix_Print(Q, T, st_MX, sp_MX);
      // // dp_matrix_Save(Q, T, st_MX, sp_MX, "output/myversion.cloud_fwd.3.tsv");
      // printf("=== CLOUD FORWARD -> END ===\n\n");

      // /* cloud backward */
      // printf("=== CLOUD BACKWARD -> START ===\n");
      // cloud_backward_Run3(query_seq, target_prof, Q, T, st_MX3, sp_MX, res, tr, edg_bck, alpha, beta);
      // // edgebounds_Print(edg_bck);
      // // dp_matrix_Print(Q, T, st_MX, sp_MX);
      // // dp_matrix_Save(Q, T, st_MX, sp_MX, "output/myversion.cloud_bck.3.tsv");
      // printf("=== CLOUD BACKWARD -> END ===\n\n");

      // // printf("=== MERGE CLOUD -> START ===\n");
      // // edgebounds_Merge(edg_bck, edg_fwd, edg);
      // // printf("=== MERGE CLOUD -> END ===\n");

      // // printf("=== REORIENT CLOUD -> START ===\n");
      // // edgebounds_Reorient(edg, edg_row);
      // // printf("=== REORIENT CLOUD -> END ===\n");

      // printf("=== REORIENT/MERGE CLOUD -> START ===\n");
      // num_cells = edgebounds_Merge_Reorient_Cloud(edg_fwd, edg_bck, edg, Q, T, st_MX, sp_MX);
      // perc_cells = (float)num_cells/(float)tot_cells;
      // printf("Cells Computed = %d/%d = %.5f\n", num_cells, tot_cells, perc_cells);
      // // edgebounds_Print(edg);
      // printf("=== REORIENT/MERGE CLOUD -> END ===\n\n");

      // /* Data viz to see cloud */
      // // printf("=== TEST CLOUD -> START ===\n");
      // // dp_matrix_Clear_X(Q, T, st_MX, sp_MX, 0);
      // // cloud_Fill(Q, T, st_MX, sp_MX, edg_fwd, 1, MODE_DIAG);
      // // cloud_Fill(Q, T, st_MX, sp_MX, edg_bck, -1, MODE_DIAG);
      // // // dp_matrix_Print(Q, T, st_MX, sp_MX);
      // // dp_matrix_Save(Q, T, st_MX, sp_MX, "output/myversion.cloud.3.tsv");
      // // printf("=== TEST CLOUD -> END ===\n");

      // /* bounded forward */
      // printf("=== BOUNDED FORWARD -> START ===\n");
      // dp_matrix_Clear_X(Q, T, st_MX, sp_MX, -INF);
      // sc = forward_bounded_Run3(query_seq, target_prof, Q, T, st_MX3, sp_MX, res, edg);
      // printf("Bounded Forward Score: %.5f\n", sc);
      // // dp_matrix_Print(Q, T, st_MX, sp_MX);
      // // dp_matrix_Save(Q, T, st_MX, sp_MX, "output/myversion.bounded_fwd.3.tsv");
      // printf("=== BOUNDED FORWARD -> END ===\n\n");

      // /* bounded backward */
      // printf("=== BOUNDED BACKWARD -> START ===\n");
      // dp_matrix_Clear_X(Q, T, st_MX, sp_MX, -INF);
      // sc = backward_bounded_Run3(query_seq, target_prof, Q, T, st_MX3, sp_MX, res, edg);
      // printf("Bounded Backward Score: %.5f\n", sc);
      // // dp_matrix_Print(Q, T, st_MX, sp_MX);
      // // dp_matrix_Save(Q, T, st_MX, sp_MX, "output/myversion.bounded_bck.3.tsv");
      // printf("=== BOUNDED BACKWARD -> END ===\n\n");
   }
   else
   if (fwdbck_mode == MODE_QUAD) 
   {
      // printf("Using quadratic space algs...\n");

      // /* cloud forward */
      // printf("=== CLOUD FORWARD -> START ===\n");
      // cloud_forward_Run(query_seq, target_prof, Q, T, st_MX, sp_MX, res, tr, edg_fwd, alpha, beta);
      // dp_matrix_Save(Q, T, st_MX, sp_MX, "output/myversion.cloud_fwd_vals.tsv");
      // // edgebounds_Print(edg_fwd);
      // // dp_matrix_Print(Q, T, st_MX, sp_MX);
      // dp_matrix_Clear_X(Q, T, st_MX, sp_MX, 0);
      // cloud_Fill(Q, T, st_MX, sp_MX, edg_fwd, 1, MODE_DIAG);
      // dp_matrix_Copy(Q, T, st_MX, sp_MX, st_MX_fwd_d, sp_MX_fwd_d);
      // dp_matrix_Save(Q, T, st_MX, sp_MX, "output/myversion.cloud_fwd_diag.tsv");
      // edgebounds_Save(edg_fwd, "output/myversion.edgebounds_fwd_diag.tsv");
      // printf("=== CLOUD FORWARD -> END ===\n\n");

      // /* cloud backward */
      // printf("=== CLOUD BACKWARD -> START ===\n");
      // cloud_backward_Run(query_seq, target_prof, Q, T, st_MX, sp_MX, res, tr, edg_bck, alpha, beta);
      // dp_matrix_Save(Q, T, st_MX, sp_MX, "output/myversion.cloud_bck_vals.tsv");
      // // edgebounds_Print(edg_bck);
      // // dp_matrix_Print(Q, T, st_MX, sp_MX);
      // dp_matrix_Clear_X(Q, T, st_MX, sp_MX, 0);
      // cloud_Fill(Q, T, st_MX, sp_MX, edg_bck, 1, MODE_DIAG);
      // dp_matrix_Copy(Q, T, st_MX, sp_MX, st_MX_bck_d, sp_MX_bck_d);
      // dp_matrix_Save(Q, T, st_MX, sp_MX, "output/myversion.cloud_bck_diag.tsv");
      // edgebounds_Save(edg_bck, "output/myversion.edgebounds_bck_diag.tsv");
      // printf("=== CLOUD BACKWARD -> END ===\n\n");

      // printf("=== MERGE & REORIENT CLOUD -> START ===\n");
      // /* merge and reorient edgebounds from by-diag to by-row */
      // edgebounds_Merge_Reorient_Cloud(edg_fwd, edg_bck, edg, Q, T, st_MX, sp_MX);
      // dp_matrix_Clear_X(Q, T, st_MX, sp_MX, 0);
      // cloud_Fill(Q, T, st_MX, sp_MX, edg, 1, MODE_ROW);
      // num_cells = cloud_Cell_Count(Q, T, st_MX, sp_MX);
      // perc_cells = (float)num_cells/(float)tot_cells;
      // printf("Total Cells Computed = %d/%d = %.4f\n", num_cells, tot_cells, perc_cells);
      // // edgebounds_Print(edg);
      // edgebounds_Save(edg, "output/myversion.edgebounds_row.tsv");
      // dp_matrix_Save(Q, T, st_MX, sp_MX, "output/myversion.cloud.tsv");
      // printf("=== MERGE & REORIENT CLOUD -> END ===\n\n");

      // /* bounded forward */
      // printf("=== BOUNDED FORWARD -> START ===\n");
      // dp_matrix_Clear(Q, T, st_MX, sp_MX);
      // sc = forward_bounded_Run(query_seq, target_prof, Q, T, st_MX, sp_MX, res, edg);
      // printf("Bounded Forward Score: %.5f\n", sc);
      // // dp_matrix_Print(Q, T, st_MX, sp_MX);
      // dp_matrix_Save(Q, T, st_MX, sp_MX, "output/myversion.bounded_fwd.tsv");
      // printf("=== BOUNDED FORWARD -> END ===\n\n");

      // /* bounded backward */
      // printf("=== BOUNDED BACKWARD -> START ===\n");
      // dp_matrix_Clear(Q, T, st_MX, sp_MX);
      // sc = backward_bounded_Run(query_seq, target_prof, Q, T, st_MX, sp_MX, res, edg);
      // printf("Bounded Backward Score: %.5f\n", sc);
      // // dp_matrix_Print(Q, T, st_MX, sp_MX);
      // dp_matrix_Save(Q, T, st_MX, sp_MX, "output/myversion.bounded_bck.tsv");
      // printf("=== BOUNDED BACKWARD -> END ===\n\n");
   }
   else 
   if (fwdbck_mode == MODE_NAIVE) 
   {
      printf("using naive algs...\n");

      /* merge forward and backward clouds, then reorient edgebounds from by-diag to by-row */
      printf("=== MERGE & REORIENT CLOUD -> START ===\n");
      edgebounds_Merge_Reorient_Cloud(edg_fwd, edg_bck, edg, Q, T, st_MX, sp_MX);
      dp_matrix_trace_Save(Q, T, st_MX, sp_MX, tr, "output/myversion.cloud.diags.tsv");
      dp_matrix_Clear_X(Q, T, st_MX, sp_MX, 0);
      cloud_Fill(Q, T, st_MX, sp_MX, edg, 1, MODE_ROW);
      edgebounds_Save(edg, "output/myversion.edgebounds_row.tsv");
      dp_matrix_trace_Save(Q, T, st_MX, sp_MX, tr, "output/myversion.cloud.rows.tsv");
      num_cells = cloud_Cell_Count(Q, T, st_MX, sp_MX);
      scores->perc_cells = (float)num_cells/(float)tot_cells;
      printf("Perc. Total Cells Computed = %d/%d = %f\n", num_cells, tot_cells, scores->perc_cells);
      scores->perc_window = (float)num_cells/(float)window_cells;
      printf("Perc. Window Cells Computed = %d/%d = %f\n", num_cells, window_cells, scores->perc_window);
      dp_matrix_Copy(Q, T, st_MX, sp_MX, st_MX_cloud, sp_MX_cloud);
      printf("=== MERGE & REORIENT CLOUD -> END ===\n\n");

      // /* create cloud that covers entire matrix (full fwd/bck) */
      // printf("=== TEST CLOUD -> START ===\n");
      // dp_matrix_Clear_X(Q, T, st_MX_cloud, sp_MX_cloud, 1);
      // edgebounds_Build_From_Cloud(edg, Q, T, st_MX_cloud, MODE_ROW);
      // edgebounds_Save(edg, "output/myversion.edgebounds_full.tsv");
      // dp_matrix_trace_Save(Q, T, st_MX_cloud, sp_MX_cloud, tr, "output/myversion.cloud.full_fwdbck.tsv");
      // printf("=== TEST CLOUD -> END ===\n\n");

      /* bounded forward */
      printf("=== BOUNDED FORWARD -> START ===\n");
      dp_matrix_Clear(Q, T, st_MX, sp_MX);
      forward_Bounded_Naive_Run(query_seq, target_prof, Q, T, st_MX, sp_MX, st_MX_cloud, &sc);
      printf("Bounded Forward Score (Naive): %f\n", sc);
      dp_matrix_trace_Save(Q, T, st_MX, sp_MX, tr, "output/myversion.bounded_fwd.naive.tsv");
      
      dp_matrix_Clear(Q, T, st_MX, sp_MX);
      forward_bounded_Run(query_seq, target_prof, Q, T, st_MX, sp_MX, edg, &sc);
      printf("Bounded Forward Score (Quadratic): %f\n", sc);
      scores->cloud_fwd_sc = sc;
      dp_matrix_trace_Save(Q, T, st_MX, sp_MX, tr, "output/myversion.bounded_fwd.quadratic.tsv");

      // backward_bounded_Run3(query_seq, target_prof, Q, T, st_MX3, sp_MX, edg, &sc);
      // printf("Bound Backward Score (Linear): %f\n", sc);
      printf("=== BOUNDED FORWARD -> END ===\n\n");

      /* bounded backward */
      printf("=== BOUNDED BACKWARD -> START ===\n");
      dp_matrix_Clear(Q, T, st_MX, sp_MX);
      backward_Bounded_Naive_Run(query_seq, target_prof, Q, T, st_MX, sp_MX, st_MX_cloud, &sc);
      printf("Bounded Backward Score (Naive): %f\n", sc);
      scores->cloud_bck_sc = sc;
      dp_matrix_trace_Save(Q, T, st_MX, sp_MX, tr, "output/myversion.bounded_bck.naive.tsv");

      dp_matrix_Clear(Q, T, st_MX, sp_MX);
      backward_bounded_Run(query_seq, target_prof, Q, T, st_MX, sp_MX, edg, &sc);
      printf("Bounded Backward Score (Quadratic): %f\n", sc);
      dp_matrix_trace_Save(Q, T, st_MX, sp_MX, tr, "output/myversion.bounded_bck.quadratic.tsv");

      // backward_bounded_Run3(query_seq, target_prof, Q, T, st_MX3, sp_MX, edg, &sc);
      // printf("Bound Backward Score (Linear): %f\n", sc);
      printf("=== BOUNDED BACKWARD -> END ===\n\n");

      /* sample stats */
      FILE *fp = fopen("scripts/cloud_stats.tsv", "a+");
      fprintf(fp, "%s\t", hmm_file);
      fprintf(fp, "%f\t", scores->viterbi_sc);
      fprintf(fp, "%f\t", scores->fwd_sc);
      fprintf(fp, "%f\t", scores->bck_sc);
      fprintf(fp, "%f\t", scores->cloud_fwd_sc);
      fprintf(fp, "%f\t", scores->cloud_bck_sc);
      fprintf(fp, "%f\t", alpha);
      fprintf(fp, "%d\t", beta);
      fprintf(fp, "%f\t", scores->perc_cells);
      fprintf(fp, "%f\n", scores->perc_window);
      fclose(fp);
   }

   /* TODO: free memory! */
   free(res);
   free(tr);
   free(edg_fwd);
   free(edg_bck);
   free(edg);
   free(st_MX);
   free(sp_MX);
   free(st_MX_cloud);
   free(sp_MX_cloud);
   free(st_MX3);

   printf("...test finished. \n");
}


/* standard pipeline */
void pipeline() 
{

}


/* tests a query/target against several possible  */
void test_barrage(char *hmm_file, char *fasta_file)
{
   int alpha = 5;
}