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

#include "testing.h"


void parse_args(int argc, char *argv, ARGS *args);
void test(char *hmm_file, char *fasta_file);

/* MAIN */
int main (int argc, char *argv[])
{
   ARGS *args = (ARGS *)malloc( sizeof(ARGS) );

   char *hmm_file, *fasta_file;

   if (argc == 3)
   {
      hmm_file = argv[1];
      fasta_file = argv[2];
   }
   else
   {
      /* TEST 1 */
      char *hmm_file = "../data/test1_2.hmm";
      char *fasta_file = "../data/test1_1.fa";

      printf("Usage: <target_hmm> <query_fa>\n");
      exit(0);
   }

   test(hmm_file, fasta_file);
}

/* Parses Arguments from the command line */
void parse_args (int argc, char *argv, ARGS *args)
{
   /* flagged args */
   int opt;
   while ((opt = getopt(argc, &argv, ":if:hoqt")) != -1)
   {
      switch(opt)
      {
         case 'h':
            printf("This is the help dialogue.\n");
            exit(EXIT_SUCCESS);
         case 'o':
            args->outfile = optarg;
            printf("Sending output to file: %s\n", optarg);
            break;
      }
   }

   /* if improper number of args */
   int numArgs = optind - argc;
   if (numArgs != 6)
   {
      printf("Usage: <query_file> <start end> <target_hmm> <start end> \n");
      exit(0);
   }

   /* non-flagged args */
   args->infile_query = &argv[optind];
   optind++;
   args->range_query->start = atoi(&argv[optind]);
   optind++;
   args->range_query->end = atoi(&argv[optind]);
   optind++;
   args->infile_target = &argv[optind];
   optind++;
   args->range_target->start = atoi(&argv[optind]);
   optind++;
   args->range_target->end = atoi(&argv[optind]);
   return;
}

/* unit test */
void test(char *hmm_file, char *fasta_file)
{
   //parse_args(argc, *argv, args);

   /* load substitution matrix */
   // char *submat_file = "../data/submat/blosum62.submat";
   // SUBMAT *submat = (SUBMAT *)malloc( sizeof(SUBMAT) );
   // submat_Create(submat, submat_file);
   // submat_Display(submat);

   #ifdef DEBUG
   printf("This is a debugging...\n");
   #endif

   printf("building hmm profile...\n");
   /* get target profile */
   HMM_PROFILE *target_prof = (HMM_PROFILE *)malloc( sizeof(HMM_PROFILE) );
   hmmprofile_Create(target_prof, hmm_file);
   // hmmprofile_Display(target_prof);
   hmmprofile_Save(target_prof, "output/myversion.pre-profile.tsv");

   printf("configuring...\n");
   /* mode choices: MODE_NONE, MODE_MULTILOCAL, MODE_MULTIGLOCAL, MODE_UNILOCAL, MODE_UNIGLOCAL */
   int mode;
   // mode = MODE_MULTILOCAL;
   mode = MODE_UNILOCAL;
   hmmprofile_Config(target_prof, mode);
   // hmmprofile_Display(target_prof);
   hmmprofile_Save(target_prof, "output/myversion.post-profile.tsv");
   int T = target_prof->leng;

   printf("building query sequence...\n");
   /* get query sequence */
   SEQ *query_seq = (SEQ *)malloc( sizeof(SEQ) );
   printf("test...\n");
   seq_Create(query_seq, fasta_file);
   seq_Display(query_seq);
   int Q = query_seq->leng;

   printf("I/O was successful!\n");


   /* allocate memory to store results */
   float sc, perc_cells;
   int num_cells; 
   int tot_cells = (Q + 1) * (T + 1);
   RESULTS *res = (RESULTS *)malloc( sizeof(RESULTS) );
   TRACEBACK *tr = (TRACEBACK *)malloc( sizeof(TRACEBACK) );
   EDGEBOUNDS *edg_fwd = (EDGEBOUNDS *)malloc( sizeof(EDGEBOUNDS) );
   EDGEBOUNDS *edg_bck = (EDGEBOUNDS *)malloc( sizeof(EDGEBOUNDS) );
   EDGEBOUNDS *edg = (EDGEBOUNDS *)malloc( sizeof(EDGEBOUNDS) );

   /* allocate memory for square matrices (for testing) */
   float st_MX[ NUM_NORMAL_STATES * (Q+1) * (T+1) ];
   float st_MX3[ NUM_NORMAL_STATES * (Q+1) * 3 ];
   float sp_MX[ NUM_SPECIAL_STATES * (Q+1) ];


   // printf("=== TEST CYCLES ===\n");
   // test_cycle(Q, T, st_MX, sp_MX, tr);
   // dp_matrix_Print(Q, T, st_MX, sp_MX);
   // rev_test_cycle(Q, T, st_MX, sp_MX, tr);
   // dp_matrix_Print(Q, T, st_MX, sp_MX);
   // exit(0);

   /* run viterbi algorithm */
   printf("=== VITERBI -> START ===\n");
   sc = viterbi_Run(query_seq, target_prof, Q, T, st_MX, sp_MX, res, tr);
   printf("Viterbi Score: %.9f\n", sc);
   // dp_matrix_Print(Q, T, st_MX, sp_MX);
   dp_matrix_Save(Q, T, st_MX, sp_MX, "output/myversion.viterbi.tsv");
   printf("=== VITERBI -> END ===\n");

   /* run traceback of viterbi */
   printf("=== TRACEBACK -> START ===\n");
   traceback_Build(query_seq, target_prof, Q, T, st_MX, sp_MX, tr);
   traceback_Print(tr);
   traceback_Show(Q, T, st_MX, sp_MX, tr);
   // dp_matrix_Print(Q, T, st_MX, sp_MX);
   dp_matrix_Save(Q, T, st_MX, sp_MX, "output/myversion.traceback.tsv");
   printf("=== TRACEBACK -> END ===\n");

   /* run forward/backward algorithms */
   printf("=== FORWARD -> START ===\n");
   init_Logsum();
   dp_matrix_Clear(Q, T, st_MX, sp_MX);
   sc = forward_Run(query_seq, target_prof, Q, T, st_MX, sp_MX, res);
   printf("Forward Score: %.9f\n", sc);
   // dp_matrix_Print(Q, T, st_MX, sp_MX);
   dp_matrix_Save(Q, T, st_MX, sp_MX, "output/myversion.forward.tsv");
   printf("=== FORWARD -> END ===\n");

   printf("=== BACKWARD -> START ===\n");
   sc = backward_Run(query_seq, target_prof, Q, T, st_MX, sp_MX, res);
   printf("Backward Score: %.9f\n", sc);
   // dp_matrix_Print(Q, T, st_MX, sp_MX);
   dp_matrix_Save(Q, T, st_MX, sp_MX, "output/myversion.backward.tsv");
   printf("=== BACKWARD -> END ===\n");
   fflush(stdout);


   /* cloud search algorithm parameters */
   float alpha = 20.0;
   int beta = 5;

   /* use linear or quadratic space version? */
   bool is_linspace = false;

   if (is_linspace) 
   {
      /* cloud forward */
      printf("=== CLOUD FORWARD -> START ===\n");
      cloud_forward_Run3(query_seq, target_prof, Q, T, st_MX3, sp_MX, res, tr, edg_fwd, alpha, beta);
      // edgebounds_Print(edg_fwd);
      // dp_matrix_Print(Q, T, st_MX, sp_MX);
      // dp_matrix_Save(Q, T, st_MX, sp_MX, "output/myversion.cloud_fwd.tsv");
      printf("=== CLOUD FORWARD -> END ===\n");

      /* cloud backward */
      printf("=== CLOUD BACKWARD -> START ===\n");
      cloud_backward_Run3(query_seq, target_prof, Q, T, st_MX3, sp_MX, res, tr, edg_bck, alpha, beta);
      // edgebounds_Print(edg_bck);
      // dp_matrix_Print(Q, T, st_MX, sp_MX);
      // dp_matrix_Save(Q, T, st_MX, sp_MX, "output/myversion.cloud_bck.tsv");
      printf("=== CLOUD BACKWARD -> END ===\n");

      // printf("=== MERGE CLOUD -> START ===\n");
      // edgebounds_Merge(edg_bck, edg_fwd, edg);
      // printf("=== MERGE CLOUD -> END ===\n");

      // printf("=== REORIENT CLOUD -> START ===\n");
      // edgebounds_Reorient(edg, edg_row);
      // printf("=== REORIENT CLOUD -> END ===\n");

      printf("=== REORIENT/MERGE CLOUD -> START ===\n");
      num_cells = edgebounds_Merge_Reorient(edg_fwd, edg_bck, edg, Q, T, st_MX, sp_MX);
      perc_cells = (float)num_cells/(float)tot_cells;
      printf("Cells Computed = %d/%d = %.5f\n", num_cells, tot_cells, perc_cells);
      // edgebounds_Print(edg);
      printf("=== REORIENT/MERGE CLOUD -> END ===\n");

      /* Data viz to see cloud */
      // printf("=== TEST CLOUD -> START ===\n");
      // dp_matrix_Clear_X(Q, T, st_MX, sp_MX, 0);
      // test_cloud(Q, T, st_MX, sp_MX, edg_fwd, 1);
      // test_cloud(Q, T, st_MX, sp_MX, edg_bck, -1);
      // // dp_matrix_Print(Q, T, st_MX, sp_MX);
      // dp_matrix_Save(Q, T, st_MX, sp_MX, "output/myversion.cloud.tsv");
      // printf("=== TEST CLOUD -> END ===\n");

      /* bounded forward */
      printf("=== BOUNDED FORWARD -> START ===\n");
      dp_matrix_Clear_X(Q, T, st_MX, sp_MX, -INF);
      sc = forward_bounded_Run3(query_seq, target_prof, Q, T, st_MX3, sp_MX, res, edg);
      printf("Bounded Forward Score: %.5f\n", sc);
      // dp_matrix_Print(Q, T, st_MX, sp_MX);
      // dp_matrix_Save(Q, T, st_MX, sp_MX, "output/myversion.bounded_fwd.tsv");
      printf("=== BOUNDED FORWARD -> END ===\n");

      /* bounded backward */
      printf("=== BOUNDED BACKWARD -> START ===\n");
      dp_matrix_Clear_X(Q, T, st_MX, sp_MX, -INF);
      // sc = backward_bounded_Run3(query_seq, target_prof, Q, T, st_MX3, sp_MX, res, edg);
      printf("Bounded Backward Score: %.5f\n", sc);
      // dp_matrix_Print(Q, T, st_MX, sp_MX);
      // dp_matrix_Save(Q, T, st_MX, sp_MX, "output/myversion.bounded_bck.tsv");
      printf("=== BOUNDED BACKWARD -> END ===\n");

      // /* display results */
      // results_Display(results1);
   }
   else
   {
      /* cloud forward */
      printf("=== CLOUD FORWARD -> START ===\n");
      cloud_forward_Run(query_seq, target_prof, Q, T, st_MX, sp_MX, res, tr, edg_fwd, alpha, beta);
      // edgebounds_Print(edg_fwd);
      // dp_matrix_Print(Q, T, st_MX, sp_MX);
      dp_matrix_Save(Q, T, st_MX, sp_MX, "output/myversion.cloud_fwd.tsv");
      printf("=== CLOUD FORWARD -> END ===\n");

      /* cloud backward */
      printf("=== CLOUD BACKWARD -> START ===\n");
      cloud_backward_Run(query_seq, target_prof, Q, T, st_MX, sp_MX, res, tr, edg_bck, alpha, beta);
      // edgebounds_Print(edg_bck);
      // dp_matrix_Print(Q, T, st_MX, sp_MX);
      dp_matrix_Save(Q, T, st_MX, sp_MX, "output/myversion.cloud_bck.tsv");
      printf("=== CLOUD BACKWARD -> END ===\n");

      // printf("=== MERGE CLOUD -> START ===\n");
      // edgebounds_Merge(edg_bck, edg_fwd, edg);
      // printf("=== MERGE CLOUD -> END ===\n");

      // printf("=== REORIENT CLOUD -> START ===\n");
      // edgebounds_Reorient(edg, edg_row);
      // printf("=== REORIENT CLOUD -> END ===\n");

      printf("=== REORIENT/MERGE CLOUD -> START ===\n");
      num_cells = edgebounds_Merge_Reorient(edg_fwd, edg_bck, edg, Q, T, st_MX, sp_MX);
      perc_cells = (float)num_cells/(float)tot_cells;
      printf("Cells Computed = %d/%d = %.2f\n", num_cells, tot_cells, perc_cells);
      // edgebounds_Print(edg);
      printf("=== REORIENT/MERGE CLOUD -> END ===\n");

      /* Data viz to see cloud */
      // printf("=== TEST CLOUD -> START ===\n");
      // dp_matrix_Clear_X(Q, T, st_MX, sp_MX, 0);
      // test_cloud(Q, T, st_MX, sp_MX, edg_fwd, 1);
      // test_cloud(Q, T, st_MX, sp_MX, edg_bck, -1);
      // // dp_matrix_Print(Q, T, st_MX, sp_MX);
      // dp_matrix_Save(Q, T, st_MX, sp_MX, "output/myversion.cloud.tsv");
      // printf("=== TEST CLOUD -> END ===\n");

      /* bounded forward */
      printf("=== BOUNDED FORWARD -> START ===\n");
      dp_matrix_Clear_X(Q, T, st_MX, sp_MX, 0);
      sc = forward_bounded_Run(query_seq, target_prof, Q, T, st_MX, sp_MX, res, edg);
      printf("Bounded Forward Score: %.2f\n", sc);
      // dp_matrix_Print(Q, T, st_MX, sp_MX);
      dp_matrix_Save(Q, T, st_MX, sp_MX, "output/myversion.bounded_fwd.tsv");
      printf("=== BOUNDED FORWARD -> END ===\n");

      /* bounded backward */
      printf("=== BOUNDED BACKWARD -> START ===\n");
      dp_matrix_Clear_X(Q, T, st_MX, sp_MX, 0);
      sc = backward_bounded_Run(query_seq, target_prof, Q, T, st_MX, sp_MX, res, edg);
      printf("Bounded Backward Score: %.2f\n", sc);
      // dp_matrix_Print(Q, T, st_MX, sp_MX);
      dp_matrix_Save(Q, T, st_MX, sp_MX, "output/myversion.bounded_bck.tsv");
      printf("=== BOUNDED BACKWARD -> END ===\n");

      // /* display results */
      // results_Display(results1);
   }

   printf("...test finished. \n");
}


/* unit test */
void test2(char *hmm_file, char *fasta_file)
{
}