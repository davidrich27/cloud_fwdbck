/*******************************************************************************
 *  @file main.c
 *  @brief Parses files, 
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

void parse_args (int argc, char *argv, ARGS *args);

/* MAIN */
int main (int argc, char *argv[])
{
   printf("Test begins... \n");

   ARGS *args = (ARGS *)malloc( sizeof(ARGS) );
   //parse_args(argc, *argv, args);

   /* load substitution matrix */
   // char *submat_file = "../data/submat/blosum62.submat";
   // SUBMAT *submat = (SUBMAT *)malloc( sizeof(SUBMAT) );
   // submat_Create(submat, submat_file);
   // submat_Display(submat);

   /* TEST 1 */
   char *hmm_file = "../data/test1_2.hmm";
   char *fasta_file = "../data/test1_1.fa";

   /* get target profile */
   HMM_PROFILE *target_prof1 = (HMM_PROFILE *)malloc( sizeof(HMM_PROFILE) );
   hmmprofile_Create(target_prof1, hmm_file);
   // hmmprofile_Display(target_prof1);
   hmmprofile_Config(target_prof1);
   // hmmprofile_Display(target_prof1);
   int T = target_prof1->leng;

   /* get query sequence */
   SEQ *query_seq1 = (SEQ *)malloc( sizeof(SEQ) );
   seq_Create(query_seq1, fasta_file);
   seq_Display(query_seq1);
   int Q = query_seq1->leng;

   /* allocate memory to store results */
   RESULTS *results1 = (RESULTS *)malloc( sizeof(RESULTS) );
   TRACEBACK *trace1 = (TRACEBACK *)malloc( sizeof(TRACEBACK) );
   float st_MX[ NUM_NORMAL_STATES * (Q+1) * (T+1) ];
   float sp_MX[ NUM_SPECIAL_STATES * (Q+1) ];

   /* run viterbi and backtrace */
   viterbi_Run(query_seq1, target_prof1, Q, T, st_MX, sp_MX, results1, trace1);
   // viterbi_Traceback(query_seq1, target_prof1, Q, T, st_MX, sp_MX, results1, trace1);
   printf("=== VITERBI RESULTS ===\n");
   dp_matrix_Print(Q, T, st_MX, sp_MX);
   dp_matrix_Save(Q, T, st_MX, sp_MX, "viterbi.tsv");

   /* run forward/backward algorithms */
   dp_matrix_Clear(Q, T, st_MX, sp_MX);
   forward_Run(query_seq1, target_prof1, Q, T, st_MX, sp_MX, results1);
   printf("=== FORWARD RESULTS ===\n");
   dp_matrix_Print(Q, T, st_MX, sp_MX);
   dp_matrix_Save(Q, T, st_MX, sp_MX, "forward.tsv");

   backward_Run(query_seq1, target_prof1, Q, T, st_MX, sp_MX, results1);
   printf("=== BACKWARD RESULTS ===\n");
   dp_matrix_Print (Q, T, st_MX, sp_MX);

   /* run cloud search algorithms */
   cloud_search_forward_Run(query_seq1, target_prof1, Q, T, st_MX, sp_MX, results1, trace1);

   /* display results */
   results_Display(results1);

   /* TEST 2 */
   char *fasta_file1 = "../data/test1_1.fa";
   char *fasta_file2 = "../data/test1_2.fa";

   /* get target sequence */
   SEQ *target_seq2 = (SEQ *)malloc( sizeof(SEQ) );
   seq_Create(target_seq2, fasta_file1);
   seq_Display(target_seq2);

   /* get query sequence */
   SEQ *query_seq2 = (SEQ *)malloc( sizeof(SEQ) );
   seq_Create(query_seq2, fasta_file2);
   seq_Display(query_seq2);

   printf("Test finished... \n");
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
void test()
{

}

