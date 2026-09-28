/**
 * @file varfreq_main.c
 * @brief entry point to varfreq
 * @author Suneth Samarasinghe (imsuneth@gmail.com)

MIT License

Copyright (c) 2026 Suneth Samarasinghe (imsuneth@gmail.com)

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.


******************************************************************************/

#include <getopt.h>
#include <pthread.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "error.h"
#include "misc.h"
#include "minimod.h"
#include "mod.h"
#include "ref.h"
#include "varmod.h"

static struct option long_options[] = {
    {"bedmethyl", no_argument, 0, 'b'},            //0 output in bedMethyl format
    {"mod_codes", required_argument, 0, 'c'},      //1 modification codes (eg. m, h or mh) [m]
    {"mod_thresh", required_argument, 0, 'm'},     //2 min modification threshold 0.0 to 1.0 [0.8]
    {"threads", required_argument, 0, 't'},        //3 number of threads [8]
    {"batchsize", required_argument, 0, 'K'},      //4 batchsize - number of reads loaded at once [512]
    {"max-bytes", required_argument, 0, 'B'},      //5 batchsize - number of bytes loaded at once
    {"verbose", required_argument, 0, 'v'},        //6 verbosity level [1]
    {"help", no_argument, 0, 'h'},                 //7
    {"version", no_argument, 0, 'V'},              //8
    {"prog-interval",required_argument, 0, 'p'},   //9 progress interval
    {"debug-break",required_argument, 0, 0},       //10 break after processing the first batch (used for debugging)
    {"output",required_argument, 0, 'o'},          //11 output file
    {"flank",required_argument, 0, 'w'},           //12 reference flank used when aligning to the ALT allele
    {"min-flank",required_argument, 0, 0},         //13 minimum reference flank a read has to provide
    {"band",required_argument, 0, 0},              //14 ksw2 band width
    {0, 0, 0, 0}};

static inline void print_help_msg(FILE *fp_help, opt_t opt){
    fprintf(fp_help,"Usage: minimod varfreq ref.fa reads.bam variants.vcf\n");
    fprintf(fp_help,"\nOutputs base modification frequencies at CpGs inside insertion structural variants.\n");
    fprintf(fp_help,"Reads are taken from the RNAMES field of each SVTYPE=INS record and are kept only when\n");
    fprintf(fp_help,"their HP tag matches the haplotype the phased genotype places the ALT on.\n");
    fprintf(fp_help,"\nbasic options:\n");
    fprintf(fp_help,"   -b                         output in bedMethyl format [%s]\n", (opt.bedmethyl_out?"yes":"not set"));
    fprintf(fp_help,"   -c STR                     modification code(s) (eg. m, h or mh or as ChEBI) [%s]\n", opt.mod_codes_str);
    fprintf(fp_help,"   -m FLOAT                   min modification threshold(s). Comma separated values for each modification code given in -c [%s]\n", opt.mod_threshes_str);
    fprintf(fp_help,"   -t INT                     number of processing threads [%d]\n",opt.num_thread);
    fprintf(fp_help,"   -K INT                     batch size (max number of reads loaded at once) [%d]\n",opt.batch_size);
    fprintf(fp_help,"   -B FLOAT[K/M/G]            max number of bases loaded at once [%.1fM]\n",opt.batch_size_bases/(float)(1000*1000));
    fprintf(fp_help,"   -w INT                     reference flank each side of the insertion when aligning [%d]\n",opt.flank);
    fprintf(fp_help,"   -h                         help\n");
    fprintf(fp_help,"   -p INT                     print progress every INT seconds (0: per batch) [%d]\n", opt.progress_interval);
    fprintf(fp_help,"   -o FILE                    output file [%s]\n", opt.output_file==NULL?"stdout":opt.output_file);
    fprintf(fp_help,"   --verbose INT              verbosity level [%d]\n",(int)get_log_level());
    fprintf(fp_help,"   --version                  print version\n");

    fprintf(fp_help,"\nadvanced options:\n");
    fprintf(fp_help,"   --min-flank INT            skip a read providing less reference flank than this [%d]\n",opt.min_flank);
    fprintf(fp_help,"   --band INT                 alignment band width on top of the length difference [%d]\n",opt.band);
    fprintf(fp_help,"   --debug-break INT          break after processing the specified no. of batches\n");
}

//function that processes a databatch - for pthreads when I/O and processing are interleaved
static void* varfreq_processor(void* voidargs) {
    pthread_arg2_t* args = (pthread_arg2_t*)voidargs;
    db_t* db = args->db;
    core_t* core = args->core;
    double realtime0=core->realtime0;

    double realtime_prog = realtime();

    process_db(core, db);

    if(core->opt.progress_interval<=0 || realtime()-realtime_prog > core->opt.progress_interval){
        fprintf(stderr, "[%s::%.3f*%.2f] %d Entries (%.1fM bytes) processed\n", __func__,
                realtime() - realtime0, cputime() / (realtime() - realtime0),
                (db->n_bam_recs), (db->total_bytes)/(1000.0*1000.0));
    }

    //need to inform the output thread that we completed the processing
    pthread_mutex_lock(&args->mutex);
    pthread_cond_signal(&args->cond);
    args->finished=1;
    pthread_mutex_unlock(&args->mutex);

    pthread_exit(0);
}

//function that merges the output and frees - for pthreads when I/O and processing are interleaved
static void* varfreq_post_processor(void* voidargs){
    pthread_arg2_t* args = (pthread_arg2_t*)voidargs;
    db_t* db = args->db;
    core_t* core = args->core;

    //wait until the processing thread has informed us
    pthread_mutex_lock(&args->mutex);
    while(args->finished==0){
        pthread_cond_wait(&args->cond, &args->mutex);
    }
    pthread_mutex_unlock(&args->mutex);

    merge_db(core, db);

    // unlike freq, varfreq is meant to skip almost every read in the bam, so there is
    // no warning here about the number of skipped reads

    free_db_tmp(core, db);
    free_db(core, db);
    free(args);
    pthread_exit(0);
}

/* varfreq finds its CpGs in the ALT sequence, so a context given in -c means nothing here */
static void warn_about_contexts(opt_t * opt){
    for (khint_t i = kh_begin(opt->modcodes_map); i < kh_end(opt->modcodes_map); ++i) {
        if (!kh_exist(opt->modcodes_map, i)) continue;
        modcodem_t * mod_code_map = kh_value(opt->modcodes_map, i);
        if (strcmp(mod_code_map->context, "CG") == 0) continue;
        WARNING("Context %s given for modification code %s is ignored. varfreq always uses CG.",
                mod_code_map->context, (char *)kh_key(opt->modcodes_map, i));
    }
}

int varfreq_main(int argc, char* argv[]) {

    double realtime0 = realtime();

    const char* optstring = "m:c:t:B:K:v:p:o:w:hVb";

    int longindex = 0;
    int32_t c = -1;

    FILE *fp_help = stderr;

    opt_t opt;
    init_opt(&opt); //initialise options to defaults
    opt.subtool = VARFREQ;

    //parse the user args
    while ((c = getopt_long(argc, argv, optstring, long_options, &longindex)) >= 0) {

        if (c == 'B') {
            opt.batch_size_bases = mm_parse_num(optarg);
            if(opt.batch_size_bases<=0){
                ERROR("%s","Maximum number of bases should be larger than 0.");
                exit(EXIT_FAILURE);
            }
        } else if (c == 'K') {
            opt.batch_size = atoi(optarg);
            if (opt.batch_size < 1) {
                ERROR("Batch size should larger than 0. You entered %d",opt.batch_size);
                exit(EXIT_FAILURE);
            }
        } else if (c == 't') {
            opt.num_thread = atoi(optarg);
            if (opt.num_thread < 1) {
                ERROR("Number of threads should larger than 0. You entered %d", opt.num_thread);
                exit(EXIT_FAILURE);
            }
        } else if (c == 'w') {
            opt.flank = atoi(optarg);
            if (opt.flank < 1) {
                ERROR("Flank should be larger than 0. You entered %d", opt.flank);
                exit(EXIT_FAILURE);
            }
        } else if (c=='v'){
            int v = atoi(optarg);
            set_log_level((enum log_level_opt)v);
        } else if (c=='p'){
            if (atoi(optarg) < 0) {
                ERROR("Progress interval should be 0 or positive. You entered %d", atoi(optarg));
                exit(EXIT_FAILURE);
            }
            opt.progress_interval = atoi(optarg);
        } else if (c=='o'){
            FILE *fp = fopen(optarg, "w");
            if (fp == NULL) {
                ERROR("Cannot open file %s for writing", optarg);
                exit(EXIT_FAILURE);
            }
            opt.output_file = optarg;
            opt.output_fp = fp;
        } else if (c=='V'){
            fprintf(stdout,"minimod %s\n",MINIMOD_VERSION);
            exit(EXIT_SUCCESS);
        } else if (c=='h'){
            fp_help = stdout;
        } else if (c=='m'){
            opt.mod_threshes_str = (char *)malloc(strlen(optarg)+1);
            MALLOC_CHK(opt.mod_threshes_str);
            strcpy(opt.mod_threshes_str,optarg);
        } else if (c=='c') {
            opt.mod_codes_str = optarg;
        } else if (c=='b'){
            opt.bedmethyl_out = 1;
        } else if(c == 0 && longindex == 10){ //debug break
            opt.debug_break = atoi(optarg);
        } else if(c == 0 && longindex == 13){ //min flank
            opt.min_flank = atoi(optarg);
            if (opt.min_flank < 1) {
                ERROR("Minimum flank should be larger than 0. You entered %d", opt.min_flank);
                exit(EXIT_FAILURE);
            }
        } else if(c == 0 && longindex == 14){ //band
            opt.band = atoi(optarg);
            if (opt.band < 1) {
                ERROR("Band width should be larger than 0. You entered %d", opt.band);
                exit(EXIT_FAILURE);
            }
        } else {
            print_help_msg(fp_help, opt);
            if(fp_help == stdout){
                exit(EXIT_SUCCESS);
            }
            exit(EXIT_FAILURE);
        }
    }

    if(opt.mod_codes_str==NULL || strlen(opt.mod_codes_str)==0){
        INFO("%s", "Modification codes not provided. Using default modification code m");
        opt.mod_codes_str = "m";
    }

    parse_mod_codes(&opt);
    warn_about_contexts(&opt);

    if(opt.mod_threshes_str==NULL || strlen(opt.mod_threshes_str)==0){
        INFO("%s", "Modification threshold not provided. Using default threshold 0.8");

        opt.mod_threshes_str = (char *)malloc(opt.n_mods * 4 * sizeof(char) + 1);
        MALLOC_CHK(opt.mod_threshes_str);
        memset(opt.mod_threshes_str,0,opt.n_mods * 4 * sizeof(char)+1);

        char * thresh_str = "0.8";

        for(int i=0;i<opt.n_mods;i++){
            strcat(opt.mod_threshes_str,thresh_str);
            if(i<opt.n_mods-1){
                strcat(opt.mod_threshes_str,",");
            }
        }
    }

    parse_mod_threshes(&opt);

    if (opt.min_flank > opt.flank) {
        ERROR("Minimum flank (%d) cannot be larger than the flank (%d)", opt.min_flank, opt.flank);
        exit(EXIT_FAILURE);
    }

    // three arguments are needed: reference, bam and vcf
    if (argc - optind != 3 || fp_help == stdout) {
        WARNING("%s","Missing arguments");
        print_help_msg(fp_help, opt);
        if(fp_help == stdout){
            exit(EXIT_SUCCESS);
        }
        exit(EXIT_FAILURE);
    }

    opt.ref_file = argv[optind];
    opt.bam_file = argv[optind+1];
    opt.vcf_file = argv[optind+2];

    if (access(opt.bam_file, F_OK) == -1) {
        ERROR("BAM file %s does not exist", opt.bam_file);
        exit(EXIT_FAILURE);
    }

    if (access(opt.vcf_file, F_OK) == -1) {
        ERROR("VCF file %s does not exist", opt.vcf_file);
        exit(EXIT_FAILURE);
    }

    //load the insertions and the read names that support them
    double realtime1 = realtime();
    fprintf(stderr, "[%s] Loading variants %s\n", __func__, opt.vcf_file);
    load_variants(opt.vcf_file);
    fprintf(stderr, "[%s] Variants loaded in %.3f sec\n", __func__, realtime()-realtime1);

    //load the reference genome. varfreq takes its CpG context from the ALT allele,
    //so there is no need to precompute the contexts in the reference
    double realtime2 = realtime();
    fprintf(stderr, "[%s] Loading reference genome %s\n", __func__, opt.ref_file);
    load_ref(opt.ref_file);
    fprintf(stderr, "[%s] Reference genome loaded in %.3f sec\n", __func__, realtime()-realtime2);

    check_variants_against_ref();

    //initialise the core data structure
    core_t* core = init_core(opt, realtime0);

    int32_t counter=0;

    print_varfreq_header(core);

    ret_status_t status = {core->opt.batch_size,core->opt.batch_size_bases};
    int8_t first_flag_p=0;
    int8_t first_flag_pp=0;
    pthread_t tid_p; //process thread
    pthread_t tid_pp; //post-process thread

    while (status.num_reads >= core->opt.batch_size || status.num_bases>=core->opt.batch_size_bases) {

        //init and load a databatch
        db_t* db = init_db(core);
        status = load_db(core, db);

        fprintf(stderr, "[%s::%.3f*%.2f] %d Entries (%.1fM bases) loaded\n", __func__,
                realtime() - realtime0, cputime() / (realtime() - realtime0),
                status.num_reads,status.num_bases/(1000.0*1000.0));

        if(first_flag_p){ //if not the first time of the "process" wait for the previous "process"
            int ret = pthread_join(tid_p, NULL);
            NEG_CHK(ret);
        }
        first_flag_p=1;

        //set up args
        pthread_arg2_t *pt_arg = (pthread_arg2_t*)malloc(sizeof(pthread_arg2_t));
        MALLOC_CHK(pt_arg);
        pt_arg->core=core;
        pt_arg->db=db;
        pthread_cond_init(&pt_arg->cond, NULL);
        pthread_mutex_init(&pt_arg->mutex, NULL);
        pt_arg->finished = 0;

        //process thread launch
        int ret = pthread_create(&tid_p, NULL, varfreq_processor, (void*)(pt_arg));
        NEG_CHK(ret);

        if(first_flag_pp){ //if not the first time of the post-process wait for the previous post-process
            int ret = pthread_join(tid_pp, NULL);
            NEG_CHK(ret);
        }
        first_flag_pp=1;

        //post-process thread launch (merging and freeing thread)
        ret = pthread_create(&tid_pp, NULL, varfreq_post_processor, (void*)(pt_arg));
        NEG_CHK(ret);

        if(opt.debug_break==counter){
            break;
        }
        counter++;
    }

    //final round
    int ret = pthread_join(tid_p, NULL);
    NEG_CHK(ret);
    ret = pthread_join(tid_pp, NULL);
    NEG_CHK(ret);

    output_core(core);

    print_varfreq_log();

    destroy_variants();
    destroy_ref(opt.n_mods);

    fprintf(stderr, "[%s] total entries: %ld", __func__,(long)core->total_reads);
    fprintf(stderr,"\n[%s] total bytes: %.1f M",__func__,core->total_bytes/(float)(1000*1000));
    fprintf(stderr,"\n[%s] total processed entries: %ld",__func__,(long)core->processed_reads);
    fprintf(stderr,"\n[%s] total processed bytes: %.1f M",__func__,(core->processed_bytes)/(float)(1000*1000));

    fprintf(stderr, "\n[%s] Data loading time: %.3f sec", __func__,core->load_db_time);
    fprintf(stderr, "\n[%s] Data processing time: %.3f sec", __func__,core->process_db_time);
    fprintf(stderr, "\n[%s] Data merging time: %.3f sec", __func__,core->merge_db_time);
    fprintf(stderr, "\n[%s] Data sorting time: %.3f sec", __func__,core->sort_time);
    fprintf(stderr, "\n[%s] Data output time: %.3f sec", __func__,core->output_time);

    fprintf(stderr,"\n");

    free_core(core,opt);
    free_opt(&opt);

    return 0;
}
