/**
 * @file varmod.h
 * @brief base modifications inside insertion structural variants
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

#ifndef VARMOD_H
#define VARMOD_H

#include "minimod.h"

/* one insertion variant taken from the VCF */
typedef struct {
    char * chrom;       // contig name
    int pos;            // 0-based position of the anchor base (VCF POS - 1)
    char * id;          // variant id from the VCF (eg. Sniffles2.INS.5B0S15)
    char * ins_seq;     // the inserted bases only (ALT without the leading anchor base)
    int ins_len;        // length of ins_seq
    char ref_base;      // the reference base at pos, as the VCF spells it
    int usable;         // 0 once the variant is dropped, eg. it disagrees with the reference
    int alt_hap;        // haplotype carrying the ALT. 1, 2, or 0 when hom-alt (both)
    int n_carrier;      // number of reads that were actually used for this variant
} var_t;

/* read the VCF and keep the insertion variants we can work with */
void load_variants(const char * vcf_file);

/* free everything load_variants allocated */
void destroy_variants();

/* drop variants whose contig or anchor base does not agree with the reference.
   call this after load_ref */
void check_variants_against_ref();

/* 1 if this read name supports at least one of the loaded variants */
int read_has_variants(const char * qname);

/* process a single read - fills db->freq_maps[bam_i] */
void varfreq_single(core_t * core, db_t * db, int32_t bam_i);

/* output */
void print_varfreq_header(core_t * core);
void print_varfreq_output(core_t * core);

/* print how many variants and reads were skipped and why */
void print_varfreq_log();

#endif
