/**
 * @file varmod.c
 * @brief base modification frequencies inside insertion structural variants
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <assert.h>
#include <pthread.h>

#include <htslib/vcf.h>
#include <htslib/sam.h>

#include "ksw2.h"
/* ksw2.h defines kmalloc/kcalloc/krealloc/kfree taking a memory pool argument.
   khash.h wants the plain malloc versions of those names, so drop ksw2's here.
   ksw2 only uses them inside ksw2_extd.c anyway. */
#undef kmalloc
#undef kcalloc
#undef krealloc
#undef kfree

#include "minimod.h"
#include "varmod.h"
#include "mod.h"
#include "ref.h"
#include "error.h"
#include "misc.h"
#include "ksort.h"

#define WILDCARD_STR "*"
#define THRESH_UINT8_TO_DBL(x) ((double)( (x + 0.5) / 256.0 ))
#define IS_DIGIT(c) ((c) >= '0' && (c) <= '9')
#define IS_ALPHA(c) (((c) >= 'a' && (c) <= 'z') || ((c) >= 'A' && (c) <= 'Z'))

/* the two-piece affine gap penalties minimap2 uses for ONT reads */
#define KSW_MATCH 2
#define KSW_MISMATCH 4
#define KSW_GAPO 4
#define KSW_GAPE 2
#define KSW_GAPO2 24
#define KSW_GAPE2 1

/* defined in mod.c */
uint8_t get_hp_tag(bam1_t *record);

static const int valid_bases[256] = { ['A'] = 1, ['C'] = 1, ['G'] = 1, ['T'] = 1, ['U'] = 1, ['N'] = 1, ['a'] = 1, ['c'] = 1, ['g'] = 1, ['t'] = 1, ['u'] = 1, ['n'] = 1 };
static const int valid_strands[256] = { ['+'] = 1, ['-'] = 1 };
static const int base_idx_lookup[256] = { ['A'] = 0, ['C'] = 1, ['G'] = 2, ['T'] = 3, ['U'] = 3, ['N'] = 4, ['a'] = 0, ['c'] = 1, ['g'] = 2, ['t'] = 3, ['u'] = 3, ['n'] = 4 };
static const char base_complement_lookup[256] = { ['A'] = 'T', ['C'] = 'G', ['G'] = 'C', ['T'] = 'A', ['U'] = 'A', ['N'] = 'N', ['a'] = 't', ['c'] = 'g', ['g'] = 'c', ['t'] = 'a', ['u'] = 'a', ['n'] = 'n' };

/* the variants we loaded from the VCF */
static var_t * variants = NULL;
static int n_variants = 0;
static int cap_variants = 0;

/* the variant indices a read name supports */
typedef struct {
    int * idx;
    int n;
    int cap;
} var_list_t;

KHASH_MAP_INIT_STR(readm, var_list_t *)
static khash_t(readm) * read_map = NULL;

/* counters filled while reading the VCF (single threaded) */
static int64_t n_vcf_total = 0;
static int64_t n_vcf_not_pass = 0;
static int64_t n_vcf_not_ins = 0;
static int64_t n_vcf_bad_alt = 0;
static int64_t n_vcf_bad_gt = 0;
static int64_t n_vcf_no_rnames = 0;
static int64_t n_vcf_no_contig = 0;
static int64_t n_vcf_ref_mismatch = 0;

/* counters filled while processing reads (many threads) */
static pthread_mutex_t var_mutex = PTHREAD_MUTEX_INITIALIZER;
static int64_t n_read_wrong_contig = 0;
static int64_t n_read_no_hp = 0;
static int64_t n_read_hp_mismatch = 0;
static int64_t n_read_short_flank = 0;
static int64_t n_read_used = 0;
static int64_t n_site_clash = 0;

/* encode a base the way ksw2 wants it: A=0 C=1 G=2 T=3 anything else=4 */
static inline uint8_t encode_base(char c) {
    if (c == 'A' || c == 'a') return 0;
    if (c == 'C' || c == 'c') return 1;
    if (c == 'G' || c == 'g') return 2;
    if (c == 'T' || c == 't') return 3;
    return 4;
}


/* Check if ALT is a  symbolic ALT or or other non-derivable format */
static int alt_not_derivable(const char * alt, int alt_len) {
    if (alt_len == 0) return 1;                             // no sequence to substitute
    if (alt_len == 1 && alt[0] == '*') return 1;            // allele missing due to an upstream deletion
    if (strchr(alt, '<') || strchr(alt, '>')) return 1;     // symbolic: <INS>, <DEL>, and the C<ctg1> shorthand
    if (strchr(alt, '[') || strchr(alt, ']')) return 1;     // mated breakend, eg. G]1:10000]
    if (alt[0] == '.' || alt[alt_len - 1] == '.') return 1; // single breakend, eg. G. or .TGCA
    return 0;
}

/* which haplotype carries the ALT?
   returns 1 or 2 for a phased het, 0 when hom-alt (both), -1 when we cannot tell */
static int get_alt_hap(bcf_hdr_t * hdr, bcf1_t * rec, int32_t ** gt, int * gt_cap) {
    int n = bcf_get_genotypes(hdr, rec, gt, gt_cap);
    if (n != 2) return -1; // no GT field (n < 0) or not diploid. the sample count is checked in load_variants

    int32_t * g = *gt;
    if (g[0] == bcf_int32_vector_end || g[1] == bcf_int32_vector_end) return -1;
    if (bcf_gt_is_missing(g[0]) || bcf_gt_is_missing(g[1])) return -1;

    int a0 = bcf_gt_allele(g[0]);
    int a1 = bcf_gt_allele(g[1]);

    if (a0 == 1 && a1 == 1) return 0; // hom-alt, both haplotypes carry the insertion

    // a het is only usable when it is phased, otherwise we cannot pick a haplotype
    if (!bcf_gt_is_phased(g[1])) return -1;
    if (a0 == 1 && a1 == 0) return 1;
    if (a0 == 0 && a1 == 1) return 2;

    return -1;
}

/* add one variant to the array and return its index */
static int add_variant(const char * chrom, int pos, const char * id, const char * ins_seq, int ins_len, int alt_hap) {
    if (n_variants == cap_variants) {
        cap_variants = cap_variants ? cap_variants * 2 : 1024;
        variants = (var_t *)realloc(variants, sizeof(var_t) * cap_variants);
        MALLOC_CHK(variants);
    }

    var_t * var = &variants[n_variants];

    var->chrom = (char *)malloc(strlen(chrom) + 1);
    MALLOC_CHK(var->chrom);
    strcpy(var->chrom, chrom);

    var->id = (char *)malloc(strlen(id) + 1);
    MALLOC_CHK(var->id);
    strcpy(var->id, id);

    var->ins_seq = (char *)malloc(ins_len + 1);
    MALLOC_CHK(var->ins_seq);
    for (int i = 0; i < ins_len; i++) { // the reference is uppercased, match it
        var->ins_seq[i] = toupper((unsigned char)ins_seq[i]);
    }
    var->ins_seq[ins_len] = '\0';

    var->pos = pos;
    var->ins_len = ins_len;
    var->alt_hap = alt_hap;
    var->n_carrier = 0;

    n_variants++;
    return n_variants - 1;
}

/* record that this read name supports this variant */
static void add_read_name(const char * name, int var_idx) {
    int ret;
    khiter_t k = kh_get(readm, read_map, name);

    if (k == kh_end(read_map)) {
        char * key = (char *)malloc(strlen(name) + 1);
        MALLOC_CHK(key);
        strcpy(key, name);

        k = kh_put(readm, read_map, key, &ret);

        var_list_t * vl = (var_list_t *)malloc(sizeof(var_list_t));
        MALLOC_CHK(vl);
        vl->cap = 2;
        vl->n = 0;
        vl->idx = (int *)malloc(sizeof(int) * vl->cap);
        MALLOC_CHK(vl->idx);

        kh_value(read_map, k) = vl;
    }

    var_list_t * vl = kh_value(read_map, k);
    if (vl->n == vl->cap) {
        vl->cap *= 2;
        vl->idx = (int *)realloc(vl->idx, sizeof(int) * vl->cap);
        MALLOC_CHK(vl->idx);
    }
    vl->idx[vl->n] = var_idx;
    vl->n++;
}

void load_variants(const char * vcf_file) {

    htsFile * fp = bcf_open(vcf_file, "r");
    if (fp == NULL) {
        ERROR("Could not open the VCF file %s", vcf_file);
        exit(EXIT_FAILURE);
    }

    bcf_hdr_t * hdr = bcf_hdr_read(fp);
    NULL_CHK(hdr);

    // we only support single sample VCFs
    if (bcf_hdr_nsamples(hdr) != 1) {
        ERROR("VCF file %s has %d samples. This program only supports a single sample VCF.", vcf_file, bcf_hdr_nsamples(hdr));
        exit(EXIT_FAILURE);
    }

    bcf1_t * rec = bcf_init();
    NULL_CHK(rec);

    read_map = kh_init(readm);

    char * svtype = NULL;
    int svtype_cap = 0;
    char * rnames = NULL;
    int rnames_cap = 0;
    int32_t * gt = NULL;
    int gt_cap = 0;

    while (bcf_read(fp, hdr, rec) >= 0) {

        n_vcf_total++;
        bcf_unpack(rec, BCF_UN_ALL);

        // 1. only PASS records
        if (bcf_has_filter(hdr, rec, "PASS") != 1) {
            n_vcf_not_pass++;
            continue;
        }

        // 2. only insertions
        if (bcf_get_info_string(hdr, rec, "SVTYPE", &svtype, &svtype_cap) < 0 || strcmp(svtype, "INS") != 0) {
            n_vcf_not_ins++;
            continue;
        }

        // 3. we need a single ALT written out as plain sequence
        if (rec->n_allele != 2) {
            n_vcf_bad_alt++;
            continue;
        }
        char * ref_allele = rec->d.allele[0];
        char * alt_allele = rec->d.allele[1];
        int ref_len = strlen(ref_allele);
        int alt_len = strlen(alt_allele);
        if (alt_not_derivable(alt_allele, alt_len) || ref_len != 1 || alt_len < 2) {
            n_vcf_bad_alt++;
            continue;
        }

        // 4. the genotype has to tell us which haplotype carries the insertion
        int alt_hap = get_alt_hap(hdr, rec, &gt, &gt_cap);
        if (alt_hap < 0) {
            n_vcf_bad_gt++;
            continue;
        }

        // 5. the supporting read names
        if (bcf_get_info_string(hdr, rec, "RNAMES", &rnames, &rnames_cap) < 0) {
            n_vcf_no_rnames++;
            continue;
        }

        // the first ALT base is the anchor base, the insertion is everything after it
        const char * chrom = bcf_hdr_id2name(hdr, rec->rid);

        // 6. the record has to agree with the reference we were given.
        ref_t * ref = get_ref(chrom);
        if (ref == NULL) {
            n_vcf_no_contig++;
            continue;
        }
        if (rec->pos < 0 || rec->pos >= ref->ref_seq_length) {
            n_vcf_ref_mismatch++;
            continue;
        }
        if (ref->forward[rec->pos] != toupper((unsigned char)ref_allele[0])) {
            n_vcf_ref_mismatch++;
            continue;
        }

        const char * id = (rec->d.id && rec->d.id[0] != '.') ? rec->d.id : "*";
        int var_idx = add_variant(chrom, rec->pos, id, alt_allele + 1, alt_len - 1, alt_hap);

        // RNAMES is a comma separated list
        char * name = strtok(rnames, ",");
        while (name != NULL) {
            add_read_name(name, var_idx);
            name = strtok(NULL, ",");
        }
    }

    free(svtype);
    free(rnames);
    free(gt);
    bcf_destroy(rec);
    bcf_hdr_destroy(hdr);
    bcf_close(fp);

    if (n_vcf_no_contig > 0) {
        WARNING("%ld variants are on a contig the reference does not have. They are skipped.", (long)n_vcf_no_contig);
    }
    if (n_vcf_ref_mismatch > 0) {
        WARNING("%ld variants do not match the reference base at their position. They are skipped. Is this the reference the VCF was called against?", (long)n_vcf_ref_mismatch);
    }

    INFO("Loaded %d insertion variants supported by %d read names", n_variants, kh_size(read_map));

    if (n_variants == 0) {
        WARNING("%s", "No usable insertion variants found in the VCF. The output will be empty.");
    }
}

int read_has_variants(const char * qname) {
    if (read_map == NULL) return 0;
    return kh_get(readm, read_map, qname) != kh_end(read_map);
}

void destroy_variants() {

    for (int i = 0; i < n_variants; i++) {
        free(variants[i].chrom);
        free(variants[i].id);
        free(variants[i].ins_seq);
    }
    free(variants);
    variants = NULL;
    n_variants = 0;
    cap_variants = 0;

    if (read_map != NULL) {
        for (khiter_t k = kh_begin(read_map); k != kh_end(read_map); ++k) {
            if (kh_exist(read_map, k)) {
                char * key = (char *)kh_key(read_map, k);
                var_list_t * vl = kh_value(read_map, k);
                free(vl->idx);
                free(vl);
                free(key);
            }
        }
        kh_destroy(readm, read_map);
        read_map = NULL;
    }
}

/* key: chrom \t pos \t var_idx \t offset \t strand \t mod_code */
static char * make_var_key(const char * chrom, int pos, int var_idx, int offset, char strand, const char * mod_code) {
    int len = strlen(chrom) + strlen(mod_code) + 48;
    char * key = (char *)malloc(len);
    MALLOC_CHK(key);
    snprintf(key, len, "%s\t%d\t%d\t%d\t%c\t%s", chrom, pos, var_idx, offset, strand, mod_code);
    return key;
}

static void decode_var_key(char * key, char ** chrom, int * pos, int * var_idx, int * offset, char * strand, char ** mod_code) {
    char * token = strtok(key, "\t");
    *chrom = (char *)malloc(strlen(token) + 1);
    MALLOC_CHK(*chrom);
    strcpy(*chrom, token);

    *pos = atoi(strtok(NULL, "\t"));
    *var_idx = atoi(strtok(NULL, "\t"));
    *offset = atoi(strtok(NULL, "\t"));
    *strand = strtok(NULL, "\t")[0];

    token = strtok(NULL, "\t");
    *mod_code = (char *)malloc(strlen(token) + 1);
    MALLOC_CHK(*mod_code);
    strcpy(*mod_code, token);
}

/* every row of one insertion shares the same contig and position, so we have to
   keep sorting on the variant, the offset, the strand and the code as well */
static int cmp_var_key(const char * a, const char * b) {

    // 1. contig
    const char * ea = strchr(a, '\t');
    const char * eb = strchr(b, '\t');
    size_t la = ea ? (size_t)(ea - a) : strlen(a);
    size_t lb = eb ? (size_t)(eb - b) : strlen(b);
    size_t min_len = la < lb ? la : lb;
    int cmp = strncmp(a, b, min_len);
    if (cmp != 0) return cmp;
    if (la != lb) return la < lb ? -1 : 1;
    if (!ea || !eb) return 0;

    // 2. position, then variant, then offset
    const char * pa = ea;
    const char * pb = eb;
    for (int i = 0; i < 3; i++) {
        long va = strtol(pa + 1, (char **)&pa, 10);
        long vb = strtol(pb + 1, (char **)&pb, 10);
        if (va != vb) return va < vb ? -1 : 1;
    }

    // 3. strand and modification code
    return strcmp(pa, pb);
}

#define var_kv_lt(a, b) (cmp_var_key((a).key, (b).key) < 0)
KSORT_INIT(varfreq, freq_kv_t, var_kv_lt)

static void update_var_map(khash_t(freqm) * freq_map, const char * chrom, int pos, int var_idx, int offset, char strand, const char * mod_code, int is_called, int is_mod) {

    char * key = make_var_key(chrom, pos, var_idx, offset, strand, mod_code);

    khiter_t k = kh_get(freqm, freq_map, key);
    if (k == kh_end(freq_map)) { // first time we see this site
        freq_t * freq = (freq_t *)malloc(sizeof(freq_t));
        MALLOC_CHK(freq);
        freq->n_called = is_called;
        freq->n_mod = is_mod;

        int ret;
        k = kh_put(freqm, freq_map, key, &ret);
        kh_value(freq_map, k) = freq;
    } else {
        free(key);
        freq_t * freq = kh_value(freq_map, k);
        freq->n_called += is_called;
        freq->n_mod += is_mod;
    }
}


/* walk the CIGAR and find the read positions aligned to ref_lo and ref_hi.
   lo gets the first aligned read base at or after ref_lo,
   hi gets the last aligned read base at or before ref_hi. both are -1 if not found. */
static void find_read_window(bam1_t * record, int ref_lo, int ref_hi, int * lo, int * hi) {

    uint32_t * cigar = bam_get_cigar(record);
    uint32_t n_cigar = record->core.n_cigar;
    int read_pos = 0;
    int ref_pos = record->core.pos;

    *lo = -1;
    *hi = -1;

    for (uint32_t ci = 0; ci < n_cigar; ci++) {
        int len = bam_cigar_oplen(cigar[ci]);
        int op = bam_cigar_op(cigar[ci]);
        int type = bam_cigar_type(op);
        int consume_read = type & 1;
        int consume_ref = type & 2;

        if (op == BAM_CHARD_CLIP) {
            ERROR("Hard clipping found in %s and they are not supported.\nTry following workarounds.\n\t01. Filter out non-primary alignments\n\t\tsamtools view -h -F 2308 reads.bam -o primary_reads.bam\n\t02. Use minimap2 with -Y to use soft clipping for suplimentary alignments.\n", bam_get_qname(record));
            exit(EXIT_FAILURE);
        }

        if (consume_read && consume_ref) { // M, = or X
            for (int j = 0; j < len; j++) {
                int r = ref_pos + j;
                if (r < ref_lo) continue;
                if (r > ref_hi) break;
                if (*lo < 0) *lo = read_pos + j;
                *hi = read_pos + j;
            }
        }

        if (consume_read) read_pos += len;
        if (consume_ref) ref_pos += len;

        if (ref_pos > ref_hi) break;
    }
}

/* build the ALT allele around one variant, align the read to it and mark on the read
   which positions fall on a CpG of the insertion.
   returns 1 when the read was used, 0 when it was skipped. */
static int mark_sites_for_variant(core_t * core, bam1_t * record, ref_t * ref, int var_idx,
                                  int * site_var, int * site_off, char * site_base) {

    var_t * var = &variants[var_idx];
    int flank = core->opt.flank;
    int min_flank = core->opt.min_flank;

    // 1. how much reference can we use either side? the read has to cover it
    int read_ref_start = record->core.pos;
    int read_ref_end = bam_endpos(record) - 1; // last reference base the read touches

    int ref_lo = var->pos - flank;
    if (ref_lo < read_ref_start) ref_lo = read_ref_start;
    if (ref_lo < 0) ref_lo = 0;

    int ref_hi = var->pos + flank;
    if (ref_hi > read_ref_end) ref_hi = read_ref_end;
    if (ref_hi > ref->ref_seq_length - 1) ref_hi = ref->ref_seq_length - 1;

    int w5 = var->pos - ref_lo; // reference bases before the anchor base
    int w3 = ref_hi - var->pos; // reference bases after the anchor base

    if (w5 < min_flank || w3 < min_flank) {
        pthread_mutex_lock(&var_mutex);
        n_read_short_flank++;
        pthread_mutex_unlock(&var_mutex);
        return 0;
    }

    // 2. the matching region of the read
    int read_lo, read_hi;
    find_read_window(record, ref_lo, ref_hi, &read_lo, &read_hi);
    if (read_lo < 0 || read_hi <= read_lo) {
        pthread_mutex_lock(&var_mutex);
        n_read_short_flank++;
        pthread_mutex_unlock(&var_mutex);
        return 0;
    }

    // 3. the ALT allele: left reference flank, the inserted bases, right reference flank
    int t_len = (w5 + 1) + var->ins_len + w3;
    char * target_char = (char *)malloc(t_len);
    MALLOC_CHK(target_char);

    for (int i = 0; i <= w5; i++) {
        target_char[i] = ref->forward[ref_lo + i];
    }
    for (int i = 0; i < var->ins_len; i++) {
        target_char[w5 + 1 + i] = var->ins_seq[i];
    }
    for (int i = 0; i < w3; i++) {
        target_char[w5 + 1 + var->ins_len + i] = ref->forward[var->pos + 1 + i];
    }

    // 4. the CpGs we report. for offset o the C sits at target index w5+o and the G right after it.
    //    o goes from 0 (C is the reference anchor base) to ins_len (C is the last inserted base).
    int * target_site = (int *)malloc(sizeof(int) * t_len);
    MALLOC_CHK(target_site);
    for (int i = 0; i < t_len; i++) target_site[i] = -1;

    for (int o = 0; o <= var->ins_len; o++) {
        int tc = w5 + o;
        int tg = tc + 1;
        if (tg >= t_len) break;
        if (target_char[tc] == 'C' && target_char[tg] == 'G') {
            target_site[tc] = o;     // the C, reported on the + strand at offset o
            target_site[tg] = o + 1; // the G, reported on the - strand at offset o+1
        }
    }

    // 5. align the read stretch to the ALT allele
    int q_len = read_hi - read_lo + 1;
    uint8_t * seq = bam_get_seq(record);

    uint8_t * target = (uint8_t *)malloc(t_len);
    MALLOC_CHK(target);
    for (int i = 0; i < t_len; i++) target[i] = encode_base(target_char[i]);

    uint8_t * query = (uint8_t *)malloc(q_len);
    MALLOC_CHK(query);
    for (int i = 0; i < q_len; i++) query[i] = encode_base(seq_nt16_str[bam_seqi(seq, read_lo + i)]);

    int8_t mat[25];
    /*
        A    C    G    T    N
    A  +2   -4   -4   -4    0
    C  -4   +2   -4   -4    0
    G  -4   -4   +2   -4    0
    T  -4   -4   -4   +2    0
    N   0    0    0    0    0
    */
    for (int i = 0; i < 5; i++) {
        for (int j = 0; j < 5; j++) {
            if (i == 4 || j == 4) mat[i * 5 + j] = 0; // N matches nothing in particular
            else mat[i * 5 + j] = (i == j) ? KSW_MATCH : -KSW_MISMATCH;
        }
    }

    // band the alignment, otherwise the backtrack matrix is q_len*t_len bytes
    int band = core->opt.band + abs(t_len - q_len);

    ksw_extz_t ez;
    memset(&ez, 0, sizeof(ksw_extz_t));
    ksw_extd(0, q_len, query, t_len, target, 5, mat, KSW_GAPO, KSW_GAPE, KSW_GAPO2, KSW_GAPE2, band, -1, 0, &ez);

    // 6. walk the alignment and carry the CpG sites over to read positions
    int used = 0;
    int clash = 0;
    int qi = 0;
    int ti = 0;
    for (int c = 0; c < ez.n_cigar; c++) {
        int len = ez.cigar[c] >> 4;
        int op = ez.cigar[c] & 0xf;

        if (op == 0) { // in both the read and the ALT
            for (int j = 0; j < len; j++) {
                int t = ti + j;
                if (target_site[t] < 0) continue; // not a CpG site
                int read_pos = read_lo + qi + j;

                // two insertions this close would both claim the same read base
                if (site_var[read_pos] >= 0 && site_var[read_pos] != var_idx) clash++;
                
                site_var[read_pos] = var_idx;
                site_off[read_pos] = target_site[t];
                site_base[read_pos] = target_char[t];
                used = 1;
            }
            ti += len;
            qi += len;
        } else if (op == 1) { // in the read but not in the ALT
            qi += len;
        } else if (op == 2) { // in the ALT but not in the read
            ti += len;
        }
    }

    free(ez.cigar);
    free(query);
    free(target);
    free(target_site);
    free(target_char);

    if (used) {
        pthread_mutex_lock(&var_mutex);
        var->n_carrier++;
        n_read_used++;
        n_site_clash += clash;
        pthread_mutex_unlock(&var_mutex);
    }

    return used;
}

/* is this modification code one the user asked for? NULL when it is not */
static modcodem_t * get_required_mod(core_t * core, char * mod_code) {
    khiter_t mk = kh_get(modcodesm, core->opt.modcodes_map, WILDCARD_STR);
    if (mk != kh_end(core->opt.modcodes_map)) { // wildcard, everything is wanted
        return kh_value(core->opt.modcodes_map, mk);
    }
    mk = kh_get(modcodesm, core->opt.modcodes_map, mod_code);
    if (mk == kh_end(core->opt.modcodes_map)) return NULL;
    return kh_value(core->opt.modcodes_map, mk);
}

/* count one base call against the site it landed on */
static void count_one(db_t * db, int32_t bam_i, char strand, int read_pos,
                      int * site_var, int * site_off,
                      char * mod_code, modcodem_t * req_mod, uint8_t mod_prob) {

    int var_idx = site_var[read_pos];
    var_t * var = &variants[var_idx];

    uint8_t is_called = 0;
    uint8_t is_mod = 0;
    double thresh = req_mod->thresh;
    double mod_prob_dbl = THRESH_UINT8_TO_DBL(mod_prob);

    if (mod_prob_dbl >= thresh) {
        is_called = 1;
        is_mod = 1;
    } else if (mod_prob_dbl <= 1 - thresh) {
        is_called = 1;
    } else {
        return; // ambiguous, not a call either way
    }

    update_var_map(db->freq_maps[bam_i], var->chrom, var->pos, var_idx, site_off[read_pos], strand, mod_code, is_called, is_mod);
}

/* walk the MM/ML tags once and count every call that landed on a marked site */
static void count_calls(core_t * core, db_t * db, int32_t bam_i, int * site_var, int * site_off, char * site_base) {

    bam1_t * record = db->bam_recs[bam_i];
    int8_t rev = bam_is_rev(record);
    uint8_t * seq = bam_get_seq(record);
    int seq_len = record->core.l_qseq;
    char strand = rev ? '-' : '+';

    const char * mm_string = db->mm[bam_i];
    uint32_t ml_len = db->ml_lens[bam_i];
    uint8_t * ml = db->ml[bam_i];

    int ** bases_pos = db->bases_pos[bam_i];
    int bases_pos_lens[N_BASES] = {0};
    int * skip_counts = db->skip_counts[bam_i];
    char * mod_codes = db->mod_codes[bam_i];

    memset(db->mod_codes[bam_i], 0, core->opt.n_mods);

    // where every base sits in the read, so we can turn an MM skip count into a read position
    for (int i = 0; i < seq_len; i++) {
        int base_char = seq_nt16_str[bam_seqi(seq, i)];
        int idx = base_idx_lookup[base_char];
        bases_pos[idx][bases_pos_lens[idx]++] = i;
    }

    int mm_str_len = strlen(mm_string);
    int i = 0;
    int ml_start_idx = 0;

    while (i < mm_str_len) {

        int skip_counts_len = 0;
        int mod_codes_len = 0;
        char status_flag = '.';
        char modbase;

        // the base this modification applies to
        ASSERT_MSG(valid_bases[(int)mm_string[i]], "Invalid base:%c\n", mm_string[i]);
        modbase = mm_string[i] == 'U' ? 'T' : mm_string[i];
        i++;

        // the strand, which we do not use
        ASSERT_MSG(valid_strands[(int)mm_string[i]], "Invalid strand:%c\n", mm_string[i]);
        i++;

        // the modification codes
        int j = 0;
        int has_nums = 0;
        int has_alpha = 0;
        while (i < mm_str_len && mm_string[i] != ',' && mm_string[i] != ';' && mm_string[i] != '?' && mm_string[i] != '.') {
            if (IS_DIGIT(mm_string[i])) {
                has_nums = 1;
            } else if (IS_ALPHA(mm_string[i])) {
                has_alpha = 1;
            } else {
                ERROR("Invalid base modification code:%c. Modification codes should be either numeric or alphabetic.\n", mm_string[i]);
                exit(EXIT_FAILURE);
            }

            if (j >= db->mod_codes_cap[bam_i]) {
                db->mod_codes_cap[bam_i] *= 2;
                db->mod_codes[bam_i] = (char *)realloc(db->mod_codes[bam_i], sizeof(char) * (db->mod_codes_cap[bam_i] + 1));
                MALLOC_CHK(db->mod_codes[bam_i]);
            }
            mod_codes = db->mod_codes[bam_i];

            mod_codes[j] = mm_string[i];
            j++;
            i++;
        }
        mod_codes[j] = '\0';
        mod_codes_len = j;

        if (has_nums) mod_codes_len = 1; // a ChEBI id is a single code

        ASSERT_MSG(mod_codes_len > 0, "Invalid modification codes:%s. Modification codes cannot be empty.\n", mod_codes);
        ASSERT_MSG((has_nums && has_alpha) == 0, "Invalid modification codes:%s. Modification codes should be either numeric or alphabetic, not both.\n", mod_codes);

        // the status flag
        if (i < mm_str_len && (mm_string[i] == '?' || mm_string[i] == '.')) {
            status_flag = mm_string[i];
            i++;
        } else {
            status_flag = '.';
        }

        // the skip counts
        int k = 0;
        while (i < mm_str_len && mm_string[i] != ';') {
            if (mm_string[i] == ',') {
                i++;
                continue;
            }

            char skip_count_str[10];
            int l = 0;
            while (i < mm_str_len && mm_string[i] != ',' && mm_string[i] != ';') {
                skip_count_str[l] = mm_string[i];
                i++;
                l++;
                assert(l < 10);
            }
            skip_count_str[l] = '\0';
            ASSERT_MSG(l > 0, "Invalid skip count:%d.\n", l);
            sscanf(skip_count_str, "%d", &skip_counts[k]);
            ASSERT_MSG(skip_counts[k] >= 0, "Skip count cannot be negative: %d.\n", skip_counts[k]);
            k++;
        }
        skip_counts_len = k;
        i++;

        // for a reverse read the MM counts run along the other strand, so look for the complement
        char mb = rev ? base_complement_lookup[(int)modbase] : modbase;
        int idx = base_idx_lookup[(int)mb];

        int ml_idx = ml_start_idx;
        int base_rank = -1;

        // the bases the MM tag actually lists
        for (int c = 0; c < skip_counts_len; c++) {
            base_rank += skip_counts[c] + 1;

            int read_pos;
            if (modbase == 'N') {
                read_pos = rev ? seq_len - base_rank - 1 : base_rank;
            } else {
                read_pos = rev ? bases_pos[idx][bases_pos_lens[idx] - base_rank - 1] : bases_pos[idx][base_rank];
            }
            ASSERT_MSG(read_pos >= 0 && read_pos < seq_len, "Read pos cannot exceed seq len. read_pos: %d seq_len: %d\n", read_pos, seq_len);

            char read_base = seq_nt16_str[bam_seqi(seq, read_pos)];

            // did this base land on a CpG of an insertion, and does it match the ALT there?
            if (site_off[read_pos] < 0 || read_base != site_base[read_pos]) {
                if (mod_codes_len > 0) ml_idx = ml_start_idx + c * mod_codes_len + mod_codes_len - 1;
                continue;
            }

            for (int m = 0; m < mod_codes_len; m++) {
                ml_idx = ml_start_idx + c * mod_codes_len + m;

                char * mod_code = has_nums ? mod_codes : &(mod_codes[m]);
                modcodem_t * req_mod = get_required_mod(core, mod_code);
                if (req_mod == NULL) continue;

                ASSERT_MSG((uint32_t)ml_idx < ml_len, "read_id:%s mod prob index mismatch. ml_idx:%d ml_len:%d \n", bam_get_qname(record), ml_idx, ml_len);
                uint8_t mod_prob = ml[ml_idx];

                count_one(db, bam_i, strand, read_pos, site_var, site_off, mod_code, req_mod, mod_prob);
            }
        }
        if (skip_counts_len > 0) ml_start_idx = ml_idx + 1;

        // with a '.' flag every base the MM tag did not list is called unmodified
        if (status_flag == '.') {
            int listed_rank = -1;
            int prev_listed_rank = -1;
            int n_bases = (modbase == 'N') ? seq_len : bases_pos_lens[idx];

            for (int c = 0; c <= skip_counts_len; c++) {

                int upto;
                if (c < skip_counts_len) {
                    listed_rank += skip_counts[c] + 1;
                    upto = listed_rank;
                } else {
                    upto = n_bases; // everything after the last listed base
                }

                for (int s = prev_listed_rank + 1; s < upto; s++) {
                    int read_pos;
                    if (modbase == 'N') {
                        read_pos = rev ? seq_len - s - 1 : s;
                    } else {
                        read_pos = rev ? bases_pos[idx][bases_pos_lens[idx] - s - 1] : bases_pos[idx][s];
                    }
                    ASSERT_MSG(read_pos >= 0 && read_pos < seq_len, "Read pos cannot exceed seq len. read_pos: %d seq_len: %d\n", read_pos, seq_len);

                    char read_base = seq_nt16_str[bam_seqi(seq, read_pos)];

                    if (site_off[read_pos] < 0 || read_base != site_base[read_pos]) continue;

                    for (int m = 0; m < mod_codes_len; m++) {
                        char * mod_code = has_nums ? mod_codes : &(mod_codes[m]);
                        modcodem_t * req_mod = get_required_mod(core, mod_code);
                        if (req_mod == NULL) continue;

                        int var_idx = site_var[read_pos];
                        update_var_map(db->freq_maps[bam_i], variants[var_idx].chrom, variants[var_idx].pos,
                                       var_idx, site_off[read_pos], strand, mod_code, 1, 0);
                    }
                }

                prev_listed_rank = listed_rank;
            }
        }
    }
}


void varfreq_single(core_t * core, db_t * db, int32_t bam_i) {

    bam1_t * record = db->bam_recs[bam_i];
    const char * qname = bam_get_qname(record);
    bam_hdr_t * hdr = core->bam_hdr;
    int32_t tid = record->core.tid;
    const char * tname = (tid >= 0) ? hdr->target_name[tid] : "*";
    int seq_len = record->core.l_qseq;

    // 1. the variants this read was listed as supporting
    khiter_t k = kh_get(readm, read_map, qname);
    if (k == kh_end(read_map)) return;
    var_list_t * vl = kh_value(read_map, k);

    ref_t * ref = get_ref(tname);
    if (ref == NULL) { // the read is aligned somewhere the reference does not cover
        pthread_mutex_lock(&var_mutex);
        n_read_wrong_contig += vl->n;
        pthread_mutex_unlock(&var_mutex);
        return;
    }
    ASSERT_MSG(ref->ref_seq_length == (int)hdr->target_len[tid],
                "Contig %s is %d bases in the reference but %d in the bam header. Is this the reference the bam was aligned to?\n",
                tname, ref->ref_seq_length, (int)hdr->target_len[tid]);


    int hp = get_hp_tag(record);

    // 2. marks on the read: which variant, which offset and which base we expect
    // TO-DO: allocate these in laod_db
    int * site_var = (int *)malloc(sizeof(int) * seq_len);
    MALLOC_CHK(site_var);
    int * site_off = (int *)malloc(sizeof(int) * seq_len);
    MALLOC_CHK(site_off);
    char * site_base = (char *)malloc(sizeof(char) * seq_len);
    MALLOC_CHK(site_base);

    for (int i = 0; i < seq_len; i++) {
        site_var[i] = -1;
        site_off[i] = -1;
        site_base[i] = 0;
    }

    // 3. place the read on each variant it supports
    int count = 0;
    for (int v = 0; v < vl->n; v++) {
        int var_idx = vl->idx[v];
        var_t * var = &variants[var_idx];

        if (strcmp(var->chrom, tname) != 0) { // read contig doesn't match variant contig
            pthread_mutex_lock(&var_mutex);
            n_read_wrong_contig++;
            pthread_mutex_unlock(&var_mutex);
            continue;
        }

        if (hp == 0) { // no haplotype tag
            pthread_mutex_lock(&var_mutex);
            n_read_no_hp++;
            pthread_mutex_unlock(&var_mutex);
            continue;
        }

        if (var->alt_hap != 0 && hp != var->alt_hap) { // read haplotype tag does not match variant haplotype
            pthread_mutex_lock(&var_mutex);
            n_read_hp_mismatch++;
            pthread_mutex_unlock(&var_mutex);
            continue;
        }

        if (mark_sites_for_variant(core, record, ref, var_idx, site_var, site_off, site_base)) {
            count = 1;
        }
    }

    // 4. one pass over the MM/ML tags for everything we marked
    if (count) {
        count_calls(core, db, bam_i, site_var, site_off, site_base);
    }

    free(site_base);
    free(site_off);
    free(site_var);
}

/* where the CpG sits relative to the insertion. */
static const char * site_region(var_t * var, int offset, char strand) {
    int o = (strand == '+') ? offset : offset - 1;
    if (o == 0) return "JUNCTION_5P";
    if (o == var->ins_len) return "JUNCTION_3P";
    return "INTERNAL";
}

static const char * hap_string(var_t * var) {
    if (var->alt_hap == 1) return "1";
    if (var->alt_hap == 2) return "2";
    return "1,2";
}

void print_varfreq_header(core_t * core) {
    if (!core->opt.bedmethyl_out) {
        fprintf(core->opt.output_fp, "contig\tstart\tend\tstrand\tn_called\tn_mod\tfreq\tmod_code\tvar_id\tins_offset\tregion\talt_hap\tsvlen\tn_carrier\n");
    }
}

void print_varfreq_output(core_t * core) {

    khash_t(freqm) * freq_map = core->freq_map;
    khint_t map_size = kh_size(freq_map);
    if (map_size == 0) return;

    double sort_start = realtime();

    freq_kv_t * sorted_arr = (freq_kv_t *)malloc(sizeof(freq_kv_t) * map_size);
    MALLOC_CHK(sorted_arr);

    int size = 0;
    for (khint_t k = kh_begin(freq_map); k != kh_end(freq_map); k++) {
        if (kh_exist(freq_map, k)) {
            sorted_arr[size].key = (char *)kh_key(freq_map, k);
            sorted_arr[size].freq = kh_value(freq_map, k);
            size++;
        }
    }
    ks_introsort_varfreq(size, sorted_arr);
    core->sort_time = realtime() - sort_start;

    double output_start = realtime();
    FILE * out_fp = core->opt.output_fp;

    for (int i = 0; i < size; i++) {
        freq_t * freq = sorted_arr[i].freq;

        char * contig = NULL;
        char * mod_code = NULL;
        int pos, var_idx, offset;
        char strand;
        decode_var_key(sorted_arr[i].key, &contig, &pos, &var_idx, &offset, &strand, &mod_code);

        var_t * var = &variants[var_idx];
        int end = pos + 1;

        if (core->opt.bedmethyl_out) {
            double freq_value = (double)freq->n_mod * 100 / freq->n_called;
            fprintf(out_fp, "%s\t%d\t%d\t%s\t%d\t%c\t%d\t%d\t255,0,0\t%d\t%f\t%s\t%d\t%s\t%s\t%d\t%d\n",
                    contig, pos, end, mod_code, freq->n_called, strand, pos, end, freq->n_called, freq_value,
                    var->id, offset, site_region(var, offset, strand), hap_string(var), var->ins_len, var->n_carrier);
        } else {
            double freq_value = (double)freq->n_mod / freq->n_called;
            fprintf(out_fp, "%s\t%d\t%d\t%c\t%d\t%d\t%f\t%s\t%s\t%d\t%s\t%s\t%d\t%d\n",
                    contig, pos, end, strand, freq->n_called, freq->n_mod, freq_value, mod_code,
                    var->id, offset, site_region(var, offset, strand), hap_string(var), var->ins_len, var->n_carrier);
        }

        free(contig);
        free(mod_code);
    }

    if (out_fp != stdout) fclose(out_fp);

    free(sorted_arr);
    core->output_time += (realtime() - output_start);
}

void print_varfreq_log() {
    fprintf(stderr, "[%s] VCF records: %ld total, %ld used\n", __func__, (long)n_vcf_total, (long)n_variants);
    fprintf(stderr, "[%s] VCF skipped: %ld not PASS, %ld not INS, %ld unusable ALT, %ld unusable GT, %ld no RNAMES, %ld contig not in reference, %ld disagree with the reference\n", __func__,
            (long)n_vcf_not_pass, (long)n_vcf_not_ins, (long)n_vcf_bad_alt, (long)n_vcf_bad_gt, (long)n_vcf_no_rnames,
            (long)n_vcf_no_contig, (long)n_vcf_ref_mismatch);
    fprintf(stderr, "[%s] read-variant pairs used: %ld\n", __func__, (long)n_read_used);
    if (n_site_clash > 0) {
        WARNING("%ld read positions were claimed by more than one insertion and only the last one counted. Try a smaller -w.", (long)n_site_clash);
    }
    fprintf(stderr, "[%s] read-variant pairs skipped: %ld no HP tag, %ld HP does not match the GT, %ld flank too short, %ld aligned to another contig\n", __func__,
            (long)n_read_no_hp, (long)n_read_hp_mismatch, (long)n_read_short_flank, (long)n_read_wrong_contig);
}
