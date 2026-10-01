#!/bin/bash
# Tests for "minimod varfreq".
# Everything it needs is in test/data/varfreq, nothing is downloaded.
# Usage: test/varfreq_test.sh [mem]

BLUE='\033[0;34m'
RED='\033[0;31m'
GREEN='\033[0;32m'
NC='\033[0m'

die() {
	echo -e "${RED}$1${NC}" >&2
	echo
	exit 1
}

if [ "$1" = 'mem' ]; then
    mem=1
else
    mem=0
fi

ex() {
    if [ $mem -eq 1 ]; then
        valgrind --leak-check=full --error-exitcode=1 "$@"
    else
        "$@"
    fi
}

DATA=test/data/varfreq
EXP=test/expected/varfreq
TMP=test/tmp/varfreq

REF=$DATA/hg03902_chm13_chr22.fa.gz
BAM=$DATA/hg03902_chm13_chr22_ins.bam
VCF=$DATA/hg03902_chm13_chr22_ins.vcf

mkdir -p $TMP || die "Creating the tmp directory failed"

testname="varfreq chr22 ins"
echo -e "${BLUE}${testname}${NC}"
ex ./minimod varfreq -t 8 $REF $BAM $VCF > $TMP/chr22_ins.m.varfreq.tsv || die "${testname} failed"
diff -q $TMP/chr22_ins.m.varfreq.tsv $EXP/chr22_ins.m.varfreq.tsv || die "${testname} failed: output does not match expected output"

testname="varfreq -b chr22 ins"
echo -e "${BLUE}${testname}${NC}"
ex ./minimod varfreq -b -t 8 $REF $BAM $VCF > $TMP/chr22_ins.m.varfreq.bedmethyl || die "${testname} failed"
diff -q $TMP/chr22_ins.m.varfreq.bedmethyl $EXP/chr22_ins.m.varfreq.bedmethyl || die "${testname} failed: output does not match expected output"

testname="varfreq -b -c m,h chr22 ins"
echo -e "${BLUE}${testname}${NC}"
ex ./minimod varfreq -b -t 8 -c "m,h" $REF $BAM $VCF > $TMP/chr22_ins.mh.varfreq.bedmethyl || die "${testname} failed"
diff -q $TMP/chr22_ins.mh.varfreq.bedmethyl $EXP/chr22_ins.mh.varfreq.bedmethyl || die "${testname} failed: output does not match expected output"

# the output must not depend on how the work is split across threads and batches
testname="varfreq single thread matches 8 threads"
echo -e "${BLUE}${testname}${NC}"
ex ./minimod varfreq -b -t 1 $REF $BAM $VCF > $TMP/chr22_ins.t1.bedmethyl || die "${testname} failed"
diff -q $TMP/chr22_ins.t1.bedmethyl $EXP/chr22_ins.m.varfreq.bedmethyl || die "${testname} failed: output depends on the thread count"

testname="varfreq small batches match large batches"
echo -e "${BLUE}${testname}${NC}"
ex ./minimod varfreq -b -t 8 -K 16 -B 1M $REF $BAM $VCF > $TMP/chr22_ins.k16.bedmethyl || die "${testname} failed"
diff -q $TMP/chr22_ins.k16.bedmethyl $EXP/chr22_ins.m.varfreq.bedmethyl || die "${testname} failed: output depends on the batch size"

# the band only limits the alignment, it must not change the result
testname="varfreq wide band matches default band"
echo -e "${BLUE}${testname}${NC}"
ex ./minimod varfreq -b -t 8 --band 1000 $REF $BAM $VCF > $TMP/chr22_ins.band1000.bedmethyl || die "${testname} failed"
diff -q $TMP/chr22_ins.band1000.bedmethyl $EXP/chr22_ins.m.varfreq.bedmethyl || die "${testname} failed: the band width changed the output"

# a vcf holding no usable insertion is not an error, it just produces nothing
testname="varfreq with no insertions in the vcf"
echo -e "${BLUE}${testname}${NC}"
grep -E "^#|SVTYPE=DEL" $VCF > $TMP/del_only.vcf || die "${testname} failed to build the vcf"
ex ./minimod varfreq -b -t 8 $REF $BAM $TMP/del_only.vcf > $TMP/del_only.bedmethyl || die "${testname} failed"
[ -s $TMP/del_only.bedmethyl ] && die "${testname} failed: expected no output"

# every ALT here is one varfreq cannot derive a sequence from, so all 8 must be skipped
# and nothing may be emitted. guards against a symbolic or breakend ALT being fed into
# the aligner as if it were insertable sequence.
testname="varfreq skips ALTs it cannot derive"
echo -e "${BLUE}${testname}${NC}"
ex ./minimod varfreq -b -t 8 $REF $BAM $DATA/pathological_alts.vcf > $TMP/pathological.bedmethyl 2> $TMP/pathological.log || die "${testname} failed"
[ -s $TMP/pathological.bedmethyl ] && die "${testname} failed: expected no output"
grep -q "8 unusable ALT" $TMP/pathological.log || die "${testname} failed: expected all 8 records to be counted as unusable ALT"

# load_variants needs exactly one sample, otherwise bcf_get_genotypes returns the
# genotypes of several samples and get_alt_hap would read two of them as one diploid call
testname="varfreq rejects a multi-sample vcf"
echo -e "${BLUE}${testname}${NC}"
tr -d '\r' < $VCF | awk -F'\t' 'BEGIN{OFS="\t"} /^##/{print; next} /^#CHROM/{print $0, "SAMPLE2"; next} {print $0, $10}' > $TMP/two_samples.vcf || die "${testname} failed to build the vcf"
./minimod varfreq -b $REF $BAM $TMP/two_samples.vcf > /dev/null 2>&1 && die "${testname} failed: expected a non-zero exit"

testname="varfreq rejects a missing vcf"
echo -e "${BLUE}${testname}${NC}"
./minimod varfreq -b $REF $BAM $TMP/does_not_exist.vcf > /dev/null 2>&1 && die "${testname} failed: expected a non-zero exit"

testname="varfreq rejects min-flank larger than flank"
echo -e "${BLUE}${testname}${NC}"
./minimod varfreq -b -w 100 --min-flank 200 $REF $BAM $VCF > /dev/null 2>&1 && die "${testname} failed: expected a non-zero exit"

echo -e "${GREEN}ALL varfreq TESTS PASSED !${NC}"
