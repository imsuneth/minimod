#!/usr/bin/env python3

# Usage: test/varfreq_summary.py variants.vcf freq.bedmethyl varfreq.bedmethyl [reference.fa]
#
# Counts what went into a varfreq run and what came out of it.
# The optional reference.fa enables the reference-CpG check, which asks whether a
# JUNCTION_5P CpG was actually gained or was already there in the reference.
#
# varfreq bedmethyl columns:
#   0 contig  1 start  2 end  3 mod_code  4 n_called  5 strand  6 thickStart
#   7 thickEnd  8 colour  9 coverage  10 freq  11 var_id  12 ins_offset
#   13 region  14 alt_hap  15 svlen  16 n_carrier

import sys
import os
import gzip

if len(sys.argv) not in (4, 5):
    print("Usage: {} variants.vcf freq.bedmethyl varfreq.bedmethyl [reference.fa]".format(sys.argv[0]))
    sys.exit(1)

vcf_file = sys.argv[1]
freq_file = sys.argv[2]
varfreq_file = sys.argv[3]
reference_file = sys.argv[4] if len(sys.argv) == 5 else None


def opener(fn):
    return gzip.open(fn, "rt") if fn.endswith(".gz") else open(fn, "rt")


def counts(title, d, sort_key=None):
    print("## {}:".format(title))
    keys = sorted(d, key=sort_key) if sort_key else sorted(d, key=lambda k: -d[k])
    total = 0
    for k in keys:
        print("{}: {}".format(k, d[k]))
        total += d[k]
    print("Total: {}".format(total))


def bump(d, k):
    d[k] = d.get(k, 0) + 1


# ---------------------------------------------------------------- vcf

def vcf_load(fn):
    # var_id -> (svtype, filter, gt, n_cpg_alt, svlen, usable)
    d = {}
    svtype_counts = {}
    filter_counts = {}
    gt_counts = {}

    with opener(fn) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            info = parts[7]

            svtype = "?"
            for kv in info.split(";"):
                if kv.startswith("SVTYPE="):
                    svtype = kv[len("SVTYPE="):]
                    break
            bump(svtype_counts, svtype)
            bump(filter_counts, parts[6])

            if svtype != "INS":
                continue

            gt = parts[9].split(":")[0]
            bump(gt_counts, gt)

            ref_allele = parts[3].upper()
            alt_allele = parts[4].upper()

            # the same rules varfreq itself applies
            usable = (parts[6] == "PASS"
                      and gt in ("1|0", "0|1", "1|1", "1/1")
                      and not alt_allele.startswith("<")
                      and len(ref_allele) == 1
                      and len(alt_allele) >= 3)

            d[parts[2]] = (svtype, parts[6], gt, alt_allele.count("CG"),
                           len(alt_allele) - 1, usable)

    print("# vcf: {}".format(fn))
    counts("Variant type summary", svtype_counts)
    counts("Filter summary", filter_counts)
    counts("Genotype summary (INS only)", gt_counts)
    print("## Usable by varfreq (PASS, phased or hom-alt, plain ALT):")
    print("{} of {} insertions".format(sum(1 for v in d.values() if v[5]), len(d)))
    return d


# ---------------------------------------------------------------- freq

def freq_summary(fn):
    mod_counts = {}
    hap_counts = {}
    n_rows = 0

    with opener(fn) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            n_rows += 1
            bump(mod_counts, parts[3])
            bump(hap_counts, parts[11] if len(parts) > 11 else "(no haplotype column)")

    print("# freq : {}".format(fn))
    print("rows: {}".format(n_rows))
    counts("Modification summary", mod_counts)
    counts("Haplotype summary", hap_counts)


# ---------------------------------------------------------------- varfreq

def varfreq_load(fn):
    # var_id -> set of CpGs, named by the offset of the C
    cpgs = {}
    mod_counts = {}
    hap_counts = {}
    region_counts = {}
    n_rows = 0
    n_called_total = 0
    svlen = {}
    carrier = {}

    with opener(fn) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            n_rows += 1
            n_called_total += int(parts[4])

            bump(mod_counts, parts[3])
            bump(hap_counts, parts[14])
            bump(region_counts, parts[13])

            var_id = parts[11]
            offset = int(parts[12])
            # on the minus strand the C of the CpG is one base earlier
            o = offset if parts[5] == "+" else offset - 1
            cpgs.setdefault(var_id, set()).add(o)
            svlen[var_id] = int(parts[15])
            carrier[var_id] = int(parts[16])

    print("# varfreq : {}".format(fn))
    print("rows: {}".format(n_rows))
    print("insertions with output: {}".format(len(cpgs)))
    print("distinct CpG sites: {}".format(sum(len(s) for s in cpgs.values())))
    print("mean n_called per row: {:.2f}".format(n_called_total / n_rows if n_rows else 0))
    counts("Modification summary", mod_counts)
    counts("ALT haplotype summary", hap_counts)
    counts("Region summary", region_counts)
    return cpgs, svlen, carrier


def cpg_recovery(cpgs, variants):
    # how many of the CG dinucleotides in each ALT did varfreq actually report?
    # a JUNCTION_3P CpG has its G in the reference, so it is not one of them
    full = partial = 0
    reported = expected = 0
    missing_examples = []

    for var_id, sites in cpgs.items():
        v = variants.get(var_id)
        if v is None:
            continue
        n_alt = v[3]
        n_seen = sum(1 for o in sites if o < v[4])  # o == svlen is the 3' junction
        reported += n_seen
        expected += n_alt
        if n_alt and n_seen >= n_alt:
            full += 1
        else:
            partial += 1
            if len(missing_examples) < 5 and n_alt:
                missing_examples.append((var_id, n_seen, n_alt))

    print("## CpG recovery (reported vs CG dinucleotides in the ALT):")
    print("insertions where every ALT CpG was reported: {}".format(full))
    print("insertions with at least one ALT CpG unreported: {}".format(partial))
    if expected:
        print("overall: {} of {} ({:.1f}%)".format(reported, expected, 100.0 * reported / expected))
    for var_id, n_seen, n_alt in missing_examples:
        print("  e.g. {} reported {} of {}".format(var_id, n_seen, n_alt))


# ---------------------------------------------------------------- reference check

class FaidxReference:
    # random access via the .fai, so a whole genome does not have to be held in memory
    def __init__(self, fn):
        self.fh = open(fn, "rb")
        self.idx = {}
        for line in open(fn + ".fai"):
            f = line.split("\t")
            self.idx[f[0]] = (int(f[1]), int(f[2]), int(f[3]), int(f[4]))

    def two_bases(self, contig, pos):
        # the bases at pos and pos+1, or None if we cannot reach them
        entry = self.idx.get(contig)
        if entry is None:
            return None
        length, off, lb, lw = entry
        if pos < 0 or pos + 1 >= length:
            return None
        out = []
        for p in (pos, pos + 1):
            self.fh.seek(off + (p // lb) * lw + (p % lb))
            out.append(self.fh.read(1).decode())
        return "".join(out).upper()


class WholeReference:
    # fallback for a compressed or unindexed fasta
    def __init__(self, seqs):
        self.seqs = seqs

    def two_bases(self, contig, pos):
        seq = self.seqs.get(contig)
        if seq is None or pos < 0 or pos + 1 >= len(seq):
            return None
        return seq[pos:pos + 2].upper()


def load_reference(fn):
    if not fn.endswith(".gz") and os.path.exists(fn + ".fai"):
        return FaidxReference(fn)

    seqs = {}
    name = None
    chunks = []
    with opener(fn) as f:
        for line in f:
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(chunks)
                name = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line.strip())
    if name is not None:
        seqs[name] = "".join(chunks)
    return WholeReference(seqs)


def reference_cpg_check(fn, reference):
    # a JUNCTION_5P CpG has its C on the reference anchor base and its G on the first
    # inserted base. if the reference already reads CG there then the insertion did not
    # create that CpG, it was there all along and freq reports it too.
    checked = 0
    already = []

    with opener(fn) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if parts[13] != "JUNCTION_5P" or parts[5] != "+":
                continue

            contig = parts[0]
            pos = int(parts[1])
            bases = reference.two_bases(contig, pos)
            if bases is None:
                continue

            checked += 1
            if bases == "CG":
                already.append((contig, pos, parts[11]))

    print("## Reference CpG check:")
    print("JUNCTION_5P '+' sites checked: {}".format(checked))
    print("  already a CpG in the reference (not gained): {}".format(len(already)))
    for contig, pos, var_id in already[:10]:
        print("  {}:{} {}".format(contig, pos, var_id))


# ---------------------------------------------------------------- main

variants = vcf_load(vcf_file)
print()
freq_summary(freq_file)
print()
cpgs, svlen, carrier = varfreq_load(varfreq_file)
print()
cpg_recovery(cpgs, variants)
print()

if reference_file is None:
    print("## Reference CpG check:")
    print("(skipped: pass reference.fa as the 4th argument to run this check)")
else:
    reference_cpg_check(varfreq_file, load_reference(reference_file))
