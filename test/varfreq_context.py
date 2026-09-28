#!/usr/bin/env python3

# Usage: test/varfreq_context.py freq.bedmethyl varfreq.bedmethyl min_window [--bed out.bed]
#
# Compares the methylation of the CpGs inside each insertion against the reference
# CpGs on either side of it, taken from a plain freq run over the same reads.
#
# varfreq bedmethyl columns:
#   0 contig  1 start  2 end  3 mod_code  4 n_called  5 strand  6 thickStart
#   7 thickEnd  8 colour  9 coverage  10 freq  11 var_id  12 ins_offset
#   13 region  14 alt_hap  15 svlen  16 n_carrier
#
# freq bedmethyl is the usual 11 columns, plus a haplotype column when the run used
# --haplotypes. In that case freq writes a row per haplotype (1, 2, and 0 for reads
# with no HP tag) plus a '*' row holding their total, so a het insertion is compared
# against its own haplotype and a hom-alt one against the total.

import sys
import gzip
from bisect import bisect_left, bisect_right
from statistics import mean, stdev

argv = sys.argv[1:]
bed_file = None
if "--bed" in argv:
    i = argv.index("--bed")
    if i + 1 >= len(argv):
        print("--bed requires an output file path")
        sys.exit(1)
    bed_file = argv[i + 1]
    del argv[i:i + 2]

if len(argv) != 3:
    print("Usage: {} freq.bedmethyl varfreq.bedmethyl min_window [--bed out.bed]".format(sys.argv[0]))
    sys.exit(1)

freq_file = argv[0]
varfreq_file = argv[1]
min_window = int(argv[2])

CATEGORIES = ("meth_in_unmeth", "unmeth_in_meth", "equal")
COLOUR = {"meth_in_unmeth": "255,0,0", "unmeth_in_meth": "0,0,255", "equal": "128,128,128"}


def opener(fn):
    return gzip.open(fn, "rt") if fn.endswith(".gz") else open(fn, "rt")


def site_pct(strand_pcts):
    # a CpG is one site, so average whatever strands were seen for it
    return mean(strand_pcts)


def stats(pcts):
    if not pcts:
        return (0, None, None)
    return (len(pcts), mean(pcts), stdev(pcts) if len(pcts) > 1 else None)


def fmt(x):
    return "NA" if x is None else "{:.1f}".format(x)


def merge_windows(wins):
    # per contig: sorted, non-overlapping (lo, hi) we need background from
    out = {}
    for contig, spans in wins.items():
        spans.sort()
        merged = []
        for lo, hi in spans:
            if merged and lo <= merged[-1][1]:
                merged[-1][1] = max(merged[-1][1], hi)
            else:
                merged.append([lo, hi])
        out[contig] = ([m[0] for m in merged], [m[1] for m in merged])
    return out


def in_window(windows, contig, pos):
    entry = windows.get(contig)
    if entry is None:
        return False
    los, his = entry
    i = bisect_right(los, pos) - 1
    return i >= 0 and pos < his[i]


def freq_load(fn, windows):
    # (contig, hap, mod_code) -> (sorted site positions, matching percentages)
    per_site = {}
    has_hap = False

    with opener(fn) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            contig, mod_code, strand = parts[0], parts[3], parts[5]
            pos, pct = int(parts[1]), float(parts[10])

            if not in_window(windows, contig, pos):
                continue

            if len(parts) > 11:
                has_hap = True
                hap = parts[11]
            else:
                hap = "*"

            # name a CpG by the position of its C on the forward strand
            site_pos = pos if strand == "+" else pos - 1
            per_site.setdefault((contig, hap, mod_code), {}).setdefault(site_pos, []).append(pct)

    out = {}
    for key, d in per_site.items():
        sites = sorted(d)
        out[key] = (sites, [site_pct(d[s]) for s in sites])
    return out, has_hap


def varfreq_load(fn):
    # (contig, var_pos, var_id, alt_hap, mod_code) -> (svlen, n_carrier, {cpg: [pcts]})
    groups = {}
    with opener(fn) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            contig, mod_code, strand = parts[0], parts[3], parts[5]
            pos, pct = int(parts[1]), float(parts[10])
            var_id, region, alt_hap = parts[11], parts[13], parts[14]
            offset, svlen, n_carrier = int(parts[12]), int(parts[15]), int(parts[16])

            # on the minus strand the C of the CpG is one base earlier
            cpg = offset if strand == "+" else offset - 1

            # the 5' junction CpG has its C on a reference base, so freq reports it too.
            # leave it out, otherwise the same site sits on both sides of the comparison.
            if region == "JUNCTION_5P":
                continue

            key = (contig, pos, var_id, alt_hap, mod_code)
            g = groups.get(key)
            if g is None:
                g = [svlen, n_carrier, {}]
                groups[key] = g
            g[2].setdefault(cpg, []).append(pct)
    return groups


def background(freqs, contig, haps, mod_code, var_pos, var_end, radius):
    # reference CpGs within radius on each side, excluding the anchor itself
    left, right = [], []
    for hap in haps:
        entry = freqs.get((contig, hap, mod_code))
        if entry is None:
            continue
        sites, pcts = entry
        lo = bisect_left(sites, var_pos - radius)
        mid_l = bisect_left(sites, var_pos)
        mid_r = bisect_left(sites, var_end)
        hi = bisect_left(sites, var_end + radius)
        left.extend(pcts[lo:mid_l])
        right.extend(pcts[mid_r:hi])
    return left, right


varfreqs = varfreq_load(varfreq_file)

# work out the background windows before reading freq, so a whole-genome freq file
# can be streamed instead of held in memory
wanted = {}
for (contig, var_pos, _var_id, _hap, _mod), (svlen, _carrier, _sites) in varfreqs.items():
    radius = int(max(min_window, svlen) / 2) + 1
    wanted.setdefault(contig, []).append([var_pos - radius, var_pos + 1 + radius])
windows = merge_windows(wanted)

freqs, freq_has_hap = freq_load(freq_file, windows)

print("# freq    : {}".format(freq_file), file=sys.stderr)
print("# varfreq : {}".format(varfreq_file), file=sys.stderr)
print("# min_window={} bp  per-CpG freq = mean of its + and - strand freqs  insertion level = unweighted mean over its CpGs".format(
    min_window), file=sys.stderr)
print("# background = reference CpGs within max(min_window, svlen)/2 bp either side of the anchor", file=sys.stderr)
if freq_has_hap:
    print("# freq has a haplotype column, so each insertion is compared against its own haplotype", file=sys.stderr)
else:
    print("# freq has no haplotype column, so the background is not haplotype matched", file=sys.stderr)
print("# meth_in_unmeth if the insertion is more methylated than its background, unmeth_in_meth if less", file=sys.stderr)
print("# meth_mean_diff = var_meth_mean - bg_meth_mean, so its sign follows the category", file=sys.stderr)

buckets = {c: [] for c in CATEGORIES}
skipped_no_bg = 0
skipped_no_cpg = 0

for key in sorted(varfreqs):
    contig, var_pos, var_id, alt_hap, mod_code = key
    svlen, n_carrier, var_sites = varfreqs[key]

    if not var_sites:
        skipped_no_cpg += 1
        continue

    var_end = var_pos + 1  # the reference span of an insertion is its single anchor base
    var_ncpg, var_mean, var_sd = stats([site_pct(v) for v in var_sites.values()])

    radius = max(min_window, svlen) / 2

    # freq --haplotypes writes a row per haplotype (1, 2, and 0 for untagged reads)
    # plus a '*' row holding their total. A het insertion is compared against its own
    # haplotype; a hom-alt one sits on both, so use the '*' total rather than pooling
    # 1 and 2, which would count every background site twice.
    if not freq_has_hap or alt_hap == "1,2":
        haps = ["*"]
    else:
        haps = [alt_hap]

    bg_l_pcts, bg_r_pcts = background(freqs, contig, haps, mod_code, var_pos, var_end, radius)
    bg_l_ncpg, bg_l_mean, _ = stats(bg_l_pcts)
    bg_r_ncpg, bg_r_mean, _ = stats(bg_r_pcts)
    bg_ncpg, bg_mean, bg_sd = stats(bg_l_pcts + bg_r_pcts)

    if bg_mean is None:
        skipped_no_bg += 1
        continue

    meth_diff = var_mean - bg_mean
    if meth_diff > 0:
        category = "meth_in_unmeth"
    elif meth_diff < 0:
        category = "unmeth_in_meth"
    else:
        category = "equal"

    buckets[category].append(
        (contig, var_pos, var_end, var_id, alt_hap, mod_code, svlen, n_carrier,
         var_ncpg, var_mean, var_sd,
         bg_ncpg, bg_mean, bg_sd, bg_l_ncpg, bg_l_mean, bg_r_ncpg, bg_r_mean,
         meth_diff, radius))

COLUMNS = ["chrom", "start", "end", "var_id", "alt_hap", "mod", "category",
           "svlen", "n_carrier", "var_ncpg", "var_meth_mean", "var_meth_sd",
           "bg_span", "bg_ncpg", "bg_meth_mean", "bg_meth_sd",
           "bg_l_ncpg", "bg_l_meth_mean", "bg_r_ncpg", "bg_r_meth_mean", "meth_mean_diff"]

tsv_rows = []
for category in CATEGORIES:
    for (contig, var_pos, var_end, var_id, alt_hap, mod_code, svlen, n_carrier,
         var_ncpg, var_mean, var_sd,
         bg_ncpg, bg_mean, bg_sd, bg_l_ncpg, bg_l_mean, bg_r_ncpg, bg_r_mean,
         meth_diff, radius) in buckets[category]:
        tsv_rows.append((contig, var_pos, var_end, var_id, alt_hap, mod_code, category,
                         svlen, n_carrier, var_ncpg, fmt(var_mean), fmt(var_sd),
                         int(2 * radius), bg_ncpg, fmt(bg_mean), fmt(bg_sd),
                         bg_l_ncpg, fmt(bg_l_mean), bg_r_ncpg, fmt(bg_r_mean),
                         fmt(meth_diff)))

tsv_rows.sort(key=lambda r: (r[0], r[1], r[3], r[5]))

print("\t".join(COLUMNS))
for r in tsv_rows:
    print("\t".join(str(f) for f in r))

print("# meth_in_unmeth={} unmeth_in_meth={} equal={} no_bg_cpg={} no_cpg={}".format(
    len(buckets["meth_in_unmeth"]), len(buckets["unmeth_in_meth"]), len(buckets["equal"]),
    skipped_no_bg, skipped_no_cpg), file=sys.stderr)


def write_bed(fn, buckets):
    rows = []
    for category in CATEGORIES:
        for (contig, var_pos, var_end, var_id, alt_hap, mod_code, svlen, n_carrier,
             var_ncpg, var_mean, var_sd,
             bg_ncpg, bg_mean, bg_sd, bg_l_ncpg, bg_l_mean, bg_r_ncpg, bg_r_mean,
             meth_diff, radius) in buckets[category]:
            name = "{}_{}_ins{}_hap{}_{}_v{:.0f}/b{:.0f}".format(
                var_id, category, svlen, alt_hap, mod_code, var_mean, bg_mean)
            score = min(1000, int(round(var_mean * 10)))
            rows.append((contig, var_pos, var_end, name, score, COLOUR[category]))

    rows.sort(key=lambda r: (r[0], r[1], r[2]))
    with open(fn, "w") as out:
        out.write('track name="varfreq_context" description="insertion vs background methylation" itemRgb="On"\n')
        for contig, start, end, name, score, colour in rows:
            out.write("{}\t{}\t{}\t{}\t{}\t.\t{}\t{}\t{}\n".format(
                contig, start, end, name, score, start, end, colour))
    return len(rows)


if bed_file:
    n = write_bed(bed_file, buckets)
    print("# wrote {} rows to {}".format(n, bed_file), file=sys.stderr)
