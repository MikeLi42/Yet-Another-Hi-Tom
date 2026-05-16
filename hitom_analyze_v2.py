"""
Hi-TOM genotyping pipeline v2 — original filtering thresholds.

Changes from v1 to match the original Perl implementation (Liu et al. 2019):
  1. Barcode quality: each of the 4 barcode bases must have Phred > 30
  2. Read-level pre-filter: discard reads with >10 % Ns OR >50 % bases Phred ≤ 5
  3. MAPQ ≥ 30  (was ≥ 20)
  4. Per-base Q25 gate on extracted window removed (original has none)
  5. Singleton recovery removed (original drops unpaired mates)
  6. Artifact scrubbing removed (original keeps all Large_Indel reads)

Output format identical to v1 / original Hi-TOM web platform.
"""

import gzip
import itertools
import re
import subprocess
from collections import Counter, defaultdict
import difflib


# ---------------------------------------------------------------------------
# Barcode definitions
# ---------------------------------------------------------------------------
F_BARCODES = {
    'GCGT': 1, 'GTAG': 2, 'ACGC': 3, 'CTCG': 4, 'GCTC': 5,
    'AGTC': 6, 'CGAC': 7, 'GATG': 8, 'ATAC': 9, 'CACA': 10,
    'GTGC': 11, 'ACTA': 12
}
R_BARCODES = {
    'GCGT': 'A', 'GTAG': 'B', 'ACGC': 'C', 'CTCG': 'D',
    'GCTC': 'E', 'AGTC': 'F', 'CGAC': 'G', 'GATG': 'H'
}


def _build_barcode_decoder(barcodes, max_mismatch=1):
    decode = {}
    for kmer in itertools.product('ACGT', repeat=4):
        kmer = ''.join(kmer)
        best = None
        best_dist = 5
        for bc in barcodes:
            d = sum(a != b for a, b in zip(kmer, bc))
            if d < best_dist:
                best_dist = d
                best = bc
            elif d == best_dist:
                best = None
        if best is not None and best_dist <= max_mismatch:
            decode[kmer] = best
    return decode


_F_DECODE = _build_barcode_decoder(F_BARCODES)
_R_DECODE = _build_barcode_decoder(R_BARCODES)

PHRED30 = 30 + 33  # '?' in Illumina 1.8+


def decode_barcode(observed, decoder, mapping):
    exact = mapping.get(observed)
    if exact is not None:
        return exact
    corrected = decoder.get(observed)
    if corrected is not None:
        return mapping[corrected]
    return None


def parse_fastq_barcodes(r1_path, r2_path):
    """
    Scan R1 / R2 FASTQ pair and return dict[qname] -> sample_id (e.g. 'A01').

    v2 change: each of the 4 barcode bases (positions 5-8) must have
    Phred quality > 30.  Read pairs whose barcode quality is ≤ 30 on
    any base are discarded.
    """
    qmap = {}
    with gzip.open(r1_path, 'rt') as f1, gzip.open(r2_path, 'rt') as f2:
        while True:
            n1 = f1.readline()
            if not n1:
                break
            s1 = f1.readline().strip()
            q1 = f1.readline().strip()  # '+'
            q1 = f1.readline().strip()  # quality
            n2 = f2.readline().strip()
            s2 = f2.readline().strip()
            q2 = f2.readline().strip()  # '+'
            q2 = f2.readline().strip()  # quality

            if len(s1) < 8 or len(s2) < 8:
                continue

            if len(q1) < 8 or len(q2) < 8:
                continue

            # ---- v2: barcode quality check (Phred > 30 on each base) ----
            bc_q1 = q1[4:8]
            bc_q2 = q2[4:8]
            if any(ord(b) < PHRED30 for b in bc_q1):
                continue
            if any(ord(b) < PHRED30 for b in bc_q2):
                continue

            qname = n1.strip().split()[0].lstrip('@').rsplit('/', 1)[0]
            col = decode_barcode(s1[4:8], _F_DECODE, F_BARCODES)
            row = decode_barcode(s2[4:8], _R_DECODE, R_BARCODES)
            if col is not None and row is not None:
                qmap[qname] = f"{row}{col:02d}"
    return qmap


# ---------------------------------------------------------------------------
# Automatic target-window detection
# ---------------------------------------------------------------------------

def auto_detect_target_window(ref_seq, bam_path, min_window=90, max_window=110,
                               log_callback=print):
    wsl_bam = bam_path.replace('\\', '/').replace('D:/', '/mnt/d/')
    cmd = ['wsl', 'bash', '-c', f'samtools depth "{wsl_bam}"']
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True)

    depth = {}
    for line in proc.stdout:
        parts = line.rstrip('\n').split('\t')
        if len(parts) < 3:
            continue
        pos = int(parts[1])
        dp = int(parts[2])
        depth[pos] = dp
    proc.wait()

    if not depth:
        raise ValueError("No depth information from BAM (empty or invalid)")

    positions = sorted(depth)
    peak_pos = max(depth, key=depth.get)
    peak_dp = depth[peak_pos]
    shoulder = max(1, int(peak_dp * 0.05))

    left = peak_pos
    while left > 1 and depth.get(left - 1, 0) >= shoulder:
        left -= 1

    right = peak_pos
    max_pos = max(depth)
    while right < max_pos and depth.get(right + 1, 0) >= shoulder:
        right += 1

    window_len = right - left + 1
    if window_len < min_window:
        expand = (min_window - window_len) // 2
        left = max(1, left - expand)
        right = min(len(ref_seq), right + expand)
        window_len = right - left + 1
        if window_len < min_window and right < len(ref_seq):
            right += min_window - window_len
    elif window_len > max_window:
        shrink = (window_len - max_window) // 2
        left += shrink
        right -= shrink
        window_len = right - left + 1
        if window_len > max_window:
            right -= 1

    target_start = left
    target_end = right
    wt_window = ref_seq[target_start - 1:target_end]
    avg_dp = sum(depth.get(i, 0) for i in range(left, right + 1)) / max(1, right - left + 1)
    log_callback(f"Auto-detected target window: {target_start}-{target_end} "
                 f"({len(wt_window)} bp, peak depth {peak_dp}, avg depth {avg_dp:.0f})")
    return target_start, target_end, wt_window


# ---------------------------------------------------------------------------
# CIGAR-based target-window extraction  (v2: no per-base Q25 gate)
# ---------------------------------------------------------------------------

def extract_sequence(pos, cigar, seq, qual, target_start, target_end):
    """
    Extract the part of ``seq`` that aligns to reference coordinates
    ``[target_start, target_end]`` (1-based, inclusive).

    v2 change: no Phred-quality gate on the extracted window.  The original
    Hi-TOM only pre-filters reads (remove >10 % Ns or >50 % bases Q ≤ 5)
    and uses MAPQ ≥ 30; it does not check per-base quality in the window.
    """
    if not qual or qual == '*' or len(qual) != len(seq):
        return None

    ref_pos = pos
    seq_pos = 0
    seq_parts = []

    for match in re.finditer(r'(\d+)([MIDNSHP=X])', cigar):
        length = int(match.group(1))
        op = match.group(2)

        if op in 'M=X':
            block_end = ref_pos + length - 1
            ol_start = max(ref_pos, target_start)
            ol_end = min(block_end, target_end)
            if ol_start <= ol_end:
                s_start = seq_pos + (ol_start - ref_pos)
                s_end = seq_pos + (ol_end - ref_pos) + 1
                seq_parts.append(seq[s_start:s_end])
            ref_pos += length
            seq_pos += length

        elif op == 'I':
            if ref_pos >= target_start and ref_pos <= target_end + 1:
                seq_parts.append(seq[seq_pos:seq_pos + length])
            seq_pos += length

        elif op in 'DN':
            ref_pos += length

        elif op == 'S':
            seq_pos += length

    extracted = ''.join(seq_parts)
    if not extracted:
        return None
    return extracted


# ---------------------------------------------------------------------------
# Read-level pre-filter (original Hi-TOM criteria)
# ---------------------------------------------------------------------------

def original_read_filter(seq, qual):
    """
    Apply the original Hi-TOM pre-filtering criteria (Liu et al. 2019,
    Materials and Methods — Filtering reads and mapping reads):

        (i)   discard if >10 % of bases are 'N'
        (ii)  discard if >50 % of bases have Phred quality ≤ 5

    Returns True if the read SHOULD be kept, False if it should be discarded.
    """
    if not seq or not qual:
        return False
    if len(seq) != len(qual):
        return False

    n_count = seq.count('N')
    if n_count / len(seq) > 0.10:
        return False

    lowq = sum(1 for q in qual if ord(q) - 33 <= 5)
    if lowq / len(qual) > 0.50:
        return False

    return True


# ---------------------------------------------------------------------------
# Variant calling  (unchanged from v1)
# ---------------------------------------------------------------------------

def call_variant(seq, wt_seq, min_overlap=50):
    if seq == wt_seq or seq in wt_seq:
        return 'WT', '-'

    sm = difflib.SequenceMatcher(None, wt_seq, seq)
    ops = list(sm.get_opcodes())

    while ops and ops[0][0] in ('delete', 'insert'):
        ops.pop(0)
    while ops and ops[-1][0] in ('delete', 'insert'):
        ops.pop()

    ops = [op for op in ops if op[0] != 'equal']

    if not ops:
        return 'WT', '-'
    if len(ops) == 1:
        tag, i1, i2, j1, j2 = ops[0]
        if tag == 'insert':
            return f'{j2 - j1}I', seq[j1:j2]
        elif tag == 'delete':
            return f'{i2 - i1}D', wt_seq[i1:i2]
        elif tag == 'replace' and i2 - i1 == 1 and j2 - j1 == 1:
            return 'SNP', f'{wt_seq[i1:i2]}->{seq[j1:j2]}'

    longest = sm.find_longest_match(0, len(wt_seq), 0, len(seq))
    if longest.size >= max(5, len(seq) // 4):
        offset = longest.a - longest.b
        wt_start = max(0, offset)
        seq_start = max(0, -offset)
        overlap = min(len(wt_seq) - wt_start, len(seq) - seq_start)

        if overlap >= min_overlap:
            mismatches = 0
            mismatch_pos = None
            for k in range(overlap):
                if wt_seq[wt_start + k] != seq[seq_start + k]:
                    mismatches += 1
                    mismatch_pos = k
                    if mismatches > 1:
                        break
            if mismatches == 1:
                w = wt_seq[wt_start + mismatch_pos]
                s = seq[seq_start + mismatch_pos]
                return 'SNP', f'{w}->{s}'

    return 'Large_Indel', '----'


# ---------------------------------------------------------------------------
# Main per-dataset processing
# ---------------------------------------------------------------------------

def process_dataset(dataset_num, r1_path, r2_path, bam_path, ref_path,
                    out_seq, out_geno, abundance_threshold=0.05,
                    log_callback=print):
    # ---- Reference & target window ----------------------------------------
    with open(ref_path) as f:
        f.readline()
        ref_seq = f.readline().strip()

    target_start, target_end, wt_window = auto_detect_target_window(
        ref_seq, bam_path, log_callback=log_callback)

    # ---- Demultiplex FASTQ (with Q30 barcode check) ----------------------
    log_callback("Parsing FASTQ barcodes (Q30 gate)...")
    qname_to_sample = parse_fastq_barcodes(r1_path, r2_path)
    log_callback(f"Assigned {len(qname_to_sample)} read pairs to samples")
    valid_qnames = set(qname_to_sample.keys())

    # ---- Parse BAM --------------------------------------------------------
    log_callback("Parsing BAM...")
    sample_haplotypes = defaultdict(Counter)
    buffer = {}
    total_processed = 0

    wsl_bam = bam_path.replace('\\', '/').replace('D:/', '/mnt/d/')
    cmd = ['wsl', 'bash', '-c', f'samtools view {wsl_bam}']
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True)

    for line in proc.stdout:
        cols = line.rstrip('\n').split('\t')
        qname = cols[0]
        if qname not in valid_qnames:
            continue

        flag = int(cols[1])
        mapq = int(cols[4])
        # v2: MAPQ >= 30 (original Hi-TOM uses samtools view -q 30)
        if flag & 0x900 or mapq < 30:
            continue

        pos = int(cols[3])
        cigar = cols[5]
        seq = cols[9]
        qual = cols[10] if len(cols) > 10 else '*'
        if cigar == '*':
            continue

        # v2: original read-level pre-filter (>10 % Ns OR >50 % Q <=5)
        if not original_read_filter(seq, qual):
            continue

        extracted = extract_sequence(pos, cigar, seq, qual,
                                     target_start, target_end)
        if not extracted:
            continue

        is_r1 = (flag & 0x40) != 0
        is_r2 = (flag & 0x80) != 0
        if not is_r1 and not is_r2:
            continue

        if qname not in buffer:
            buffer[qname] = [None, None]
        idx_read = 0 if is_r1 else 1
        buffer[qname][idx_read] = extracted

        if buffer[qname][0] is not None and buffer[qname][1] is not None:
            sample_id = qname_to_sample[qname]
            left_seq, right_seq = buffer[qname]
            sample_haplotypes[sample_id][(left_seq, right_seq)] += 1
            total_processed += 1
            del buffer[qname]

    proc.wait()

    # v2: NO singleton recovery — original drops unpaired mates
    # Unpaired reads remaining in buffer are simply discarded.

    total_pairs = sum(sum(c.values()) for c in sample_haplotypes.values())
    log_callback(f"Processed {total_pairs} read pairs")

    # ---- Write outputs ----------------------------------------------------
    log_callback("Writing outputs...")
    with open(out_seq, 'w') as fs, open(out_geno, 'w') as fg:
        fg.write('\t' + '\t'.join(str(i) for i in range(1, 13)) + '\n')
        seq_header = ('Sort\tReads number\tRatio\tLeft variation type\t'
                      'Right variation type\tLeft variation\tRight variation\t'
                      'Left reads seq\tRight reads seq\n')
        fs.write(seq_header)

        for row_letter in 'ABCDEFGH':
            geno_row = [row_letter]
            for col in range(1, 13):
                sample_id = f"{row_letter}{col:02d}"
                hap_counter = sample_haplotypes.get(sample_id, Counter())
                total = sum(hap_counter.values())

                fs.write(f"{sample_id}\n")

                if total == 0:
                    geno_row.append("-")
                    continue

                # ---- abundance filter (v2: NO artifact scrubbing) --------
                filtered = []
                for pair, cnt in hap_counter.items():
                    if (cnt / total) < abundance_threshold:
                        continue
                    filtered.append((pair, cnt))
                filtered.sort(key=lambda x: -x[1])

                if not filtered:
                    geno_row.append("-")
                    continue

                for rank, ((left_seq, right_seq), cnt) in enumerate(filtered, 1):
                    ratio = cnt / total * 100
                    l_var, l_det = call_variant(left_seq, wt_window)
                    r_var, r_det = call_variant(right_seq, wt_window)
                    fs.write(f"{rank}\t{cnt}\t{ratio:.2f}%\t{l_var}\t{r_var}\t"
                             f"{l_det}\t{r_det}\t{left_seq}\t{right_seq}\n")

                # ---- Genotype classification --------------------------------
                variant_types = set()
                has_wt = False
                has_snp = False
                has_inframe = False
                for (l, r), cnt in filtered:
                    l_var, _ = call_variant(l, wt_window)
                    r_var, _ = call_variant(r, wt_window)
                    variant_types.add((l_var, r_var))
                    if l_var == 'WT' and r_var == 'WT':
                        has_wt = True
                    if l_var == 'SNP' or r_var == 'SNP':
                        has_snp = True
                    for var in (l_var, r_var):
                        if var.endswith('D') or var.endswith('I'):
                            num = int(var[:-1])
                            if num % 3 == 0:
                                has_inframe = True
                                break
                n_types = len(variant_types)

                if n_types >= 3:
                    genotype = 'chimeric'
                elif n_types == 1:
                    if has_wt:
                        genotype = 'AA'
                    else:
                        genotype = 'aa'
                        if has_inframe:
                            genotype += '*'
                        if has_snp:
                            genotype += '#'
                elif n_types == 2:
                    if has_wt:
                        genotype = 'Aa'
                    else:
                        genotype = 'aa'
                        if has_inframe:
                            genotype += '*'
                        if has_snp:
                            genotype += '#'
                else:
                    genotype = '-'

                geno_row.append(genotype)

            fg.write('\t'.join(geno_row) + '\n')

        fg.write("Notice: * In-frame mutation\n")
        fg.write("\t# SNP\n")
        fg.write("\t'- Missing data\n")

    log_callback(f"Done: {out_seq}, {out_geno}")


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main_cli():
    datasets = []
    for ds in datasets:
        process_dataset(*ds)
        print()


if __name__ == '__main__':
    main_cli()
