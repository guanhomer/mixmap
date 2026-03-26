#!/usr/bin/env python3
"""
per_read_per_cpg.py

Extract per-read, per-CpG methylation calls from a BAM file and collapse
overlapping observations within each inferred fragment.

Overview
--------
This script processes bisulfite sequencing alignments in two stages:

1. Stage 1: iterate through aligned reads and emit candidate per-read
   methylation calls at CpG positions.
2. Stage 2: externally sort those calls, then collapse duplicate
   observations from the same fragment at the same CpG site.

The output is a tab-separated table with one row per retained
fragment-by-CpG observation.

Assumptions
-----------
- Input alignments are stored in BAM format and are coordinate-valid enough
  for reference-based sequence lookup.
- The reference FASTA is indexed with samtools faidx / pysam (.fai required).
- Strand interpretation follows MethylDackel-like conventions.
- Chromosome filtering is restricted to canonical human chromosomes:
  1-22, X, Y, MT (with chr-prefixed names also accepted).

Example
-------
    python per_read_per_cpg.py ref.fa input.bam output.tsv

    python per_read_per_cpg.py ref.fa input.bam output.tsv \
        --min-phred 10 \
        --min-mapq 20 \
        --overlap-policy best

Output columns
--------------
fragment_id     Integer fragment index assigned after sorting/collapse
chrom           Chromosome name
ref_pos         1-based reference coordinate of the CpG-associated base
strand          Genomic strand of the methylation call: '+' or '-'
call            'm' for methylated, 'u' for unmethylated
fragment_span   Span of the inferred fragment in reference bases
"""

import argparse
import csv
import os
import sys
import tempfile
import subprocess
from typing import Optional, Tuple, List

import pysam


# ---------------------------------------------------------------------------
# Strand and CpG logic
# ---------------------------------------------------------------------------

def get_strand(read: pysam.AlignedSegment) -> int:
    """
    Infer the bisulfite strand class for a read.

    This follows MethylDackel-like logic and supports both:
    - aligners that provide XG tags (e.g. Bismark-like output)
    - directional aligners without XG tags (e.g. bwa-meth-like behavior)

    Returns
    -------
    int
        Strand class code:
            1 = OT
            2 = OB
            3 = CTOT
            4 = CTOB
            0 = undetermined
    """
    xg = None
    try:
        xg = read.get_tag("XG")
        if xg not in ("CT", "GA"):
            xg = None
    except KeyError:
        xg = None

    is_paired = read.is_paired
    is_reverse = read.is_reverse
    is_read1 = read.is_read1
    is_read2 = read.is_read2

    if xg is None:
        # Directional aligners such as bwa-meth.
        # Strand is inferred from read orientation and mate identity.
        if is_paired:
            if is_read1 and is_reverse:
                return 2
            elif is_read1:
                return 1
            elif is_read2 and is_reverse:
                return 1
            elif is_read2:
                return 2
            return 0
        else:
            return 2 if is_reverse else 1

    # Bismark-like handling when XG is available.
    if xg == "CT":
        if is_paired:
            if is_read1 and not is_reverse:
                return 1
            elif is_read1 and is_reverse:
                return 3
            elif is_read2 and not is_reverse:
                return 3
            elif is_read2 and is_reverse:
                return 1
            return 0
        else:
            return 3 if is_reverse else 1

    # xg == "GA"
    if is_paired:
        if is_read1 and not is_reverse:
            return 4
        elif is_read1 and is_reverse:
            return 2
        elif is_read2 and not is_reverse:
            return 2
        elif is_read2 and is_reverse:
            return 4
        return 0
    else:
        return 2 if is_reverse else 4


def is_cpg(ref_seq: str, idx: int) -> int:
    """
    Determine whether a reference position belongs to a CpG dinucleotide.

    This mirrors MethylDackel-style CpG interpretation:
      +1 if ref[idx] is C and ref[idx + 1] is G
      -1 if ref[idx] is G and ref[idx - 1] is C
       0 otherwise

    Parameters
    ----------
    ref_seq : str
        Reference sequence string for the local interval.
    idx : int
        Zero-based index into `ref_seq`.

    Returns
    -------
    int
        +1 for forward CpG cytosine, -1 for reverse CpG guanine, 0 otherwise.
    """
    if idx < 0 or idx >= len(ref_seq):
        return 0

    base = ref_seq[idx]
    if base in ("C", "c"):
        if idx + 1 >= len(ref_seq):
            return 0
        return 1 if ref_seq[idx + 1] in ("G", "g") else 0
    elif base in ("G", "g"):
        if idx == 0:
            return 0
        return -1 if ref_seq[idx - 1] in ("C", "c") else 0
    return 0


def informative_call_and_output_strand(
    read_base: str, direction: int, strand_class: int
) -> Optional[Tuple[str, str]]:
    """
    Convert an aligned base into a methylation call, when informative.

    Calling logic:
      - direction == +1 and strand class is odd:
          C = methylated, T = unmethylated, output strand '+'
      - direction == -1 and strand class is even:
          G = methylated, A = unmethylated, output strand '-'

    Parameters
    ----------
    read_base : str
        Base observed in the read at the aligned position.
    direction : int
        CpG orientation returned by `is_cpg()`.
    strand_class : int
        Strand class returned by `get_strand()`.

    Returns
    -------
    Optional[Tuple[str, str]]
        (genomic_strand, call), where:
            genomic_strand is '+' or '-'
            call is 'm' or 'u'

        Returns None if the observed base is not informative.
    """
    rb = read_base.upper()

    if direction == 1 and (strand_class & 1) == 1:
        if rb == "C":
            return ("+", "m")
        elif rb == "T":
            return ("+", "u")
    elif direction == -1 and (strand_class & 1) == 0:
        if rb == "G":
            return ("-", "m")
        elif rb == "A":
            return ("-", "u")

    return None


# ---------------------------------------------------------------------------
# Input helpers
# ---------------------------------------------------------------------------

def allowed_chrom(chrom: str) -> bool:
    """
    Restrict processing to canonical chromosomes.

    Accepted names:
    - 1-22, X, Y, MT
    - chr1-chr22, chrX, chrY, chrM, chrMT

    Parameters
    ----------
    chrom : str
        Reference name from the BAM header.

    Returns
    -------
    bool
        True if the chromosome is retained.
    """
    if chrom.startswith("chr"):
        c = chrom[3:]
    else:
        c = chrom

    if c == "M":
        c = "MT"

    return c in {str(i) for i in range(1, 23)} | {"X", "Y", "MT"}


def open_bam(path: str) -> pysam.AlignmentFile:
    """
    Open a BAM file, or stdin if '-' is provided.
    """
    if path == "-":
        return pysam.AlignmentFile("-", "rb")
    return pysam.AlignmentFile(path, "rb")


def shell_quote(s: str) -> str:
    """
    Single-quote a string for safe use in a shell command.
    """
    return "'" + s.replace("'", "'\"'\"'") + "'"


# ---------------------------------------------------------------------------
# Stage 1: per-read CpG extraction
# ---------------------------------------------------------------------------

def stage1_extract(
    bam_path: str,
    fasta_path: str,
    tmp_tsv: str,
    min_phred: int,
    min_mapq: int,
    require_flags: int,
    ignore_flags: int,
) -> None:
    """
    Scan reads and write candidate per-read per-CpG calls to a temporary TSV.

    The temporary output contains one row per informative aligned base that:
    - maps to a CpG position
    - passes base-quality and mapping-quality filters
    - passes BAM flag filters
    - yields an interpretable methylation state

    Output columns
    --------------
    fragment_name, chrom, ref_pos, strand, call, base_qual, read_start, read_end
    """
    bam = open_bam(bam_path)
    fasta = pysam.FastaFile(fasta_path)

    with open(tmp_tsv, "w", newline="") as out:
        writer = csv.writer(out, delimiter="\t", lineterminator="\n")
        writer.writerow(
            [
                "fragment_name",
                "chrom",
                "ref_pos",
                "strand",
                "call",
                "base_qual",
                "read_start",
                "read_end",
            ]
        )

        for read in bam.fetch(until_eof=True):
            # Basic read-level filters.
            if read.is_unmapped or read.reference_id < 0:
                continue
            if read.mapping_quality < min_mapq:
                continue
            if require_flags and (read.flag & require_flags) != require_flags:
                continue
            if ignore_flags and (read.flag & ignore_flags) != 0:
                continue
            if read.query_sequence is None or read.query_qualities is None:
                continue

            chrom = bam.get_reference_name(read.reference_id)
            if chrom is None or not allowed_chrom(chrom):
                continue

            strand_class = get_strand(read)
            if strand_class == 0:
                continue

            ref_start = read.reference_start
            ref_end = read.reference_end
            if ref_start is None or ref_end is None or ref_end <= ref_start:
                continue

            # Output interval is expressed in 1-based inclusive coordinates.
            read_start = ref_start + 1
            read_end = ref_end

            ref_seq = fasta.fetch(chrom, ref_start, ref_end).upper()
            qseq = read.query_sequence
            quals = read.query_qualities

            # Iterate only over aligned query/reference pairs.
            # `matches_only=True` excludes insertions, deletions, and soft clips.
            for qpos, rpos in read.get_aligned_pairs(matches_only=True):
                if qpos is None or rpos is None:
                    continue

                q = quals[qpos]
                if q < min_phred:
                    continue

                local_idx = rpos - ref_start
                if local_idx < 0 or local_idx >= len(ref_seq):
                    continue

                direction = is_cpg(ref_seq, local_idx)
                if direction == 0:
                    continue

                res = informative_call_and_output_strand(
                    qseq[qpos], direction, strand_class
                )
                if res is None:
                    continue

                out_strand, call = res
                writer.writerow(
                    [
                        read.query_name,
                        chrom,
                        rpos + 1,   # convert to 1-based reference coordinate
                        out_strand,
                        call,
                        q,
                        read_start,
                        read_end,
                    ]
                )

    bam.close()
    fasta.close()


# ---------------------------------------------------------------------------
# Sorting and duplicate collapse
# ---------------------------------------------------------------------------

def run_external_sort(input_tsv: str, sorted_tsv: str, tempdir: str) -> None:
    """
    Sort the stage-1 TSV using the system `sort` command.

    Sort order:
      1. fragment_name
      2. chrom
      3. ref_pos (numeric ascending)
      4. strand
      5. base_qual (numeric descending)

    The header row is removed before sorting.

    Notes
    -----
    External sorting is used to keep memory usage manageable on large inputs.
    """
    cmd = (
        f"tail -n +2 {shell_quote(input_tsv)} | "
        f"LC_ALL=C sort "
        f"-T {shell_quote(tempdir)} "
        f"-S 2G "
        f"-t $'\\t' "
        f"-k1,1 -k2,2 -k3,3n -k4,4 -k6,6nr "
        f"> {shell_quote(sorted_tsv)}"
    )
    subprocess.run(["bash", "-lc", cmd], check=True)


def choose_collapsed(
    rows: List[Tuple[str, str, str, str, str, str, str, str]],
    overlap_policy: str,
) -> Optional[Tuple[str, str, str, str, str, str, str, str]]:
    """
    Resolve multiple observations from the same fragment at the same CpG site.

    Parameters
    ----------
    rows : list of tuple
        Rows with identical (fragment_name, chrom, ref_pos, strand), already
        sorted by decreasing base quality.
    overlap_policy : {"agree", "best"}
        Conflict-resolution policy:
        - "agree": retain only if all calls agree
        - "best": retain the highest-quality observation even if calls conflict

    Returns
    -------
    Optional[tuple]
        The retained row, or None if the site should be discarded.
    """
    if len(rows) == 1:
        return rows[0]

    calls = {r[4] for r in rows}
    if len(calls) == 1:
        return rows[0]

    if overlap_policy == "best":
        return rows[0]

    return None


def stage2_stream_collapse(
    sorted_tsv: str,
    output_tsv: str,
    overlap_policy: str,
) -> None:
    """
    Stream through the sorted TSV and collapse within-fragment duplicate sites.

    A fragment is defined here by the sorted key (fragment_name, chrom).
    Each unique fragment is assigned a sequential integer `fragment_id`.

    For each fragment:
    - rows are grouped by (fragment_name, chrom, ref_pos, strand)
    - duplicate observations are collapsed via `choose_collapsed()`
    - fragment span is reported as:
          max(read_end) - min(read_start) + 1

    Final output columns
    --------------------
    fragment_id, chrom, ref_pos, strand, call, fragment_span
    """
    with open(sorted_tsv) as inp, open(output_tsv, "w", newline="") as out:
        writer = csv.writer(out, delimiter="\t", lineterminator="\n")
        writer.writerow(
            ["fragment_id", "chrom", "ref_pos", "strand", "call", "fragment_span"]
        )

        current_key = None
        current_rows: List[Tuple[str, str, str, str, str, str, str, str]] = []

        current_fragment_name = None
        current_fragment_id = 0
        current_fragment_chrom = None
        current_fragment_min_start = None
        current_fragment_max_end = None
        pending_rows = []

        def flush_site():
            """
            Collapse the currently buffered site and store it for fragment-level
            output if retained.
            """
            nonlocal current_rows, pending_rows
            if not current_rows:
                return
            chosen = choose_collapsed(current_rows, overlap_policy)
            if chosen is not None:
                pending_rows.append(chosen)
            current_rows = []

        def flush_fragment():
            """
            Write all retained collapsed rows for the current fragment.
            """
            nonlocal pending_rows
            if not pending_rows:
                return

            span = current_fragment_max_end - current_fragment_min_start + 1
            for row in pending_rows:
                writer.writerow(
                    [
                        current_fragment_id,
                        row[1],   # chrom
                        row[2],   # ref_pos
                        row[3],   # strand
                        row[4],   # call
                        span,
                    ]
                )
            pending_rows = []

        for line in inp:
            line = line.rstrip("\n")
            if not line:
                continue

            row = tuple(line.split("\t"))
            frag_name, chrom, ref_pos, strand, call, base_qual, read_start, read_end = row

            fragment_key = (frag_name, chrom)
            site_key = (frag_name, chrom, ref_pos, strand)

            if current_fragment_name is None:
                # Initialize the first fragment/site block.
                current_fragment_id = 1
                current_fragment_name = frag_name
                current_fragment_chrom = chrom
                current_fragment_min_start = int(read_start)
                current_fragment_max_end = int(read_end)
                current_key = site_key
                current_rows = [row]
                continue

            if fragment_key != (current_fragment_name, current_fragment_chrom):
                # Finished one fragment: finalize the last site, then write fragment.
                flush_site()
                flush_fragment()

                # Start a new fragment block.
                current_fragment_id += 1
                current_fragment_name = frag_name
                current_fragment_chrom = chrom
                current_fragment_min_start = int(read_start)
                current_fragment_max_end = int(read_end)
                current_key = site_key
                current_rows = [row]
                continue

            # Same fragment: update span bounds.
            rs = int(read_start)
            re = int(read_end)
            if rs < current_fragment_min_start:
                current_fragment_min_start = rs
            if re > current_fragment_max_end:
                current_fragment_max_end = re

            # If the CpG site changes, collapse the previous site group.
            if site_key != current_key:
                flush_site()
                current_key = site_key

            current_rows.append(row)

        # Flush the final buffered site and fragment.
        flush_site()
        flush_fragment()


# ---------------------------------------------------------------------------
# Command-line interface
# ---------------------------------------------------------------------------

def main() -> int:
    """
    Parse arguments and run the two-stage extraction/collapse workflow.

    Returns
    -------
    int
        Process exit code:
            0 on success
            1 on error
    """
    parser = argparse.ArgumentParser(
        description="Extract per-read per-CpG methylation calls with external sorting."
    )
    parser.add_argument("reference_fasta", help="Reference FASTA with .fai index")
    parser.add_argument("input_bam", help="Input BAM or - for stdin")
    parser.add_argument("output_tsv", help="Output TSV")
    parser.add_argument(
        "-p", "--min-phred", type=int, default=5, help="Minimum base quality"
    )
    parser.add_argument(
        "-q", "--min-mapq", type=int, default=10, help="Minimum mapping quality"
    )
    parser.add_argument(
        "--require-flags",
        type=lambda x: int(x, 0),
        default=0,
        help="Require these BAM flags (decimal or 0xHEX)",
    )
    parser.add_argument(
        "--ignore-flags",
        type=lambda x: int(x, 0),
        default=0x904,
        help="Ignore reads with any of these BAM flags set (default: 0x904)",
    )
    parser.add_argument(
        "--overlap-policy",
        choices=["agree", "best"],
        default="agree",
        help="How to resolve same-fragment same-CpG conflicts",
    )
    parser.add_argument(
        "--tempdir",
        default=None,
        help="Temporary directory for stage1 and sort files",
    )
    args = parser.parse_args()

    tempdir = args.tempdir or tempfile.gettempdir()
    os.makedirs(tempdir, exist_ok=True)

    stage1_tsv = None
    sorted_tsv = None

    try:
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".stage1.tsv", delete=False, dir=tempdir
        ) as tf1:
            stage1_tsv = tf1.name

        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".sorted.tsv", delete=False, dir=tempdir
        ) as tf2:
            sorted_tsv = tf2.name

        stage1_extract(
            bam_path=args.input_bam,
            fasta_path=args.reference_fasta,
            tmp_tsv=stage1_tsv,
            min_phred=args.min_phred,
            min_mapq=args.min_mapq,
            require_flags=args.require_flags,
            ignore_flags=args.ignore_flags,
        )

        run_external_sort(stage1_tsv, sorted_tsv, tempdir)

        stage2_stream_collapse(
            sorted_tsv=sorted_tsv,
            output_tsv=args.output_tsv,
            overlap_policy=args.overlap_policy,
        )

    except Exception as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 1

    finally:
        # Best-effort cleanup of temporary files.
        for path in (stage1_tsv, sorted_tsv):
            if path:
                try:
                    os.unlink(path)
                except OSError:
                    pass

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
