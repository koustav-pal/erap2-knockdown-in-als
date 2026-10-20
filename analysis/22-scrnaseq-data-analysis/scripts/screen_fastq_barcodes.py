#!/usr/bin/env python3
"""Annotate an nf-core/scrnaseq CSV with exact barcode-match percentages.

Python >=3.8; standard library only. Local plain/gzip FASTQs and barcode lists.
Default: screen the first 100,000 records of EVERY distinct FASTQ, matching
listed barcodes at position 1. --max-reads 0 scans every record. This is prefix
sampling, not a random sample; use a full scan for definitive percentages.

Examples:
  python3 screen_fastq_barcodes.py -i samples.csv -b whitelist.txt.gz -o screened.csv
  python3 screen_fastq_barcodes.py -i samples.csv -b whitelist.txt -o all.csv --max-reads 0
  python3 screen_fastq_barcodes.py -i samples.csv -b list.tsv -o scan.csv --anywhere --orientation both
  python3 screen_fastq_barcodes.py -i samples.csv -b list.csv -o scan.csv --barcode-column 2 --barcode-skip-rows 1
  python3 screen_fastq_barcodes.py -i samples.csv -b list.txt -o scan.csv --positions 1,39,77

Positions are 1-based starts in the ORIGINAL read. Reverse matching compares
the reverse complement of each candidate window against the list. 'both'
reports forward, reverse, and their union, counting each read once per metric.
Mixed barcode lengths are supported. Matching is exact, case-insensitive;
N/IUPAC symbols in reads do not match A/C/G/T barcode windows. No mismatch
correction, indels, trimming, chemistry inference, read swapping or cell calling.

Whitelist: one barcode per line, or select a 1-based CSV/TSV/whitespace column.
Use --barcode-skip-rows for headers. Blank lines and # comments are ignored.
Only A/C/G/T sequences are accepted (no automatic suffix removal). FASTA is
not accepted. Gzip is detected by magic bytes, not by filename.

For Parse or other split-barcode assays, supply an individual component list
and the correct component positions, or --anywhere for exploratory screening.
A hit means ANY listed contiguous component occurs at ANY tested position;
it does NOT validate a complete combinatorial cell barcode or its linkers.
Run separately for distinct round-specific lists. Short barcodes tested anywhere
can match by chance. Use known positions for interpretable assay-specific QC.

CSV rows and original columns are preserved. Additional columns are prefixed
barcode_screen_ and fastq_1_/fastq_2_. Percentages are 0-100, not fractions,
and use all screened reads as denominator, including short/ambiguous reads.
'eligible_reads' counts reads long enough for at least one requested window.
An empty or malformed FASTQ yields status=error and blank counts/percentages;
other files continue. Exit codes: 0 success, 1 input/config error, 2 file errors
with an annotated report saved. Missing fastq_2 yields status=missing, not 0%.
Only conventional four-line FASTQ records are supported; headers, separator,
sequence/quality lengths, and quality characters are checked in scanned records.

Relative FASTQ paths resolve against the samplesheet directory unless --base-dir
is provided. URLs and globs are rejected. --workers uses threads sharing one
whitelist in memory; speedup depends on storage/CPU. Repeated paths are scanned
once. A multi-million-entry whitelist can use several hundred MB of RAM.
Output must be a NEW file and cannot overwrite the input or previous report.
"""
import argparse
import csv
import gzip
import os
import re
import sys
import tempfile
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

COMP = str.maketrans('ACGTN', 'TGCAN')
METRICS = ('screened_reads', 'eligible_reads', 'matched_reads', 'matched_pct',
           'forward_matched_reads', 'forward_matched_pct',
           'reverse_matched_reads', 'reverse_matched_pct',
           'min_read_length', 'max_read_length', 'status', 'error')


def open_text(path):
    with open(path, 'rb') as handle:
        compressed = handle.read(2) == b'\x1f\x8b'
    return (gzip.open if compressed else open)(path, 'rt', encoding='utf-8-sig', newline='')


def load_barcodes(args):
    grouped = {}
    with open_text(args.barcodes) as handle:
        for line_no, raw in enumerate(handle, 1):
            if line_no <= args.barcode_skip_rows or not raw.strip() or raw.lstrip().startswith('#'):
                continue
            if '\t' in raw:
                fields = next(csv.reader([raw], delimiter='\t'))
            elif ',' in raw:
                fields = next(csv.reader([raw]))
            else:
                fields = raw.split()
            try:
                barcode = fields[args.barcode_column - 1].strip().upper()
            except IndexError:
                raise ValueError('Barcode line {} lacks column {}'.format(line_no, args.barcode_column))
            if not re.fullmatch('[ACGT]+', barcode):
                raise ValueError('Invalid barcode at line {}: {!r}; use a sequence column and skip headers'.format(line_no, barcode))
            grouped.setdefault(len(barcode), set()).add(barcode)
    if not grouped:
        raise ValueError('Barcode list is empty')
    return grouped


def fastq_sequences(handle):
    record = 0
    while True:
        header = handle.readline()
        if not header:
            return
        record += 1
        seq, plus, qual = handle.readline(), handle.readline(), handle.readline()
        if not seq or not plus or not qual:
            raise ValueError('Truncated FASTQ record {}'.format(record))
        seq, qual = seq.rstrip('\r\n'), qual.rstrip('\r\n')
        if not header.startswith('@') or not plus.startswith('+'):
            raise ValueError('Invalid FASTQ header/separator at record {}'.format(record))
        if not seq or len(seq) != len(qual) or any(c.isspace() for c in seq):
            raise ValueError('Invalid sequence/quality lengths at record {}'.format(record))
        if any(ord(c) < 33 or ord(c) > 126 for c in qual):
            raise ValueError('Invalid quality characters at record {}'.format(record))
        yield seq.upper()


def match_read(seq, grouped, positions, anywhere, orientation):
    forward = reverse = eligible = False
    for length, whitelist in grouped.items():
        starts = range(len(seq) - length + 1) if anywhere else positions
        for start in starts:
            if start + length > len(seq):
                continue
            eligible = True
            window = seq[start:start + length]
            if orientation != 'reverse' and not forward:
                forward = window in whitelist
            if orientation != 'forward' and not reverse:
                reverse = window.translate(COMP)[::-1] in whitelist
            if ((orientation == 'forward' and forward) or
                (orientation == 'reverse' and reverse) or (forward and reverse)):
                return forward, reverse, eligible
    return forward, reverse, eligible


def screen(path, grouped, args):
    n = eligible = hits = forward = reverse = 0
    shortest = None
    longest = 0
    try:
        with open_text(path) as handle:
            for seq in fastq_sequences(handle):
                f, r, e = match_read(seq, grouped, args.positions, args.anywhere, args.orientation)
                n += 1
                eligible += e
                hits += f or r
                forward += f
                reverse += r
                shortest = len(seq) if shortest is None else min(shortest, len(seq))
                longest = max(longest, len(seq))
                if args.max_reads and n >= args.max_reads:
                    break
        if not n:
            raise ValueError('Empty FASTQ')
        pct = lambda k: '{:.6f}'.format(100 * k / n)
        return dict(zip(METRICS, (n, eligible, hits, pct(hits),
                    forward if args.orientation != 'reverse' else '',
                    pct(forward) if args.orientation != 'reverse' else '',
                    reverse if args.orientation != 'forward' else '',
                    pct(reverse) if args.orientation != 'forward' else '',
                    shortest, longest, 'ok', '')))
    except (OSError, EOFError, UnicodeError, ValueError) as exc:
        result = dict.fromkeys(METRICS, '')
        result.update(status='error', error=str(exc))
        return result


def resolve_fastq(value, base):
    if not value or not value.strip():
        return None
    value = os.path.expandvars(os.path.expanduser(value.strip()))
    if '://' in value or any(c in value for c in '*?[]'):
        raise ValueError('Only explicit local FASTQ paths are supported: {}'.format(value))
    path = Path(value)
    return str((path if path.is_absolute() else base / path).resolve())


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-i', '--input', required=True, type=Path)
    parser.add_argument('-b', '--barcodes', required=True, type=Path)
    parser.add_argument('-o', '--output', required=True, type=Path)
    parser.add_argument('--max-reads', type=int, default=100000, help='Records per FASTQ; 0=all (default: 100000)')
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--positions', default='1', help='Comma-separated 1-based barcode starts (default: 1)')
    group.add_argument('--anywhere', action='store_true', help='Test every possible start; beware chance matches')
    parser.add_argument('--orientation', choices=['forward', 'reverse', 'both'], default='forward')
    parser.add_argument('--barcode-column', type=int, default=1)
    parser.add_argument('--barcode-skip-rows', type=int, default=0)
    parser.add_argument('--base-dir', type=Path)
    parser.add_argument('--workers', type=int, default=1)
    args = parser.parse_args()
    try:
        args.positions = sorted(set(int(x) - 1 for x in args.positions.split(',')))
        if min(args.positions) < 0 or args.max_reads < 0 or args.workers < 1 or args.barcode_column < 1 or args.barcode_skip_rows < 0:
            raise ValueError('Invalid negative/zero option value')
        if args.output.exists() or args.output.is_symlink():
            raise ValueError('Output already exists; choose a new output filename')
        with args.input.open(newline='', encoding='utf-8-sig') as handle:
            reader = csv.DictReader(handle)
            fields = reader.fieldnames
            if not fields or len(fields) != len(set(fields)):
                raise ValueError('CSV requires unique column names')
            if not {'sample', 'fastq_1', 'fastq_2'}.issubset(fields):
                raise ValueError('CSV must contain sample, fastq_1, fastq_2')
            rows = list(reader)
        if not rows or any(None in row or any(v is None for v in row.values()) for row in rows):
            raise ValueError('Empty samplesheet or rows with incorrect column counts')
        metadata = {
            'barcode_screen_list': str(args.barcodes.resolve()),
            'barcode_screen_positions': 'anywhere' if args.anywhere else ','.join(str(x+1) for x in args.positions),
            'barcode_screen_orientation': args.orientation,
            'barcode_screen_max_reads': args.max_reads,
            'barcode_screen_matching': 'exact',
        }
        extra = list(metadata) + ['{}_{}'.format(m, k) for m in ('fastq_1', 'fastq_2') for k in METRICS]
        if set(extra) & set(fields):
            raise ValueError('Input already contains screening columns; use the original samplesheet')
        base = (args.base_dir or args.input.resolve().parent).resolve()
        pairs = [(resolve_fastq(r['fastq_1'], base), resolve_fastq(r['fastq_2'], base)) for r in rows]
        paths = sorted({p for pair in pairs for p in pair if p})
        grouped = load_barcodes(args)
        print('Loaded {} unique barcodes; lengths {}; {} distinct FASTQs'.format(
            sum(map(len, grouped.values())), sorted(grouped), len(paths)), file=sys.stderr)
        results = {}
        with ThreadPoolExecutor(max_workers=args.workers) as pool:
            futures = {pool.submit(screen, p, grouped, args): p for p in paths}
            for future in as_completed(futures):
                path = futures[future]
                result = future.result()
                results[path] = result
                print('[{}/{}] {}: {} {}'.format(len(results), len(paths), path,
                      result['status'], result['matched_pct'] if result['status'] == 'ok' else result['error']), file=sys.stderr)
        errors = sum(r['status'] == 'error' for r in results.values())
        for row, pair in zip(rows, pairs):
            row.update(metadata)
            for mate, path in zip(('fastq_1', 'fastq_2'), pair):
                result = results[path] if path else dict(dict.fromkeys(METRICS, ''), status='missing')
                if path is None and mate == 'fastq_1':
                    errors += 1
                row.update({'{}_{}'.format(mate, k): v for k, v in result.items()})
        # Write completely, then publish without overwriting an existing destination.
        temp_name = None
        try:
            with tempfile.NamedTemporaryFile(mode='w', newline='', encoding='utf-8',
                    dir=args.output.resolve().parent, delete=False) as handle:
                temp_name = handle.name
                writer = csv.DictWriter(handle, fieldnames=fields + extra)
                writer.writeheader()
                writer.writerows(rows)
            os.link(temp_name, args.output)
        finally:
            if temp_name:
                os.unlink(temp_name)
        print('Saved {}; {} file errors/missing required mates'.format(args.output, errors), file=sys.stderr)
        return 2 if errors else 0
    except (OSError, ValueError, csv.Error) as exc:
        print('ERROR: {}'.format(exc), file=sys.stderr)
        return 1


if __name__ == '__main__':
    sys.exit(main())
