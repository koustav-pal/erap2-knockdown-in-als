"""Stage manifest inputs from a source tree or SSH host; never overwrite source data."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import shlex
import shutil
import subprocess

REPO = Path(__file__).resolve().parents[2]


def copy_stream(source, target, gene_rows=False):
    if gene_rows:
        for line in source:
            fields = line.split(b"\t", 3)
            if len(fields) >= 3 and fields[2] == b"gene":
                target.write(line)
    else:
        shutil.copyfileobj(source, target)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--source-root", required=True, help="Existing project tree")
    p.add_argument("--host", help="Optional user@SSH-host")
    p.add_argument("--data-root", default="workflow/data")
    p.add_argument("--files", nargs="+", help="Only these manifest filenames")
    args = p.parse_args()
    data = Path(args.data_root).resolve()
    data.mkdir(parents=True, exist_ok=True)
    records = []
    with (REPO / "workflow/config/inputs.tsv").open() as stream:
        rows = list(csv.DictReader(stream, delimiter="\t"))
    if args.files and set(args.files) - {row["file"] for row in rows}:
        p.error("Unknown filename in --files")
    for row in rows:
        if args.files and row["file"] not in args.files:
            continue
        source = row["source"]
        if not source.startswith("/"):
            source = str(Path(args.source_root) / source)
        target = data / row["file"]
        temporary = target.with_suffix(target.suffix + ".partial")
        gene_rows = row["transform"] == "gtf_gene_rows"
        print(f"Stage {row['file']} from {args.host or 'local'}:{source}", flush=True)
        try:
            if args.host:
                # Filter annotation on server to avoid transferring a 1.4 GB GTF.
                remote = ("awk -F '\t' '$3 == \"gene\"' < " if gene_rows else "cat -- ") + shlex.quote(source)
                with temporary.open("wb") as out:
                    subprocess.run(["ssh", "-o", "BatchMode=yes", "-o", "ConnectTimeout=10",
                                    args.host, remote], stdout=out, check=True)
            else:
                with open(source, "rb") as src, temporary.open("wb") as out:
                    copy_stream(src, out, gene_rows)
            if not temporary.stat().st_size:
                raise ValueError(f"Source produced empty file: {source}")
            temporary.replace(target)
        finally:
            temporary.unlink(missing_ok=True)
        h = hashlib.sha256()
        with target.open("rb") as stream:
            for chunk in iter(lambda: stream.read(8 * 1024 * 1024), b""):
                h.update(chunk)
        records.append(dict(row, host=args.host, bytes=target.stat().st_size, sha256=h.hexdigest()))
    # Merge subset receipts instead of losing provenance from earlier transfers.
    receipt = data / "materialization.json"
    previous = json.loads(receipt.read_text()) if receipt.exists() else []
    current = {row["file"]: row for row in previous}
    current.update({row["file"]: row for row in records})
    receipt.write_text(json.dumps(list(current.values()), indent=2) + "\n")


if __name__ == "__main__":
    main()
