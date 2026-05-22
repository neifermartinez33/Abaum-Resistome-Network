import subprocess
import pandas as pd
from pathlib import Path

GENOMES_DIR = Path("/home/abaum-gcloud/genomes")
ANN_DIR     = Path("/home/abaum-gcloud/annotation/amrfinderplus")
ANN_DIR.mkdir(parents=True, exist_ok=True)

sources = {"ncbi": GENOMES_DIR/"ncbi", "patric": GENOMES_DIR/"patric", "pubmlst": GENOMES_DIR/"pubmlst"}

failed = []
total = 0
done = 0

for source, fna_dir in sources.items():
    fna_files = list(fna_dir.glob("*.fna"))
    total += len(fna_files)
    for i, fna in enumerate(fna_files):
        name    = fna.stem
        out_tsv = ANN_DIR / f"{name}.tsv"
        if out_tsv.exists() and out_tsv.stat().st_size > 0:
            done += 1
            continue
        cmd = ["amrfinder", "--nucleotide", str(fna),
               "--organism", "Acinetobacter_baumannii",
               "--output", str(out_tsv), "--threads", "8", "--plus"]
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            failed.append(str(fna))
        done += 1
        if done % 500 == 0:
            print(f"Progress: {done}/{total} done, {len(failed)} failed", flush=True)

print(f"Done. Total: {total}, Failed: {len(failed)}", flush=True)
