import subprocess
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor

GENOMES_DIR = Path("/home/abaum-gcloud/genomes/pubmlst_missing")
ANN_DIR     = Path("/home/abaum-gcloud/annotation/amrfinderplus")
ANN_DIR.mkdir(parents=True, exist_ok=True)

all_fna = list(GENOMES_DIR.glob("*.fna"))
print(f"Total genomes to annotate: {len(all_fna)}", flush=True)

def annotate(fna):
    name    = fna.stem
    out_tsv = ANN_DIR / f"{name}.tsv"
    if out_tsv.exists() and out_tsv.stat().st_size > 0:
        return "skipped"
    cmd = ["amrfinder", "--nucleotide", str(fna),
           "--organism", "Acinetobacter_baumannii",
           "--output", str(out_tsv), "--threads", "1", "--plus"]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        return f"failed:{fna}"
    return "ok"

failed = []
done = 0

with ProcessPoolExecutor(max_workers=32) as executor:
    for result in executor.map(annotate, all_fna):
        done += 1
        if "failed" in result:
            failed.append(result)
        if done % 200 == 0:
            print(f"Progress: {done}/{len(all_fna)} done, {len(failed)} failed", flush=True)

print(f"Done. Total: {len(all_fna)}, Failed: {len(failed)}", flush=True)
