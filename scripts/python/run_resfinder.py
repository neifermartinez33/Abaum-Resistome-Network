import subprocess
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor

GENOMES_DIR = Path("/home/abaum-gcloud/genomes")
ANN_DIR     = Path("/home/abaum-gcloud/annotation/resfinder")
RES_DB      = Path("/home/abaum-gcloud/resfinder_db")
POINT_DB    = Path("/home/abaum-gcloud/pointfinder_db")
ANN_DIR.mkdir(parents=True, exist_ok=True)

all_fna = []
for source in ["ncbi", "patric", "pubmlst"]:
    all_fna.extend(list((GENOMES_DIR / source).glob("*.fna")))

print(f"Total genomes: {len(all_fna)}", flush=True)

def annotate(fna):
    name    = fna.stem
    out_dir = ANN_DIR / name
    if out_dir.exists() and any(out_dir.iterdir()):
        return "skipped"
    out_dir.mkdir(parents=True, exist_ok=True)
    cmd = ["python3", "-m", "resfinder",
           "-ifa", str(fna),
           "-o", str(out_dir),
           "-s", "Acinetobacter baumannii",
           "-l", "0.6", "-t", "0.8",
           "--acquired",
           "-db_res", str(RES_DB),
           "-db_point", str(POINT_DB)]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        return f"failed:{name}"
    return "ok"

failed = []
done = 0

with ProcessPoolExecutor(max_workers=32) as executor:
    for result in executor.map(annotate, all_fna):
        done += 1
        if "failed" in result:
            failed.append(result)
        if done % 500 == 0:
            print(f"Progress: {done}/{len(all_fna)} done, {len(failed)} failed", flush=True)

print(f"Done. Total: {len(all_fna)}, Failed: {len(failed)}", flush=True)
