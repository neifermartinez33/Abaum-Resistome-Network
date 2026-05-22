import subprocess
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor

ANN_DIR = Path("/home/abaum-gcloud/annotation/card_rgi")
ANN_DIR.mkdir(parents=True, exist_ok=True)

all_fna = []
for folder in ["ncbi", "patric", "pubmlst"]:
    all_fna.extend(list(Path(f"/home/abaum-gcloud/genomes/{folder}").glob("*.fna")))

print(f"Total genomes to annotate: {len(all_fna)}", flush=True)

def annotate(fna):
    name = fna.stem
    out_json = ANN_DIR / f"{name}.json"
    if out_json.exists() and out_json.stat().st_size > 0:
        return "skipped"
    cmd = [
        "rgi", "main",
        "--input_sequence", str(fna),
        "--output_file", str(ANN_DIR / name),
        "--input_type", "contig",
        "--alignment_tool", "BLAST",
        "--num_threads", "1",
        "--include_loose",
        "--clean"
    ]
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
