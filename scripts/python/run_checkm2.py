import os
import shutil
import subprocess
from pathlib import Path

NCBI_FASTA_DIR = Path("/media/miguel/MEMORY_333/Abaum_Resistome_Network/data/processed/ncbi_fasta")
BATCH_OUTPUT_DIR = Path("/media/miguel/MEMORY_333/Abaum_Resistome_Network/data/processed/checkm2_batch_output")
BATCH_INPUT_DIR = Path("/tmp/checkm2_batch_input")
BATCH_OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
BATCH_INPUT_DIR.mkdir(exist_ok=True)

BATCH_SIZE = 500
THREADS = 6

ncbi_genomes = sorted(NCBI_FASTA_DIR.glob("*.fna"))
batches = [ncbi_genomes[i:i+BATCH_SIZE] for i in range(0, len(ncbi_genomes), BATCH_SIZE)]

print(f"Total genomes : {len(ncbi_genomes)}", flush=True)
print(f"Total batches : {len(batches)}", flush=True)

failed = []

for i, batch in enumerate(batches):
    batch_out = BATCH_OUTPUT_DIR / f"batch_{i:03d}"
    report = batch_out / "quality_report.tsv"
    if report.exists():
        print(f"Batch {i:03d} already done, skipping", flush=True)
        continue
    print(f"Batch {i:03d} processing {len(batch)} genomes...", flush=True)
    batch_in = BATCH_INPUT_DIR / f"batch_{i:03d}"
    if batch_in.exists():
        shutil.rmtree(batch_in)
    batch_in.mkdir()
    for fna in batch:
        os.symlink(fna, batch_in / fna.name)
    batch_out.mkdir(exist_ok=True)
    cmd = ["checkm2", "predict", "--input", str(batch_in),
           "--output-directory", str(batch_out),
           "--threads", str(THREADS), "--force", "-x", "fna"]
    result = subprocess.run(cmd)
    if result.returncode != 0:
        print(f"Batch {i:03d} FAILED", flush=True)
        failed.append(i)
    else:
        print(f"Batch {i:03d} OK", flush=True)
    shutil.rmtree(batch_in)

print(f"Done. Failed batches: {failed}", flush=True)
