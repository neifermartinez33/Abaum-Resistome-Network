#!/bin/bash
conda run -n abaum_qc fastANI \
    --ql /media/miguel/MEMORY_333/Abaum_Resistome_Network/data/processed/ani/query_list.txt \
    --ref /media/miguel/MEMORY_333/Abaum_Resistome_Network/data/processed/ncbi_fasta/GCF_000018445.1.fna \
    --output /media/miguel/MEMORY_333/Abaum_Resistome_Network/data/processed/ani/ani_results.txt \
    --threads 8 \
    --minFraction 0.5 \
    > /home/miguel/Abaum_Resistome_Network/fastani.log 2>&1
echo "FastANI finished"
