for file in ./01_raw/*.fq.gz; do
    base=$(basename "$file" .fq.gz)
    gunzip -k "$file"                # -k keeps the original .gz just in case
    mv "./01_raw/$base.fq" "./01_raw/$base.fastq"
done