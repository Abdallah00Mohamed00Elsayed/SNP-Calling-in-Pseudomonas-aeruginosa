# SNP-Calling-in-Pseudomonas-aeruginosa


This pipeline identifies single nucleotide polymorphisms (SNPs) in paired-end sequencing data using standard bioinformatics tools.

---

## 1. Reference Genome Download

Download the reference genome (NC_002516.2) in FASTA format from NCBI.

---

## 2. Sample Data Download

Use Galaxy (tool: Faster Download) to retrieve paired-end reads for:

- SRR8737536 (forward and reverse)
- SRR8737537 (forward and reverse)

---

## 3. Create and Activate Conda Environment

Activate the environment that includes all required tools (e.g., `bwa`, `samtools`, `bcftools`, `fastqc`, `trimmomatic`, `snpEff`).

```bash
conda activate SNP_calling
```

---

## 4. Quality Control (FastQC)

Run FastQC to assess raw read quality:

```bash
fastqc SRR*.fastqsanger.gz
```

Check metrics like:

- Per base sequence quality
- GC content
- Sequence length distribution

SRR8737537 reverse reads often require trimming due to low quality.

---

## 5. Trimming (Trimmomatic)

Apply trimming only on SRR8737537 to remove low-quality bases and adapter sequences:

```bash
java -jar $(find ~ -name "trimmomatic*.jar" | head -1) PE \
  -phred33 \
  SRR37_f.fastqsanger.gz SRR37_r.fastqsanger.gz \
  SRR37_trimmed_f.fastq.gz /dev/null \
  SRR37_trimmed_r.fastq.gz /dev/null \
  ILLUMINACLIP:$(find ~ -name "TruSeq3-PE.fa" | head -1):2:30:10 \
  LEADING:30 TRAILING:30 SLIDINGWINDOW:4:30 MINLEN:36
```

---

## 6. Quality Control After Trimming

Re-run FastQC on the trimmed reads to confirm quality improvement.

---

## 7. Index the Reference Genome (BWA)

```bash
bwa index Reference.fasta
```

---

## 8. Align Reads to Reference (BWA-MEM)

```bash
# SRR36 (raw)
bwa mem Reference.fasta SRR36_f.fastqsanger.gz SRR36_r.fastqsanger.gz > SRR36_aligned.sam

# SRR37 (trimmed)
bwa mem Reference.fasta SRR37_trimmed_f.fastq.gz SRR37_trimmed_r.fastq.gz > SRR37_aligned.sam
```

---

## 9. Convert SAM to BAM (Samtools)

```bash
samtools view -S -b -q 20 SRR36_aligned.sam > SRR36_aligned.bam
samtools view -S -b -q 20 SRR37_aligned.sam > SRR37_aligned.bam
```

---

## 10. Sort, Index and Generate Stats

```bash
parallel "samtools sort {} -o {.}.sorted.bam && samtools index {.}.sorted.bam && samtools stats {.}.sorted.bam > {.}.stats.txt" ::: SRR36_aligned.bam SRR37_aligned.bam
```

---

## 11. Variant Calling (BCFtools)

```bash
bcftools mpileup -Ou -f Reference.fasta SRR36_aligned.sorted.bam | bcftools call -mv -Ob -o variants_36.bcf
bcftools mpileup -Ou -f Reference.fasta SRR37_aligned.sorted.bam | bcftools call -mv -Ob -o variants_37.bcf
```

---

## 12. Convert BCF to VCF

```bash
bcftools view variants_36.bcf > variants_36.vcf
bcftools view variants_37.bcf > variants_37.vcf
```

---

## 13. Filter Variants

Keep only high-quality variants:

```bash
bcftools filter -i 'QUAL>20 && DP>10' variants_36.vcf > filtered_variants_36.vcf
bcftools filter -i 'QUAL>20 && DP>10' variants_37.vcf > filtered_variants_37.vcf
```

---

## 14. Rename Chromosome Field (for SnpEff compatibility)

```bash
sed 's/^NC_002516.2/Chromosome/' filtered_variants_36.vcf > output_36.vcf
sed 's/^NC_002516.2/Chromosome/' filtered_variants_37.vcf > output_37.vcf
```

---

## 15. Variant Annotation (SnpEff)

```bash
java -Xmx8g -jar snpEff.jar -v Pseudomonas_aeruginosa_pao1_gca_000006765 output_36.vcf > annotated_variants_36.vcf
java -Xmx8g -jar snpEff.jar -v Pseudomonas_aeruginosa_pao1_gca_000006765 output_37.vcf > annotated_variants_37.vcf
```

---

## Output

Final output includes:

- Annotated VCFs
- BAM files
- Summary statistics
- Variant impact classification (MODIFIER, LOW, MODERATE, HIGH)

This pipeline enables SNP discovery and functional analysis in bacterial clinical isolates.
