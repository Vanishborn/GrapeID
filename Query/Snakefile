# Turns query sequence reads into consensus sequences

import os

configfile: "config.yml"

SAMPLES = config["samples"]
REF = config["reference"]
REF_PREFIX = "Genome/Index/ref"
QUERY_DIR = config["query_dir"]
SAMPLE_NAMES = list(SAMPLES.keys())

# Helper function to resolve FASTQ paths
def fq_path(sample):
	fq = SAMPLES[sample]
	if isinstance(fq, list):
		return [os.path.join(QUERY_DIR, f) for f in fq]
	else:
		return os.path.join(QUERY_DIR, fq)

rule all:
	input:
		expand("Results/Consensus/{sample}.fa", sample=SAMPLE_NAMES)

rule hisat2_build:
	input: REF
	output:
		expand("Genome/Index/ref.{i}.ht2", i=range(1, 9))
	threads: config["threads"]["hisat2_build"]
	shell:
		"hisat2-build -p {threads} {input} Genome/Index/ref"

rule align_all:
	input:
		expand("Results/Align/{sample}.sam", sample=SAMPLE_NAMES)

rule hisat2_align:
	input:
		index = expand("Genome/Index/ref.{i}.ht2", i=range(1, 9))
	output:
		temp("Results/Align/{sample}.sam")
	wildcard_constraints:
		sample="|".join(SAMPLE_NAMES)
	threads: config["threads"]["hisat2_align"]
	run:
		fq = fq_path(wildcards.sample)
		if isinstance(fq, list):
			shell("hisat2 -p {threads} -x Genome/Index/ref -1 {fq[0]} -2 {fq[1]} -S {output}")
		else:
			shell("hisat2 -p {threads} -x Genome/Index/ref -U {fq} -S {output}")

rule samtools_view:
	input: "Results/Align/{sample}.sam"
	output: temp("Results/BAM/{sample}.bam")
	wildcard_constraints:
		sample="|".join(SAMPLE_NAMES)
	threads: config["threads"]["samtools_view"]
	shell:
		"samtools view -@ {threads} -b -o {output} {input}"

rule samtools_fixmate:
	input: "Results/BAM/{sample}.bam"
	output: temp("Results/BAM/{sample}.fixmate.bam")
	wildcard_constraints:
		sample="|".join(SAMPLE_NAMES)
	threads: config["threads"]["samtools_fixmate"]
	shell:
		"samtools fixmate -@ {threads} -m {input} {output}"

rule samtools_sort:
	input: "Results/BAM/{sample}.fixmate.bam"
	output: temp("Results/BAM/{sample}.sorted.bam")
	wildcard_constraints:
		sample="|".join(SAMPLE_NAMES)
	threads: config["threads"]["samtools_sort"]
	shell:
		"samtools sort -@ {threads} {input} -o {output}"

rule samtools_markdup:
	input: "Results/BAM/{sample}.sorted.bam"
	output: temp("Results/BAM/{sample}.markdup.bam")
	wildcard_constraints:
		sample="|".join(SAMPLE_NAMES)
	threads: config["threads"]["samtools_markdup"]
	shell:
		"samtools markdup -@ {threads} {input} {output}"

rule samtools_consensus:
	input: "Results/BAM/{sample}.markdup.bam"
	output: "Results/Consensus/{sample}.fa"
	wildcard_constraints:
		sample="|".join(SAMPLE_NAMES)
	threads: config["threads"]["samtools_consensus"]
	shell:
		"samtools consensus -@ {threads} -m simple --show-del yes --show-ins no --ambig -d 0 -l 80 -f fasta {input} -o {output}"
