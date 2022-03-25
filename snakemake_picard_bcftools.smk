configfile: "config.yaml"
import glob, os

files=glob.glob(config["Reads"])
print(files)

SAMPLENAME=[os.path.basename(x) for x in files]
DIR=[os.path.dirname(files[0])]

SAMPLES=[]
SAMPLENUM=[]
LANE=[]

for x in SAMPLENAME:
	if re.match('.*_?.*_L\d+.*',x) is not None:
		SAMPLES.append(re.match('(.*)_.*_L\d+.*',x).group(1))
		# SAMPLENUM.append(re.match('.*_(.*)_L\d+.*',x).group(1))
		LANE.append(re.match('.*_(.*_L\d+)_.*',x).group(1))

SAMPLES_UNIQ_L=[]
SAMPLES_UNIQ=set(SAMPLES)
for x in SAMPLES_UNIQ:
	SAMPLES_UNIQ_L.append(x)

DIR_UNIQ_L=[]
DIR_UNIQ=set(DIR)
for x in DIR_UNIQ:
	DIR_UNIQ_L.append(x)
print(SAMPLES_UNIQ_L, "SAMPLES_UNIQ_L", LANE, "LANE")

rule all:
	input:
		expand("{dirname}/bcftools_norm/{samplename}.la.md.bqsr.bcftools.normed.bcf", samplename=SAMPLES_UNIQ_L,dirname=DIR_UNIQ_L)

# TODO: fastqc report, summary stats

rule bwa_mem_map:
	input:
		 "{dirname}/{samplename}_{lane}_R1_001.fastq.gz"
	output:
		temp("{dirname}/sorted_bams/{samplename}_{lane}_sorted.bam")
	wildcard_constraints:
		lane=".*_L\d+"
	log:
		"{dirname}/logs/sorted_bams/{samplename}_{lane}_bwa_mem.log"
	run:
		read1=input[0]
		read2=input[0].replace('_R1_','_R2_')
		print(input[0], "input[0]")
		print(read1, read2, "read1 and 2")
		shell(
			"""
			header=$(zcat {input} | head -1)
			seqID=$(echo $header | cut -d':' -f 1 | cut -d' ' -f 1 | cut -c 2-)
			id=$(echo $header | cut -d':' -f 2-4,10 --output-delimiter='_')
			sm=$(echo {wildcards.samplename} | cut -d'_' -f 1)
			pu=$(echo $header | cut -d':' -f 3,4 --output-delimiter='_')
			rg=$(echo "@RG\tID:$seqID"_"$id\tSM:$sm\tLB:$sm\tPL:ILLUMINA")
			bwa mem -R $(echo "@RG\\tID:$sm"_"$seqID"_"$id\\tSM:$sm\\tLB:$sm\\tPL:ILLUMINA\\tPU:$pu") \
			-t {config[bwamemParams][bwaThreads]} \
			-M {config[refs][Reference]} \
			{read1} {read2} | \
			samtools view -bh - | \
			samtools sort -@ 2 -m 50G - > {output} && samtools index {output} 2> {log}
			"""
		)

def merge_inputs(wildcards):
	files = expand("{dirname}/sorted_bams/{samplename}_{lane}_sorted.bam", samplename=SAMPLES_UNIQ_L, lane=LANE,
			   dirname=DIR_UNIQ_L)
	print(files, "files")
	return files

rule merged_bams:
	input:
		merge_inputs
	output:
		"{dirname}/merged_bams/{samplename}.bam"
	run:
		inputs_all= " ".join(["I=" + repr(i) for i in input])
		shell("java -jar  /usr/local/src/picard.jar MergeSamFiles CREATE_INDEX=true O={output} {inputs_all}")

rule leftalign_bams:
	input:
		"{dirname}/merged_bams/{samplename}.bam"
	output:
		temp("{dirname}/leftalign_bams/{samplename}.la.bam")
	log:
		"{dirname}/logs/leftalign_bams/{samplename}.la.log"
	shell:
		"gatk LeftAlignIndels --verbosity ERROR --QUIET true -R {config[refs][Reference]} -I {input} -O {output} 2> {log}"


rule markdup_bams:
	input:
		"{dirname}/leftalign_bams/{samplename}.la.bam"
	output:
		bam="{dirname}/markdup_bams/{samplename}.la.md.bam",
		metrics="{dirname}/markdup_bams/{samplename}.la.md.metrics"
	shell:
		"java -Xmx225G -jar  /usr/local/src/picard.jar  MarkDuplicates COMPRESSION_LEVEL=9 CREATE_INDEX=true "
		 "VERBOSITY=ERROR QUIET=true I={input} O={output.bam} M={output.metrics} "
		"TMP_DIR=/mnt/scratch0/samtmp"

rule bqsr_baserecalib:
	input:
		bam="{dirname}/markdup_bams/{samplename}.la.md.bam",
		metrics="{dirname}/markdup_bams/{samplename}.la.md.metrics"
	output:
		"{dirname}/happy_bqsr_bams/{samplename}.recal.csv"
	shell:
		"gatk --java-options {config[bqsrParam][javaoptions]} BaseRecalibrator -R {config[refs][Reference]} -O {output} "
		"--known-sites {config[bqsrParam][dbSNP]} "
		"--known-sites {config[bqsrParam][indelFile]} "
		"--known-sites {config[bqsrParam][thousandGenomesSite]} "
		"-I {input.bam}"

rule bqsr_apply:
	input:
		bam="{dirname}/markdup_bams/{samplename}.la.md.bam",
		recal_file="{dirname}/happy_bqsr_bams/{samplename}.recal.csv"
	output:
		"{dirname}/happy_bqsr_bams/{samplename}.la.md.bqsr.bam"
	shell:
		"gatk ApplyBQSR --QUIET true --create-output-bam-index true -R {config[refs][Reference]} "
		"--bqsr-recal-file {input.recal_file} -I {input.bam} "
		"-O {output} "


rule bcftools_call:
	input:
		"{dirname}/happy_bqsr_bams/{samplename}.la.md.bqsr.bam"
	output:
		"{dirname}/bcfs/{samplename}.la.md.bqsr.bcftools.bcf"
	shell:
		"bcftools mpileup -d 1023 -q 9 "
		"-a 'AD,ADF,ADR,DP,SP,INFO/AD,INFO/ADF,INFO/ADR' "
		"-Ou -f {config[refs][Reference]} "
		"{input} "
		 "| bcftools call --threads {config[bcftoolsParams][bcfThreads]} "
		 "-g 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30 -f GP -m -Ob -o {output} "
		"&& bcftools index {output}"

rule bcftools_norm:
	input:
		"{dirname}/bcfs/{samplename}.la.md.bqsr.bcftools.bcf"
	output:
		"{dirname}/bcftools_norm/{samplename}.la.md.bqsr.bcftools.normed.bcf"
	shell:
		"bcftools norm -m -snps -c x -f {config[refs][Reference]} {input} "
		 "-Ob -o {output} && "
		 "bcftools index {output}"
