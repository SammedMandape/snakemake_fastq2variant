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

		#print(SAMPLES, LANE)
#SAMPLES=

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
         expand("{dirname}/deepvariant_vcf/{samplename}.output.vcf.gz", samplename=SAMPLES_UNIQ_L,dirname=DIR_UNIQ_L)

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
			seqID=$(echo $header | cut -d':' -f 1 | cut -c 2-)
			id=$(echo $header | cut -d':' -f 2-4,10 --output-delimiter='_')
			sm=$(echo {wildcards.samplename} | cut -d'_' -f 1)
			pu=$(echo $header | cut -d':' -f 3,4 --output-delimiter='_')
			rg=$(echo "@RG\tID:$seqID"_"$id\tSM:$sm\tLB:$sm\tPL:ILLUMINA")
			bwa mem -R $(echo "@RG\\tID:$sm"_"$seqID"_"$id\\tSM:$sm\\tLB:$sm\\tPL:ILLUMINA\\tPU:$pu") \
			-t {config[bwamemParams][bwaThreads]} \
			-M {config[refs][Reference]} \
			{read1} {read2} | \
			samtools view -bh - | \
			samtools sort -m 50G - > {output} && samtools index {output} 2> {log}
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
		shell("java -jar /usr/local/src/picard.jar MergeSamFiles CREATE_INDEX=true O={output} {inputs_all}")


# TODO: take deepvariant bin version from config.yaml file.
rule deepvariant:
    input:
         "{dirname}/merged_bams/{samplename}.bam"
    output:
          vcf="{dirname}/deepvariant_vcf/{samplename}.output.vcf.gz",
          gvcf="{dirname}/deepvariant_vcf/{samplename}.output.g.vcf.gz"
    run:
        shell(
            """
            singularity run -B /usr/lib/locale,/mnt/ docker://google/deepvariant:"1.2.0" \
            /opt/deepvariant/bin/run_deepvariant --model_type=WGS \
            --ref={config[refs][Reference]} \
            --reads={input} \
            --output_vcf={output.vcf} \
            --output_gvcf={output.gvcf} \
            --num_shards=128
            """
        )
