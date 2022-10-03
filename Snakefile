rule all:  
	input: "output/regressionAnalysis"

# STAR
rule index:
	input:
             fa = "inputs/STAR_input/Homo_sapiens.GRCh38.dna.chromosome.10.fa",
	     gtf = "inputs/STAR_input/gencode.v29.annotation_chr10.gtf"
	output: directory('output/STAR_output/indexes')
	shell: "mkdir {output} && STAR --runMode genomeGenerate --genomeDir {output} --genomeFastaFiles {input.fa} --sjdbGTFfile {input.gtf} --sjdbOverhang 50 --outFileNamePrefix indexes/chr10"

rule alignments:
	input: 
		output_index= 'output/STAR_output/indexes',
		een = 'inputs/STAR_input/A549_25_3chr10_2.fastq.gz',
		twee = 'inputs/STAR_input/A549_25_3chr10_1.fastq.gz'
	output: directory('output/STAR_output/alignments')
        shell: ' mkdir {output} && STAR --genomeDir output/STAR_output/indexes --readFilesIn {input.een} {input.twee} --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --quantMode GeneCounts --outFileNamePrefix output/STAR_output/alignments/A549_0_2'


rule normalisation:
	input:"output/STAR_output/alignments/A549_0_2ReadsPerGene.out.tab"
	output: "output/STAR_output/normalised_counts_A549_matrix.txt"
        shell: "python3 src/normalizing.py {input} {output}"

#the files had version numbers these have been removed: final output of star
rule remove_version:
	input: "output/STAR_output/alignments/normalised_counts_A549_matrix.txt"
	output:"output/STAR_output/normalized_values.txt"
	shell: "python3 src/version.py {input} {output}"

#Executes the following 6 GATK commands and collects the Clinvar files in the meantime

#Assigns all the reads in a file to a single new read-group
rule AddOrReplaceReadGroups:
	input: 'inputs/GATK_inputs/A549_0_2Aligned.sortedByCoord.out.bam'
	output: 'output/gatk_outputs/A549_0_2Aligned.sortedByCoord.out.bam'
	shell: 'gatk AddOrReplaceReadGroups I={input} O={output} RGID=4 RGLB=lib1 RGPL=ILLUMINA RGPU=unit1 RGSM=20'

#Splits reads that contains Ns in their cigar string
rule SplitNCigarReads:
	input:
		bam = "output/gatk_outputs/A549_0_2Aligned.sortedByCoord.out.bam",
		fasta = "inputs/STAR_input/Homo_sapiens.GRCh38.dna.chromosome.10.fa"
	output: "output/gatk_outputs/A549_0_2Aligned-SNCR.sortedByCoord.out.bam"
	shell: "gatk SplitNCigarReads -R {input.fasta} -I {input.bam} -O {output}"

#Retrieves the .vcf.gz file from the Clinvar database
rule get_clinvar:
	shell: 'wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/archive_2.0/2022/clinvar_20220213.vcf.gz'

#Retrieves the .vcf.gz.tbi file from the Clinvar database
rule get_tbi:
	shell: 'wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/archive_2.0/2022/clinvar_20220213.vcf.gz.tbi'

#'10' will be replaced with 'chr10' for comparing the 2 files in BaseRecalibrator
rule convert:
	shell: 'echo "10 chr10" >> chr_name_conv.txt'

rule annotation:
	shell: 'bcftools annotate -rename-chrs chr_name_conv.txt clinvar_20220213.vcf.gz | bgzip > clinvar_20220213_name.vcf.gz'

#The Clinvar file will be indexed
rule indexing:
	shell: 'gatk IndexFeatureFile -I clinvar_20220213_name.vcf.gz'

#Will create a .vcf file with only chr10 in it
rule replacing:
	shell: 'bcftools view -r chr10 clinvar_20220213_name.vcf.gz -o chr10.vcf'

#The file will be zipped
rule zip_file:
	shell: 'bgzip -c chr10.vcf > chr10.vcf.gz'

#The zipped file will be indexed
rule index_zip:
	shell: 'gatk IndexFeatureFile -I chr10.vcf.gz'

#Generates recalibration table for BQSR (next step) with the Clinvar file
rule BaseRecalibrator:
	input:
		bam = 'output/gatk_outputs/A549_0_2Aligned-SNCR.sortedByCoord.out.bam',
		fasta = 'inputs/STAR_input/Homo_sapiens.GRCh38.dna.chromosome.10.fa'
	output: 'output/gatk_outputs/recalibrator_data.table'
	shell: 'gatk BaseRecalibrator -I {input.bam} -R {input.fasta} --known-sites chr10.vcf.gz -O {output}'

#It recalibrates the base qualities of the input reads based on the recalibration table (made by the BaseRecalibrator)
rule ApplyBQSR:
	input:
		bam = 'output/gatk_outputs/A549_0_2Aligned-SNCR.sortedByCoord.out.bam',
		fasta = 'inputs/STAR_input/Homo_sapiens.GRCh38.dna.chromosome.10.fa',
		recalibrated = 'output/gatk_outputs/recalibrator_data.table'
	output: 'output/gatk_outputs/A549_0_2Aligned-recalibrated.sortedByCoord.out.bam'
	shell: 'gatk ApplyBQSR -I {input.bam} -R {input.fasta} -bqsr {input.recalibrated} -O {output}'

#It will call the potential variant sited per sample
rule HaplotypeCaller:
	input:
		fasta = 'inputs/STAR_input/Homo_sapiens.GRCh38.dna.chromosome.10.fa',
		recalibrated = 'output/gatk_outputs/A549_0_2Aligned-recalibrated.sortedByCoord.out.bam'
	output: 'output/gatk_outputs/Haplo.g.vcf.gz'
	shell: 'gatk --java-options "-Xmx4g" HaplotypeCaller -R {input.fasta} -I {input.recalibrated} -O {output} -ERC GVCF'

#Performs joint genotyping
rule GenotypeGVCFs:
	input:
		fasta = 'inputs/STAR_input/Homo_sapiens.GRCh38.dna.chromosome.10.fa',
		haplo = 'output/gatk_outputs/Haplo.g.vcf.gz'
	output: 'output/gatk_outputs/final_output_genotypegvcfs.vcf.gz'
	shell: 'gatk GenotypeGVCFs -R {input.fasta} -V {input.haplo} -O {output}'

#SNP's
rule snp_sel:
	input: een="inputs/SNPS_input/SNPs_chr10.tsv",
	       twee="output/gatk_outputs/final_output_genotypegvcfs.vcf.gz"
	output: "output/SNP_output/snpselnew.txt"
	shell: "python3 src/snps_sel.py {input.een} {input.twee} {output}"
# regression analysis
rule regressie: 
	input: sel="output/SNP_output/snpselnew.txt",
		counts="output/STAR_output/normalized_values.txt"
	output:	directory("output/regressionAnalysis")
	shell: "mkdir {output} && python3 src/regressie_analyse_final.py {input.counts} {input.sel} {output}"
