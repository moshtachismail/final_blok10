#rule all: moet nog gemaakt worden 

# STAR
rule index:
	input:
             fa = "/Users/mushtaaqismail/Documents/course_9_10/STAR_in_snakemake/STAR_inputs/Homo_sapiens.GRCh38.dna.chromosome.10.fa",
	     gtf = "/Users/mushtaaqismail/Documents/course_9_10/STAR_in_snakemake/input_files/STAR_inputs/gencode.v29.annotation_chr10.gtf"
	output: directory('output/STAR_output/indexes')
	shell: "mkdir {output} && STAR --runMode genomeGenerate --genomeDir {output} --genomeFastaFiles {input.fa} --sjdbGTFfile {input.gtf} --sjdbOverhang 50 --outFileNamePrefix indexes/chr10"

rule alignments:
	input: 
		output_index= '/Users/mushtaaqismail/Documents/course_9_10/STAR_in_snakemake/indexes/chr10',
		een = '/Users/mushtaaqismail/Documents/course_9_10/STAR_in_snakemake/A549_25_3chr10_2.fastq.gz',
		twee = '/Users/mushtaaqismail/Documents/course_9_10/STAR_in_snakemake/A549_25_3chr10_1.fastq.gz'
	output: directory('output/STAR_output/alignments')
 	shell: ' mkdir {output} && STAR --genomeDir /Users/mushtaaqismail/Documents/course_9_10/STAR_in_snakemake/output/STAR_output/indexes --readFilesIn {input.een} {input.twee} --readFilesCommand gunzip -c --outSAMtype BAM SortedByCoordinate --quantMode GeneCounts --outFileNamePrefix output/STAR_output/alignments/A549_0_2'


rule normalisation:
	input:"/Users/mushtaaqismail/Documents/course_9_10/STAR_in_snakemake/output/STAR_output/alignments/A549_0_2ReadsPerGene.out.tab"
	output: "/Users/mushtaaqismail/Documents/course_9_10/STAR_in_snakemake/output/STAR_output/alignments/normalised_counts_A549_matrix.txt"
        shell: "python3 /Users/mushtaaqismail/Documents/course_9_10/STAR_in_snakemake/src/normalizing.py {input} {output}"

#the files had version numbers these have been removed: final output of star
rule remove_version:
	input: "/Users/mushtaaqismail/Documents/course_9_10/STAR_in_snakemake/output/STAR_output/alignments/normalised_counts_A549_matrix.txt"
	output:"/Users/mushtaaqismail/Documents/course_9_10/STAR_in_snakemake/output/STAR_output/normalized_values.txt"
	shell: "python3 /Users/mushtaaqismail/Documents/course_9_10/STAR_in_snakemake/src/version.py {input} {output}"

#GATK Rule moet hier!


#start getting SNP's
rule get_snps:
	input: "/Users/mushtaaqismail/Documents/course_9_10/STAR_in_snakemake/input_files/SNPS_input/test_deelchr10.txt"
	output: file= "/Users/mushtaaqismail/Documents/course_9_10/STAR_in_snakemake/output/SNP_output/get_snps/snps.tsv",
		makedir= directory('output/SNP_output/get_snps')
	shell: "mkdir {output.makedir} && python3 /Users/mushtaaqismail/Documents/course_9_10/STAR_in_snakemake/src/SNP_data.py {input} {output.file}"

#SNP selectie rule hier!
#regression analysis rule!
