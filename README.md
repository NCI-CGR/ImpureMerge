# ImputeMerge: Merging TOPMed Imputation Server imputation data

## About
Due to the input sample size limit by TOPMed imputation server, input genomic data is processed in batches and merged later.
Goal of this pipeline is to take data imputed in batches and merge to produce a single file.

This pipeline filters for SNPs R2 > 0.25 (in atleast on one batch) and filters for R2 > 0.3 in merged file.
bcftools is used for merging.  

Note: Due to the file size limit, input data sample has been uploaded into this repository.
## Command line arguments

Copy config.yaml and Snakefile.sh to your directory, edit the config file with input and output directory as needed

```bash
module load python3
module load slurm
cd ImputeMerge
sbatch --partition=cgrq -o temp.stdout Snakefile.sh
```

## Config file details:
	* directory_in: /full/path/to/input/data_sample/input
	* directory_out: /full/path/to/output/data_sample/output
	* batches: number of batches separated by comma (Example: "1,2,3")
 	* rsq_1: rsq first filter cut-off (filters for rsq in each batch; default 0.25)
  	& rsq_2: rsq second filter cut-off (filters for merged average rsq ; default 0.3)

## Output Folders/Files

	1. info_Reformat: Reformated original info file to below format

	* Example reformated info file:  

	 |  CHROM | POS   |     ID       | REF | ALT |   AF    |    MAF  |    R2   | IMPUTED |
	 | -----  |-----  | ------------ | --- | --- | ------- | ------- | ------  | ------- |
	 | chr1   | 14675 | rs1339357485 | C   | A   | 0.00003 | 0.00003 | 0.44840 | IMPUTED |
	 | chr1   | 14766 | rs1420833025 | T   | G   | 0.00004 | 0.00004 | 0.46070 | IMPUTED |
	 | chr1   | 14808 | .            | A   | G   | 0.00002 | 0.00002 | 0.65686 | IMPUTED |
	 | chr1   | 14838 | rs1401618782 | C   | T   | 0.00001 | 0.00001 | 0.18849 | IMPUTED |


	2. data – Intermediate files used to merge the final data
	* R2>0.25 filtered variants text file used to filter vcf files before merging:

	3. VCF_filter - R2>0.25 filtered vcf files for each batch

	4. merged – Final merged vcf files for each chromosome
