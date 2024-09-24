from getpass import getuser
from pathlib import Path
import pandas as pd
import numpy as np
import os
import socket
import re
import glob
import shutil
import gzip
import csv

CHRS=[str(item+1) for item in range(22)]+['X']

RWD = os.getcwd()

configfile: "config/config.yaml"

out_dir=config["directory_out"]

INPATH = config["directory_in"]
OUTPATH = config["directory_out"]
BATCHES = config["batches"].split(",")

ruleorder: INFO_Reformat > join > extend > ID_INFO > filter_bcftools > prepare_file > combine



rule all:
    input:
        expand(OUTPATH + "/info_Reformat/batch{BATCH}/chr{CHR}.Reformat_info.gz",CHR=CHRS,BATCH=BATCHES),
        expand(OUTPATH + "/data/chr{CHR}.Reformat_info.gz",CHR=CHRS),
        expand(OUTPATH + "/data/chr{CHR}.info.txt.gz",CHR=CHRS),
        expand(OUTPATH + "/data/chr{CHR}.info.pick.txt.gz",CHR=CHRS),
        expand(OUTPATH + "/data/chr{CHR}.info.pick.True.txt",CHR=CHRS),
        expand(OUTPATH + "/data/chr{CHR}.info.pick.rsq03.txt",CHR=CHRS),
        expand(OUTPATH + "/data/chr{CHR}.info.pick.True.first3columns.txt",CHR=CHRS),
        expand(OUTPATH + "/VCF_filter/batch{BATCH}/chr{CHR}.dose.filter.vcf.gz",CHR=CHRS,BATCH=BATCHES),
        expand(OUTPATH + "/data/chr{CHR}.dose.filter.vcf.gz.csv",CHR=CHRS),
        expand(OUTPATH + "/merged/chr{CHR}.dose.filter.vcf.gz",CHR=CHRS)


rule INFO_Reformat:
    input:
        infoFile = INPATH + "/" + "batch{BATCH}/chr{CHR}.info.gz"
    output:
        outFile = OUTPATH + "/info_Reformat/batch{BATCH}/chr{CHR}.Reformat_info.gz"
    #priority: 59
    run:
        with gzip.open(input.infoFile,'rb') as vcffile:
            vcfheader=vcffile.readline()
            vcfheaderDecode=vcfheader.decode('ascii')
            while vcfheaderDecode.startswith("##"):
                vcfheader=vcffile.readline()
                vcfheaderDecode=vcfheader.decode('ascii')
            vcfheaderSplit=vcfheaderDecode.rstrip("\n\r").split("\t")
            IDs=[]
            for line in vcffile:
                lineDecode=line.decode('ascii')
                if not lineDecode.startswith("#"):
                    lineSplit=lineDecode.rstrip("\n\r").split("\t")
                    IDs.append(lineSplit)
            df = pd.DataFrame(data=IDs, columns=['CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO'])
            df['AF'] = df['INFO'].str.split(';', expand=True)[0].map(lambda x: x.lstrip('AF='))
            df['MAF'] = df['INFO'].str.split(';', expand=True)[1].map(lambda x: x.lstrip('MAF='))
            df['AVG_CS'] = df['INFO'].str.split(';', expand=True)[2].map(lambda x: x.lstrip('AVG_CS='))
            df['R2'] = df['INFO'].str.split(';', expand=True)[3].map(lambda x: x.lstrip('R2='))
            df['IMPUTED'] = df['INFO'].str.split(';').str[-1]
            df.to_csv(output.outFile, sep='\t', columns = ['CHROM','POS','ID','REF','ALT','AF','MAF','R2','IMPUTED'] , index=False, compression='gzip')



rule join:
#    input:
#        OUTPATH + "/info_Reformat/batch{BATCH}/chr{CHR}.Reformat_info.gz"
    output:
        OUTPATH + "/data/chr{CHR}.Reformat_info.gz"
    threads:
          5
    #priority: 58
    run:
        root = OUTPATH + "/info_Reformat"
        CONCAT_DIR = OUTPATH + "/data/"
        files = pd.DataFrame([file for file in sorted(glob.glob(root + "/*/*info.gz"))], columns=["fullpath"])
        files_split = files['fullpath'].str.rsplit("/", n=1, expand=True).rename(columns={0: 'path', 1:'filename'})
        files = files.join(files_split)
        for f in files['filename'].unique():
          paths = files[files['filename'] == f]['fullpath']
          dfs = [pd.read_table(path,compression="gzip") for path in paths]
          concat_df = pd.concat(dfs,axis=1)
          cols = pd.Series(concat_df.columns)
          dup_count = cols.value_counts()
          for dup in cols[cols.duplicated()].unique():
            cols[cols[cols == dup].index.values.tolist()] = [dup +'_b'+ str(i) for i in range(1, dup_count[dup]+1)]
          concat_df.columns = cols
          concat_df.to_csv(CONCAT_DIR + f, sep='\t', index=False)
#          fileList = os.listdir(CONCAT_DIR)
#          for ii in fileList:
#              newName = ii.replace('Reformat_info','info_pre.txt')
#              if newName != ii:
#                  os.rename(CONCAT_DIR+ii,CONCAT_DIR+newName)

rule extend:
    input:
        ipath = OUTPATH + "/data/chr{CHR}.Reformat_info.gz"
    output:
        opath = OUTPATH + "/data/chr{CHR}.info.txt.gz"
    threads:
        3
    #priority: 57
    run:
        df = pd.read_table(input.ipath, compression='gzip')
        df["AF_avg"] = df[df.columns[df.columns.str.startswith('AF')]].mean(axis=1)
        df["MAF_avg"] = df[df.columns[df.columns.str.startswith('MAF')]].mean(axis=1)
        df["R2_avg"] = df[df.columns[df.columns.str.startswith('R2')]].mean(axis=1)
        df.to_csv(output.opath, sep='\t', index=False, compression='gzip')

rule pick:
    input:
        inp = OUTPATH + "/data/chr{CHR}.info.txt.gz"
    output:
        one = OUTPATH + "/data/chr{CHR}.info.pick.txt.gz",
        two = OUTPATH + "/data/chr{CHR}.info.pick.True.txt",
    params:
        rsq = config["rsq_1"]
    #priority: 56
    run:
        df = pd.read_table(input.inp, compression='gzip')
        r2_filter_col = [col for col in df if col.startswith('R2_b')]
        r2_col = df[r2_filter_col]
        r2_array = np.where(df[df.columns[df.columns.str.startswith('R2_b')]] >= params.rsq, True, False)
        r2_array_df = pd.DataFrame(r2_array)
        combined_df = pd.concat([df, r2_array_df], axis=1)
        combined_df.to_csv(output.one, sep='\t', index=False, compression='gzip')
        pick_df = combined_df[combined_df.apply(lambda row: row.astype(str).str.contains('True').any(), axis=1)]
        pick_df.to_csv(output.two, sep='\t', index=False)


rule pick_rsq03:
    input:
        inp = OUTPATH + "/data/chr{CHR}.info.pick.True.txt"
    output:
        one = OUTPATH + "/data/chr{CHR}.info.pick.rsq03.txt"
    params:
        rsq = config["rsq_2"]
    #priority: 56
    run:
        df = pd.read_table(input.inp)
        df['Avg_Rsq_gt03'] = np.where(df['R2_avg'] >= params.rsq, True, False)
        out_df = df[df['Avg_Rsq_gt03'].astype(str).str.contains("True")]
        out_df.to_csv(output.one, sep='\t', index=False)


rule ID_INFO:
    input:
        OUTPATH + "/data/chr{CHR}.info.pick.rsq03.txt"
    output:
        OUTPATH + "/data/chr{CHR}.info.pick.True.first3columns.txt"
    #priority: 55
    shell:
        """
        cat {input} | cut -f1-3 | tail -n +2 > {output}
        """

rule filter_bcftools:
    input:
        pick=OUTPATH + "/data/chr{CHR}.info.pick.True.first3columns.txt",
        dose=INPATH + "/batch{BATCH}/chr{CHR}.dose.vcf.gz"
    output:
        out=OUTPATH + "/VCF_filter/batch{BATCH}/chr{CHR}.dose.filter.vcf.gz",
        tbi=OUTPATH + "/VCF_filter/batch{BATCH}/chr{CHR}.dose.filter.vcf.gz.tbi"
    threads:
          5
    #priority: 54
    shell:
        """
        module load bcftools;
        bcftools view {input.dose} -Oz --targets-file {input.pick} -o {output.out}
        module load tabix; tabix -f -p vcf {output.out}
        """
'''
rule rename_chr23:
    input:
        vcf = OUTPATH + "/VCF_filter/batch{BATCH}/chrX.dose.filter.vcf.gz",
    #    tbi = OUTPATH + "/VCF_filter/batch{BATCH}/chrX.dose.filter.vcf.gz.tbi",
    output:
        vcf = OUTPATH + "/VCF_filter/batch{BATCH}/chr23.dose.filter.vcf.gz",
    #    tbi = OUTPATH + "/VCF_filter/batch{BATCH}/chr23.dose.filter.vcf.gz.tbi",
    run:
        indir = OUTPATH
        for f in Path(indir + '/VCF_filter').rglob('chrX.dose.filter.vcf.gz'):
            os.rename(f, os.path.join(f.parent,'chr23.dose.filter.vcf.gz'))
'''
rule prepare_file:
    input:
        OUTPATH + "/VCF_filter/batch1/chr1.dose.filter.vcf.gz"
    output:
        one = OUTPATH + "/data/chr{CHR}.dose.filter.vcf.gz.csv",
#    priority: 0
    run:
        indir = OUTPATH
        outdir = OUTPATH + "/data/"
        file_paths = {}
        for file_path in Path(indir + '/VCF_filter').rglob('*.dose.filter.vcf.gz'):
          file_name = os.path.basename(file_path)
          if file_name not in file_paths:
            file_paths[file_name] = []
          file_paths[file_name].append(file_path)

        for index, (file_name, paths) in enumerate(file_paths.items()):
          csv_file_name = f'{outdir}{file_name}.csv'
        #  csv_file_name = f'{output.one}'
          with open(csv_file_name, mode='w', newline='') as csvfile:
            csv_writer = csv.writer(csvfile)
            for word in paths:
              csv_writer.writerow([word])


rule combine:
    input:
        OUTPATH + "/data/chr{CHR}.dose.filter.vcf.gz.csv",
    output:
        OUTPATH + "/merged/chr{CHR}.dose.filter.vcf.gz"
    threads:
       10
#    priority: 0
    shell:
        """
        module load bcftools;
        bcftools merge --merge none --info-rules R2:avg,AF:avg,MAF:avg -l {input} -Oz -o {output}
        """
