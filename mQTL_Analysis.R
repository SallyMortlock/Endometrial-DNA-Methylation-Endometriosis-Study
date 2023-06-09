################################################################
# This script runs matrixEQTL to identify endometrial mQTLs.

# Author: Sally Mortlock
# July 2022
################################################################

library("MatrixEQTL")

## Location of the package with the data files.
base.dir = "/scratch/Methylation";
# base.dir = '.';

## Settings

# Linear model to use, modelANOVA, modelLINEAR, or modelLINEAR_CROSS
useModel = modelLINEAR; # modelANOVA, modelLINEAR, or modelLINEAR_CROSS

# Genotype file name
SNP_file_name = "/scratch/Methylation/Euro_GenotypeData_July2022_mQTL.csv";
snps_location_file_name = "/scratch/Methylation/Euro_SNPPosData_July2022_mQTL.csv";

# Gene expression file name
expression_file_name = "/scratch/Methylation/M_batch1_batch2_July2022_mQTL.csv";
gene_location_file_name = "/scratch/Methylation/M_batch1_batch2_July2022_mQTL_pos.csv";

# Covariates file name
# Set to character() for no covariates
covariates_file_name = "/scratch/Methylation/M_batch1_batch2_July2022_mQTL_cov.csv";

# Output file name
output_file_name_cis = tempfile();
output_file_name_tra = tempfile();

# Only associations significant at this level will be saved
pvOutputThreshold_cis = 1;
pvOutputThreshold_tra = 0;

# Error covariance matrix
# Set to numeric() for identity.
errorCovariance = numeric();
# errorCovariance = read.table("Sample_Data/errorCovariance.txt");

# Distance for local gene-SNP pairs
cisDist = 1e6;

## Load genotype data

snps = SlicedData$new();
snps$fileDelimiter = ",";      # the TAB character
snps$fileOmitCharacters = "NA"; # denote missing values;
snps$fileSkipRows = 1;          # one row of column labels
snps$fileSkipColumns = 1;       # one column of row labels
snps$fileSliceSize = 2000;      # read file in slices of 2,000 rows
snps$LoadFile(SNP_file_name);

## Load gene expression data

gene = SlicedData$new();
gene$fileDelimiter = ",";      # the TAB character
gene$fileOmitCharacters = "NA"; # denote missing values;
gene$fileSkipRows = 1;          # one row of column labels
gene$fileSkipColumns = 1;       # one column of row labels
gene$fileSliceSize = 2000;      # read file in slices of 2,000 rows
gene$LoadFile(expression_file_name);

## Load covariates

cvrt = SlicedData$new();
cvrt$fileDelimiter = ",";      # the TAB character
cvrt$fileOmitCharacters = "NA"; # denote missing values;
cvrt$fileSkipRows = 1;          # one row of column labels
cvrt$fileSkipColumns = 1;       # one column of row labels
if(length(covariates_file_name)>0) {
  cvrt$LoadFile(covariates_file_name);
}

## Run the analysis
snpspos = read.csv(snps_location_file_name, header = TRUE, stringsAsFactors = FALSE);
genepos = read.csv(gene_location_file_name, header = TRUE, stringsAsFactors = FALSE);

me = Matrix_eQTL_main(
  snps = snps, 
  gene = gene, 
  cvrt = cvrt,
  output_file_name     = output_file_name_tra,
  pvOutputThreshold     = pvOutputThreshold_tra,
  useModel = useModel, 
  errorCovariance = errorCovariance, 
  verbose = TRUE, 
  output_file_name.cis = output_file_name_cis,
  pvOutputThreshold.cis = pvOutputThreshold_cis,
  snpspos = snpspos, 
  genepos = genepos,
  cisDist = cisDist,
  pvalue.hist = "qqplot",
  min.pv.by.genesnp = FALSE,
  noFDRsaveMemory = FALSE);

unlink(output_file_name_tra);
unlink(output_file_name_cis);

## Results:

cat('Analysis done in: ', me$time.in.sec, ' seconds', '\n');
cat('Detected local eQTLs:', '\n');
#show(me$cis$eqtls)
dim(me$cis$eqtls)
cat('Detected distant eQTLs:', '\n');
#show(me$trans$eqtls)
dim(me$trans$eqtls)

mqtl <- as.data.frame(me$cis$eqtls)
saveRDS(mqtl, "R01_Euro_cis_mQTL_Endometrium_July2022_1mb_Stage_Endo_SVA.rds")
