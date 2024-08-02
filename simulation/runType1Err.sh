## To get type I error
## 1. we simulated a GWAS summmary statistics considering LD information.
## 2. we used the real x, y, cov matrix, and the simulated GWAS summary statistics to explore the type I error.
## 3. For each gene, we ran MTSG of 10,000 times.
## 4. For all the testable gene in the genome, the simulation achieved 10,000,000 times in total.
## 5. Be aware, the example_data only includes one gene. Our type I error simulation were done on all genes.

parallel -j 30 -q echo 'Rscript simulation/Type1Err_0226.R --x example_data/{}.select.rsid.vcf.ac.x_all --y example_data/{}.all.final.matrix.decomp   --cov example_data/{}.select.rsid.vcf.ac.cov.RData --nsims 10000 --dir_out Type1Err/ --verbose TRUE' ::: $(ls /nobackup/cgg/yany14/MSG/test_rtensor/test_final/ | grep decomp| perl -lane 'print $1 if /(ENSG\d+)/' | sort | uniq ) > stimulation/Type1Err.cmd

## Run the cmd file using computing cluster. Here is an excutable command line example.
bash Type1Err.cmd
