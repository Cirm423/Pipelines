from snakemake.shell import shell
import os

for file in snakemake.input.treat:
    file_name = os.path.basename(os.path.splitext(file)[0])
    control_name = os.path.basename(os.path.splitext(snakemake.input.control[0])[0])
    shell(r"""paste {file} {snakemake.input.control[0]}  | awk -F"\t" 'BEGIN{{OFS="\t"}} {{print $1,$2,$3,$4,$5,$12,$6,$13,$7,$14,$5/$12}}' > results/single/counts.temp 2>{snakemake.log[0]}""")
    shell("Rscript workflow/scripts/pois-test_RSEM.R {snakemake.params[0]} 2>>{snakemake.log[0]}")
    shell(r"""echo -e "GeneName\tGeneID\tlength\teffective_length\t{file_name}_NormalizedCounts\t{control_name}_NormalizedCounts\t{file_name}_TPM\t{control_name}_TPM\t{file_name}_FPKM\t{control_name}_FPKM\tFoldChange\tp-value\tFDR" > results/single/{file_name}_vs_{control_name}.tsv 2>>{snakemake.log[0]}""")
    shell("cat results/single/counts_diff.temp >> results/single/{file_name}_vs_{control_name}.tsv 2>>{snakemake.log[0]}")