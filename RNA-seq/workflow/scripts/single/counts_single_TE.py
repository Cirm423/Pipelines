from snakemake.shell import shell
import os

for file in snakemake.input.treat:
    file_name = os.path.basename(os.path.dirname(file))
    control_name = os.path.basename(os.path.dirname(snakemake.input.control[0]))
    shell(r"""paste {file} {snakemake.input.control[0]}  | awk 'BEGIN{{OFS="\t"}} {{if($15>0) {{print $1,$2,$3,$4,$5,$6,$7,$15,$8,$16,$7/$15}} else {{print $1,$2,$3,$4,$5,$6,$7,"1",$8,$16,$7/1}}}}' > results/TE_single/counts.temp 2>{snakemake.log[0]}""")
    shell("Rscript workflow/scripts/pois-test_uniq-TE-individual.R {snakemake.params[0]} 2>>{snakemake.log[0]}")
    shell(r"""echo -e "Chr\tStart\tEnd\tRegion_ID\tMappability\tStrand\t{file_name}_NormalizedCounts\t{control_name}_NormalizedCounts\t{file_name}_rpkm\t{control_name}_rpkm\tFoldChange\tp-value\tFDR" > results/TE_single/{file_name}_vs_{control_name}.tsv 2>>{snakemake.log[0]}""")
    shell("cat results/TE_single/counts_diff.temp >> results/TE_single/{file_name}_vs_{control_name}.tsv 2>>{snakemake.log[0]}")