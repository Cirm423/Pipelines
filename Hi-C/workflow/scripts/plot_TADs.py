#Basic TADs + loops plot taken from TADlib documentation.

from tadlib.visualize.heatmaps import *

vis = Triangle(snakemake.params["uri"], snakemake.params["coords"][0], int(snakemake.params["coords"][1]), int(snakemake.params["coords"][2]))
vis.matrix_plot()
vis.plot_loops(snakemake.input["loops"])
vis.plot_TAD(snakemake.input["tad"], linewidth=1.5)
vis.outfig(snakemake.output[0])