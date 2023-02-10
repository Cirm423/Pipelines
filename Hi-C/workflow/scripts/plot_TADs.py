#Basic TADs + loops plot taken from TADlib documentation.

from tadlib.visualize.heatmaps import *

#Process the extra params
extra_dict = {"--dpi":300,"--vmin":None,"--vmax":None}

extra = snakemake.params["extra"].split(" ")

for key in extra_dict.keys():
    if key in extra:
        print(key)
        value = float(extra[extra.index(key) + 1])
        extra_dict[key] = value

vis = Triangle(snakemake.params["uri"], snakemake.params["coords"][0], int(snakemake.params["coords"][1]), int(snakemake.params["coords"][2]))
vis.matrix_plot(vmin=extra_dict['--vmin'],vmax=extra_dict['--vmax'])
vis.plot_loops(snakemake.input["loops"])
vis.plot_TAD(snakemake.input["tad"], linewidth=1.5)
vis.outfig(snakemake.output[0],dpi=extra_dict['--dpi'])