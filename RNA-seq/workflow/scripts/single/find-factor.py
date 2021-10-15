for file in snakemake.input:
    res=[]
    with open(file,"r") as f:
        for line in f:
            if line.strip().startswith("Number of input reads"):
                res.append(line.split("|")[1].strip())
                continue
else:
    final = min(res)
    with open(snakemake.output[0],"w") as f:
        f.write(final)