import csv

targets_genera = ["Panicum", "Setaria", "Zea", "Tripsacum", "Sorghastrum",
                  "Andropogon", "Miscanthus", "Sorghum", "Saccharum",
                  "Themeda", "Cymbopogon", "Hyparrhenia", "Coix"]

infile = "/Users/eo235/Library/CloudStorage/OneDrive-Personal/Projects/BucklerLab/pacmad_life_history/bien_envData.txt"  # adjust filename

found = {}
with open(infile, "r") as fin:
    reader = csv.DictReader(fin, delimiter="\t")
    for row in reader:
        sp = row.get("scientificName", "")
        for g in targets_genera:
            if sp.startswith(g):
                found.setdefault(g, set()).add(sp)

for g in targets_genera:
    names = found.get(g, set())
    print(f"{g}: {sorted(names)[:10]}")