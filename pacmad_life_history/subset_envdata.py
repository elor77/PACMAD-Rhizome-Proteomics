import csv

infile = "/Users/eo235/Library/CloudStorage/OneDrive-Personal/Projects/BucklerLab/pacmad_life_history/bien_envData.txt"
outfile = "envdata_subset_13species.csv"

targets = {
    "Panicum virgatum", "Setaria viridis",
    "Zea mays", "Tripsacum dactyloides",
    "Sorghastrum nutans", "Sorghastrum fuscescens",
    "Andropogon gerardii", "Miscanthus sinensis",
    "Miscanthus sacchariflorus",
    "Sorghum bicolor", "Saccharum officinarum",
    "Themeda triandra", "Cymbopogon nardus",
    "Hyparrhenia bracteata", "Coix lacryma-jobi"
}

# Header has N columns but data rows have N+1 (extra row-number column)
with open(infile, "r") as f:
    header_line = f.readline().strip()
    headers = ["rownum"] + header_line.split("\t")

keep_cols = [
    "scientificName", "decimalLatitude", "decimalLongitude",
    "bio01_Annual_Mean_Temperature",
    "bio05_Max_Temperature_Warmest_Month",
    "bio06_Min_Temperature_Coldest_Month",
    "bio10_Mean_Temperature_Warmest_Quarter",
    "bio11_Mean_Temperature_Coldest_Quarter",
    "bio04_Temperature_Seasonality",
    "bio12_Annual_Precipitation",
    "Elevation_m"
]

counts = {}
with open(infile, "r") as fin, open(outfile, "w", newline="") as fout:
    fin.readline()  # skip header
    reader = csv.DictReader(fin, fieldnames=headers, delimiter="\t")
    writer = csv.DictWriter(fout, fieldnames=keep_cols)
    writer.writeheader()
    for row in reader:
        sp = row["scientificName"]
        if sp in targets:
            writer.writerow({c: row[c] for c in keep_cols})
            counts[sp] = counts.get(sp, 0) + 1

print(f"Total rows: {sum(counts.values())}")
for sp, n in sorted(counts.items()):
    print(f"  {sp}: {n}")

# Check for missing species and scan for genus-level near-misses
targets_genera = ["Panicum", "Setaria", "Zea", "Tripsacum", "Sorghastrum",
                  "Andropogon", "Miscanthus", "Sorghum", "Saccharum",
                  "Themeda", "Cymbopogon", "Hyparrhenia", "Coix"]
missing = [t for t in targets if t not in counts]
if missing:
    print(f"\nMissing species: {missing}")
    print("Scanning for genus matches...")
    genus_found = {}
    with open(infile, "r") as fin:
        fin.readline()
        reader = csv.DictReader(fin, fieldnames=headers, delimiter="\t")
        for row in reader:
            sp = row["scientificName"]
            for g in targets_genera:
                if sp and sp.startswith(g):
                    genus_found.setdefault(g, set()).add(sp)
    for g in targets_genera:
        if any(t.startswith(g) and t in missing for t in targets):
            names = sorted(genus_found.get(g, set()))[:10]
            print(f"  {g}: {names}")
else:
    print("\nAll species found!")