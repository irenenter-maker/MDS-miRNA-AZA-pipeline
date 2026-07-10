import pandas as pd
import sys
from pathlib import Path
from functools import reduce

if len(sys.argv) != 3:
    print("Usage: python count_combination.py <path_to_directory> <output_file_name>")
    sys.exit(1)

path = Path(sys.argv[1])
output_file_name = sys.argv[2]

if not path.exists():
    print(f"Error: Path '{path}' does not exist")
    sys.exit(1)

if not path.is_dir():
    print(f"Error: '{path}' is not a directory")
    sys.exit(1)

files = list(path.glob("*_miRNA_counts.txt"))

if not files:
    print(f"Warning: No files matching '*_miRNA_counts.txt' found in {path}")
    sys.exit(1)

dataframes = []
for file in files:
    df = pd.read_csv(
        file, 
        sep="\t",
        header=None,
        names=["miRNA", file.name.split("_miRNA_counts.txt")[0]]
    )
    dataframes.append(df)

merged_df = reduce(lambda left, right: pd.merge(left, right, on='miRNA', how='outer'), dataframes)

merged_df = merged_df.fillna(0)
merged_df.set_index('miRNA', inplace=True)

merged_df.to_csv(path / output_file_name, sep="\t")