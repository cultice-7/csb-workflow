import pandas as pd
from pathlib import Path
import pyarrow as pa
import pyarrow.parquet as pq

# Directory containing the CSV files and output directories
dataset = "CSB"
input_folder = Path(f"data/edited/{dataset}")
output_folder = input_folder / "parquet"
output_folder.mkdir(exist_ok=True)

# Find all CSV files containing "_table"
csv_files = list(input_folder.glob("*supplement_*_table.csv"))


for csv_file in csv_files:
    print(f"Processing {csv_file.name}...")

    # Create output file name
    output_file = output_folder / (csv_file.stem + ".parquet")
    
    # Read input file
    df = pd.read_csv(csv_file)
    
    # Convert CSBID to string
    if dataset == "CSB":
        df['CSBID'] = df['CSBID'].astype("string")
    
    # Convert only float64 columns to float32
    float64_cols = df.select_dtypes(include=["float64"]).columns
    df[float64_cols] = df[float64_cols].astype("float32")
    
    # Save dataset to parquet format
    df.to_parquet(output_file, compression="zstd")

    print(f"Saved → {output_file.name}")

print("All files converted")


# Process the files using chunks
"""
for csv_file in csv_files:
    print(f"Processing {csv_file.name}...")

    output_file = output_folder / (csv_file.stem + ".parquet")
    
    writer = None
    
    # Read everything as string and detect true column nature
    for chunk in pd.read_csv(csv_file, chunksize=100_000, dtype=str, low_memory=False):
    
        # Convert only float64 columns to float32
        float64_cols = chunk.select_dtypes(include=["float64"]).columns
        chunk[float64_cols] = chunk[float64_cols].astype("float32")
    
        # Convert pandas DataFrame to Arrow table
        table = pa.Table.from_pandas(chunk, preserve_index=False)

        # Initialize writer once
        if writer is None:
            writer = pq.ParquetWriter(
                output_file,
                table.schema,
                compression="zstd"
            )

        writer.write_table(table)

    if writer:
        writer.close()

    print(f"Saved → {output_file.name}")

print("All files converted")
"""