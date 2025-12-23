#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import polars as pl
import pandas as pd
import os
from pathlib import Path
from typing import List
import math


# # Split Files

# In[ ]:


def split_gwas_file(
    input_file: str,
    output_dir: str = "gwas_chunks",
    target_size_mb: int = 100
) -> List[str]:
    """
    Split GWAS file using pandas with automatic compression detection.
    Handles both .tsv.gz and .tsv files.
    """
    
    # Create output directory
    Path(output_dir).mkdir(exist_ok=True)
    
    # Get base filename and detect compression
    input_path = Path(input_file)
    is_gzipped = input_path.suffix == '.gz'
    
    if is_gzipped:
        base_name = input_path.stem.replace('.tsv', '')
    else:
        base_name = input_path.stem
    
    print(f"Reading file: {input_file}")
    print(f"Compression detected: {'gzip' if is_gzipped else 'none'}")
    
    # Read the entire file with automatic compression detection
    compression = 'gzip' if is_gzipped else None
    gwas_df = pd.read_csv(input_file, sep='\t', compression=compression)
    total_rows = len(gwas_df)
    
    print(f"Total rows: {total_rows:,}")
    print(f"Columns: {list(gwas_df.columns)}")
    print(f"Memory usage: {gwas_df.memory_usage(deep=True).sum() / 1024**2:.1f} MB")
    
    # Estimate compressed size per row by testing a sample
    sample_size = min(10000, total_rows // 10)
    sample_df = gwas_df.head(sample_size)
    
    # Write sample to get actual compressed size
    test_file = os.path.join(output_dir, "test_sample.tsv.gz")
    sample_df.to_csv(test_file, sep='\t', index=False, compression='gzip')
    
    sample_size_mb = os.path.getsize(test_file) / (1024 * 1024)
    os.remove(test_file)  # Clean up
    
    # Calculate actual compression ratio
    bytes_per_row_compressed = (sample_size_mb * 1024 * 1024) / sample_size
    
    print(f"Sample: {sample_size:,} rows = {sample_size_mb:.2f} MB compressed")
    print(f"Estimated: {bytes_per_row_compressed:.1f} bytes per row (compressed)")
    
    # Calculate chunk size based on actual compression
    target_bytes = target_size_mb * 1024 * 1024
    chunk_size = int(target_bytes / bytes_per_row_compressed)
    
    # Add safety margin
    chunk_size = int(chunk_size * 0.9)  # 10% safety margin
    
    num_chunks = math.ceil(total_rows / chunk_size)
    
    print(f"Adjusted chunk size: {chunk_size:,} rows")
    print(f"Expected chunks: {num_chunks}")
    print("\nProcessing chunks...")
    
    chunk_files = []
    
    for chunk_num in range(num_chunks):
        start_idx = chunk_num * chunk_size
        end_idx = min((chunk_num + 1) * chunk_size, total_rows)
        
        # Extract chunk
        chunk = gwas_df.iloc[start_idx:end_idx]
        
        if len(chunk) == 0:
            break
        
        # Create output filename (always .tsv.gz for consistency)
        chunk_filename = f"{base_name}_chunk_{chunk_num + 1:03d}.tsv.gz"
        chunk_path = os.path.join(output_dir, chunk_filename)
        
        # Write compressed chunk
        chunk.to_csv(chunk_path, sep='\t', index=False, compression='gzip')
       
        chunk_files.append(chunk_path)
        
        # Check actual file size
        actual_size_mb = os.path.getsize(chunk_path) / (1024 * 1024)
        progress = end_idx / total_rows * 100
        
        print(f"Chunk {chunk_num + 1}: {len(chunk):,} rows, {actual_size_mb:.1f} MB ({progress:.1f}% complete)")
    
    print(f"\n✓ Created {len(chunk_files)} chunks in '{output_dir}'")
    
    # Summary of file sizes
    total_size_mb = sum(os.path.getsize(f) / (1024 * 1024) for f in chunk_files)
    print(f"Total output size: {total_size_mb:.1f} MB")
    print(f"Average chunk size: {total_size_mb / len(chunk_files):.1f} MB")
    
    # Check if any chunks exceed target size
    oversized = [f for f in chunk_files if os.path.getsize(f) / (1024 * 1024) > target_size_mb]
    if oversized:
        print(f"!  Warning: {len(oversized)} chunks exceed {target_size_mb} MB")
        for f in oversized:
            size_mb = os.path.getsize(f) / (1024 * 1024)
            print(f"  {Path(f).name}: {size_mb:.1f} MB")
    else:
        print("✓ All chunks are within size limit")
    
    return chunk_files


# # Sarcoid

# In[ ]:


chunk_files = split_gwas_file_pandas(
    'afr_gwas_results.tsv.gz', 
    output_dir="gwas_chunks"
)


# In[ ]:


chunk_files = split_gwas_file(
    'eur_gwas_results.tsv.gz', 
    output_dir="gwas_chunks"
)


# In[ ]:


chunk_files = split_gwas_file('condition__hpv_afr_amr_eur1.tsv')


# # Check Files

# In[ ]:


afr_df = pl.read_csv('afr_gwas_results.tsv.gz', separator='\t')


# In[ ]:


afr_df.height


# In[ ]:


afr_df.slice(1000,5)


# # Split Files (HPV)

# In[ ]:


#     'metal', '{bucket or my_bucket}/saige_gwas/v1_redo/metal/condition__hpv/condition__hpv_afr_amr_eur1.tsv',


# In[ ]:


# files_to_download = {
#     'eur': '{bucket or my_bucket}/saige_gwas/v1_redo/eur/condition__hpv/gwas_results.tsv.gz',
#     'afr': '{bucket or my_bucket}/saige_gwas/v1_redo/afr/condition__hpv/gwas_results.tsv.gz',
#     'amr': '{bucket or my_bucket}/saige_gwas/v1_redo/amr/condition__hpv/gwas_results.tsv.gz',
# }


# In[ ]:


# import subprocess
# for name, file in files_to_download.items():
#     # copy csv file to the bucket
#     args = ["gsutil", "cp", file, f"./hpv_download/{name}_gwas_results.tsv.gz"]
#     output = subprocess.run(args, capture_output=True)

#     # print output from gsutil
#     output.stderr


# In[ ]:


chunk_files = split_gwas_file(
    'hpv_download/afr_gwas_results.tsv.gz', 
    output_dir="hpv_download"
)


# In[ ]:


chunk_files = split_gwas_file(
    'hpv_download/eur_gwas_results.tsv.gz', 
    output_dir="hpv_download"
)


# In[ ]:


chunk_files = split_gwas_file(
    'hpv_download/amr_gwas_results.tsv.gz', 
    output_dir="hpv_download"
)


# In[ ]:


chunk_files = split_gwas_file('hpv_download/condition__hpv_afr_amr_eur1.tsv', output_dir='hpv_download')

