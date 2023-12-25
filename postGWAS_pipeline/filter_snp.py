import pandas as pd
import os
import concurrent.futures

# I/O parameters
in_dir = os.path.expanduser('~/saige_pipe_test/')
out_dir = os.path.expanduser('~/bugbank_data/postgwas_regenie/sig_snp')
stem = 'your_file_stem'

# Filtering thresholds
min_mac = 300
min_info = 0.3
min_log10p = 7.301

def apply_filter(input_file, output_file):
    # Read the file
    df = pd.read_csv(input_file, sep=' ', compression='infer')

    # Calculate Minor Allele Frequency (MAF)
    df['MAF'] = df['A1FREQ'].apply(lambda x: min(x, 1-x))

    # Calculate MAC
    df['MAC'] = 2 * df['N'] * df['MAF']

    # Apply filters
    filtered_df = df[(df['MAC'] >= min_mac) & 
                     (df['INFO'] >= min_info) & 
                     (df['PVAL'] >= 10 ** (-min_log10p))]

    # Write to file
    filtered_df.to_csv(output_file, sep=' ', index=False)
    print(f'{output_file} done')

# Function to get output file path
def get_output_file(input_file):
    return input_file.replace('summary', 'filtered').replace('.txt.gz', f'_mac{min_mac}info{min_info}p{min_log10p}.txt')

# Process each file using ThreadPoolExecutor
with concurrent.futures.ThreadPoolExecutor(max_workers=20) as executor:
    # Create a list of futures
    futures = [executor.submit(apply_filter, summary_file, get_output_file(summary_file)) for summary_file in summary_files]

    # Wait for all futures to complete
    for future in concurrent.futures.as_completed(futures):
        print(future.result())

