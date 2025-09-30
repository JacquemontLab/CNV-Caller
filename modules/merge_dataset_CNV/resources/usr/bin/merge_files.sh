#!/bin/bash

###############################################################################
# Script Name: merge_files.sh
# Description: This script merges multiple CNV files (gzipped or uncompressed)
#              by splitting them into smaller chunks for parallel processing,
#              while preserving the header from the first file. The results
#              are stored in a tab-delimited, gzipped TSV file for downstream
#              analysis.
#
# Usage: ./merge_files.sh --directory <directory> --chunk_size <size_of_chunks> 
#                              --num_cores <number_of_cores> --output_file <output_filename> 
#                              --file_extension <file_extension>
#
# Author: Florian Bénitière
# Date: April 2025
###############################################################################

set -euo pipefail

# Parse command line arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --directory) directory="$2"; shift ;;
        --chunk_size) chunk_size="$2"; shift ;;
        --num_cores) num_cores="$2"; shift ;;
        --output_file) output_file="$2"; shift ;;
        --file_extension) file_extension="$2"; shift ;;
        *) echo "Unknown parameter: $1"; exit 1 ;;
    esac
    shift
done

# Validate input arguments
if [[ -z "$directory" || -z "$chunk_size" || -z "$num_cores" || -z "$output_file" || -z "$file_extension" ]]; then
    echo "Error: Missing required arguments."
    echo "Usage: ./merge_files.sh --directory <directory> --chunk_size <size_of_chunks> --num_cores <number_of_cores> --output_file <output_filename> --file_extension <file_extension>"
    exit 1
fi

if [ "$directory" == "." ]; then
    directory=""
fi

# List all the files with the specified extension in the directory and split them into chunks
files_list=$(ls ${directory}*${file_extension})
files_array=($files_list)


# Function to process a chunk of files
process_chunk() {
    chunk_index=$1
    start_index=$((chunk_index * chunk_size))
    end_index=$((start_index + chunk_size - 1))

    # Create a chunk of files
    chunk_files=("${files_array[@]:$start_index:$chunk_size}")
    cat "${chunk_files[@]}" > "merged_chunk_${chunk_index}"
}

# Define how many files per chunk
num_chunks=$(( (${#files_array[@]} - 1) / $chunk_size ))

# Run the merging process in parallel using background jobs with a limit of available cores
current_jobs=0
for chunk_index in $(seq 0 $((num_chunks))); do
    # Start a background job
    process_chunk $chunk_index &
    
    # Increment the number of current jobs
    current_jobs=$((current_jobs + 1))
    
    # Wait if we've reached the maximum number of concurrent jobs (based on available cores)
    if [ $current_jobs -ge $num_cores ]; then
        wait -n  # Wait for any background job to finish
        current_jobs=$((current_jobs - 1))  # Decrease the job count
    fi
done

# Wait for all remaining background jobs to finish
wait

# Merge all chunk files into the final output
cat merged_chunk_* > global_merged

# Clean up the chunk files
rm merged_chunk_*

# Capture the header from the first file and remove it from the global merged file
header=$(head -n 1 "${files_array[0]}")  # Capture the header from the first file
{ 
  echo "$header"  # Output the header once
  grep -v "$header" "global_merged"  # Skip the header in subsequent files and concatenate
} > $output_file

# Clean up the intermediate merged file
rm global_merged

echo "Merging complete. Final output saved as $output_file."
