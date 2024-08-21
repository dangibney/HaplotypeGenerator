#!/bin/bash

# File to store recombination counts
output_file="recombination_counts.txt"

# Clear the file if it exists
> $output_file

# Temporary file to store individual recombination values per batch
temp_file="temp_recombination_counts.txt"
> $temp_file

# Total number of runs
total_runs=10

# Number of parallel jobs (adjust based on your CPU cores)
num_parallel_jobs=6

# Variables to keep track of cumulative sum and count across all batches
cumulative_sum=0
total_count=0

# Function to run the program and capture recombinations
run_program() {
    local run_id=$1
    local output=$(./haplotype_generator test_data/MHC-CHM13.0.gfa.gz generated_hap.fa)

    # Extract the number of recombinations
    recombinations=$(echo "$output" | grep "Number of recombinations:" | awk '{print $4}')

    # Append the number of recombinations to a temporary file and the main output file
    echo "$recombinations" >> $temp_file
    echo "$recombinations" >> $output_file

    # Print progress to terminal
    echo "Run $run_id complete with $recombinations recombinations."
}

# Run the C++ program in parallel
for ((i=1; i<=total_runs; i++)); do
    run_program $i &  # Run in background
    
    # Limit the number of parallel jobs
    if (( i % num_parallel_jobs == 0 )); then
        wait  # Wait for current batch of jobs to finish
        
        # Process the current batch and update the cumulative sum and count
        while read -r recombinations; do
            cumulative_sum=$((cumulative_sum + recombinations))
            total_count=$((total_count + 1))
            
            # Calculate and print the cumulative running average
            running_average=$(echo "$cumulative_sum / $total_count" | bc -l)
            echo "After $total_count runs, cumulative running average: $running_average"
        done < $temp_file
        
        # Clear temp file for the next batch
        > $temp_file
    fi
done

# Wait for any remaining background jobs to complete
wait

# Process any leftover runs from the final batch
while read -r recombinations; do
    cumulative_sum=$((cumulative_sum + recombinations))
    total_count=$((total_count + 1))
    
    # Calculate and print the cumulative running average
    running_average=$(echo "$cumulative_sum / $total_count" | bc -l)
    echo "After $total_count runs, cumulative running average: $running_average"
done < $temp_file

# Now call a Python script to plot the distribution
python3 src/plot_recombinations.py
