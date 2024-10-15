#!/bin/bash

# Check if the number of sources is provided
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <number_of_sources>"
    exit 1
fi

num_sources=$1
downsampling_results="downsampling_results"

# Create the downsampling_results directory if it doesn't exist
mkdir -p $downsampling_results

# Function to run experiments with downsampling
run_experiment() {
    source_number=$1
    input_file=$2
    output_prefix=$3
    read_count=$4
    record_count=$((read_count * 4))  # Multiply by 4 for each read line in FASTQ
    temp_file="samples/${output_prefix}_${read_count}.fastq"  # Temporary file for each iteration
    echo "Running downsampling for $source_number sources with $read_count reads..."

    # Shuffle, take the top 'n' reads, and write to a temporary file

    zcat $input_file | \
    ./shuf.sh | \
    head -n $record_count > $temp_file

    # --contamination-mode --source-assignments 0,0
    # Run vgan using the temporary file as input
    nice -19 ./../bin/vgan trailmix -t 85 -fq1 $temp_file --iter 100000 --burnin 1 -k $source_number --chains 5 \
                           -o $downsampling_results/${output_prefix}_$read_count \
                           --deam3p ../share/damageProfiles/dmid3p.prof --deam5p ../share/damageProfiles/dmid5p.prof \
                           -z tempdir --depth 10 \

    # Cleanup: Remove the temporary file
    #rm $temp_file
}

# Run the experiments based on the number of sources
case $num_sources in
    1)
        for read_count in {700..700..-100}
        do
            run_experiment 1 "../../TrailMix_Experiment_Data/simulations/1/endo/Node2605_gen_0_dhigh_bact.fq.gz" "Node2605_dhigh" $read_count
        done
        ;;
    2)
        for read_count in {30000..30000..-100}
        do
            run_experiment 2 "/home/projects2/hominin/TrailMix_Experiment_Data/simulations/2/0.2_0.8/C4e_V2a1a_gen_0_dmid_dmid.fq.gz" "CV_mid_mid" $read_count
        done
        ;;
    3)
        for read_count in {16000..16000..-100}
        do
            run_experiment 3 "../../TrailMix_Experiment_Data/simulations/3/0.3_0.3_0.4/L2a1d_S3_P9a_gen_0_dmid_dmid.fq.gz" "LSP_dmid_16000" $read_count
        done
        ;;
    *)
        echo "Invalid number of sources. Please enter a number between 1 and 3."
        exit 1
        ;;
esac

