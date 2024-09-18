#!/bin/bash
# Usage: ./run_schmutzi.sh <bam_dir> <library_type>
# Make sure this script has executable permissions: chmod +x run_schmutzi.sh

if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <bam_dir> <library_type>"
    exit 1
fi

bam_dir=$(realpath "$1")
library_type="$2"
reference_genome="/home/projects2/hominin/benchmarking/samples/human_MT.fa"
freqs_dir="/home/ctools/schmutzi-1.5.6/share/schmutzi/alleleFreqMT/197/freqs/"
output_dir=$(realpath "schmutzi_out_CV")
mkdir -p "$output_dir"

# Check if bam_dir is a directory
if [ ! -d "$bam_dir" ]; then
    echo "$bam_dir is not a directory."
    exit 1
fi

# Find all .bam files
bam_files=$(find "$bam_dir" -name "*.bam")

# Check if there are .bam files in the directory
if [ -z "$bam_files" ]; then
    echo "No .bam files found in the directory."
    exit 1
fi

contamination_file="$output_dir/contamination_estimates.tsv"
echo -e "Sample Name\tEndogenous/Total ratio" > "$contamination_file"

for bam_file in $bam_files; do
    base_name=$(basename "${bam_file%.*}")
    output_dir_sample="$output_dir/${base_name}_schmutzi_output"
    mkdir -p "$output_dir_sample"

    # Check if the BAM file is sorted and indexed
    if ! samtools view -H "$bam_file" | grep -q "SO:coordinate"; then
        echo "Sorting BAM file: $bam_file"
        samtools sort "$bam_file" -o "$output_dir_sample/$base_name.sorted.bam"
        bam_file="$output_dir_sample/$base_name.sorted.bam"
    fi

    if [ ! -f "${bam_file}.bai" ]; then
        echo "Indexing BAM file: $bam_file"
        samtools index "$bam_file"
    fi

    # Contamination estimation
    cont_deam_pl="/home/ctools/schmutzi-1.5.6/src/contDeam.pl"
    cont_deam_out="$output_dir_sample/contDeam"
    cont_deam_cmd="$cont_deam_pl -r $reference_genome --library $library_type --lengthDeam 5 --out $cont_deam_out $bam_file"
    eval "$cont_deam_cmd"

    # Parse contDeam output and write to file
    endogenous_ratio=$(head -1 "$output_dir_sample/contDeam.cont.est" | awk '{print $1}')
    echo -e "$base_name\t$endogenous_ratio" >> "$contamination_file"

    # Run schmutzi
    schmutzi_caller="/home/ctools/schmutzi-1.5.6/src/schmutzi.pl"
    schmutzi_caller_out="$output_dir_sample/schmutzi"
    schmutzi_caller_cmd="$schmutzi_caller --iterations 100 --t 45 --ref $reference_genome --uselength $cont_deam_out $freqs_dir $bam_file"
    echo $schmutzi_caller_cmd
    eval "$schmutzi_caller_cmd"
done

# Display all results
cat "$contamination_file"
