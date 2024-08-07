import argparse
import random
import gzip

def count_reads(file):
    with gzip.open(file, 'rt') if file.endswith('.gz') else open(file, 'r') as f:
        return sum(1 for line in f) // 4

def mix_reads(fastq_files, props):
    all_reads = []
    total_reads = [count_reads(file) for file in fastq_files]

    for i, fastq_file in enumerate(fastq_files):
        prop = float(props[i])
        num_reads = int(prop * total_reads[i])

        with gzip.open(fastq_file, 'rt') if fastq_file.endswith('.gz') else open(fastq_file, 'r') as fastq:
            reads = [next(fastq) + next(fastq) + next(fastq) + next(fastq) for _ in range(total_reads[i])]
            all_reads.extend(random.sample(reads, num_reads))

    random.shuffle(all_reads)
    for read in all_reads:
        print(read, end='')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Mix reads from multiple fastq files')
    parser.add_argument('fastq_files', metavar='FASTQ_FILE', nargs='+', help='Fastq files to be mixed')
    parser.add_argument('--props', metavar='PROP', nargs='+', help='Proportions of reads to be taken from each file')
    args = parser.parse_args()
    mix_reads(args.fastq_files, args.props)

