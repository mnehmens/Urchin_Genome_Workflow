import sys
import re

def read_genome_assembly(assembly_file):
    genome = {}
    with open(assembly_file, 'r') as af:
        chrom = ""
        seq = []
        for line in af:
            if line.startswith(">"):
                if chrom:
                    genome[chrom] = ''.join(seq)
                chrom = line[1:].strip().split()[0]
                seq = []
            else:
                seq.append(line.strip())
        if chrom:
            genome[chrom] = ''.join(seq)
    return genome

def write_genome_assembly(genome, output_file):
    with open(output_file, 'w') as of:
        for chrom, seq in genome.items():
            of.write(f">{chrom}\n")
            for i in range(0, len(seq), 60):
                of.write(seq[i:i+60] + "\n")

def trim_long_stretches_of_ns(genome, max_n_length=5000, replace_length=10000):
    for chrom in genome:
        seq = genome[chrom]
        # Regular expression to find stretches of N longer than max_n_length
        pattern = f'N{{{max_n_length},}}'
        genome[chrom] = re.sub(pattern, 'N' * replace_length, seq)
    return genome

def main(assembly_file, output_file, max_n_length):
    max_n_length = int(max_n_length)
    genome = read_genome_assembly(assembly_file)
    genome = trim_long_stretches_of_ns(genome, max_n_length)
    write_genome_assembly(genome, output_file)

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python script.py <genome.fasta> <output.fasta> <max_n_length>")
        sys.exit(1)
    assembly_file = sys.argv[1]
    output_file = sys.argv[2]
    max_n_length = sys.argv[3]
    main(assembly_file, output_file, max_n_length)

