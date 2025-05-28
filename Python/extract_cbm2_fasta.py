import argparse
from Bio import SeqIO

def parse_tsv(tsv_file):
    accession_to_matches = {}
    with open(tsv_file, 'r') as tsv:
        next(tsv)  # Skip the header
        for line in tsv:
            fields = line.strip().split('\t')
            accession = fields[0]
            matches = fields[-1]
            accession_to_matches[accession] = matches
    return accession_to_matches

def extract_sequences(fasta_file, accession_to_matches):
    sequences = {}
    with open(fasta_file, 'r') as fasta:
        for record in SeqIO.parse(fasta, 'fasta'):
            accession = record.id.split('|')[0]
            if accession in accession_to_matches:
                matches = accession_to_matches[accession]
                match_ranges = matches.split(',')
                for match_range in match_ranges:
                    start, end = map(int, match_range.split('..'))
                    sequence = record.seq[start - 1:end]
                    sequence_id = f"{record.id}_{start}-{end}"
                    sequences[sequence_id] = sequence
    return sequences

def write_fasta(output_file, sequences):
    with open(output_file, 'w') as fasta_out:
        for seq_id, sequence in sequences.items():
            fasta_out.write(f">{seq_id}\n{sequence}\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate a new FASTA file based on TSV and FASTA inputs.")
    parser.add_argument("tsv_file", help="Input TSV file")
    parser.add_argument("fasta_file", help="Input FASTA file")
    parser.add_argument("output_file", help="Output FASTA file")
    args = parser.parse_args()

    accession_to_matches = parse_tsv(args.tsv_file)
    sequences = extract_sequences(args.fasta_file, accession_to_matches)
    write_fasta(args.output_file, sequences)
    print(f"New FASTA file '{args.output_file}' generated successfully.")

