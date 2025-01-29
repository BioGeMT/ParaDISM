import argparse
from Bio import SeqIO
from collections import OrderedDict

def map_reference_to_msa(reference_fasta_path, msa_fasta_path, output_tsv_path):
  
    # Load reference and MSA sequences into ordered dictionaries to maintain consistent column order
    reference_sequences = OrderedDict((record.id, str(record.seq)) 
                                      for record in SeqIO.parse(reference_fasta_path, "fasta"))
    msa_sequences = OrderedDict((record.id, str(record.seq)) 
                                 for record in SeqIO.parse(msa_fasta_path, "fasta"))

    with open(output_tsv_path, "w") as output_file:
        # Write header for output TSV
        sequence_ids = list(msa_sequences.keys())
        position_columns = [f"{seq_id}_Position" for seq_id in sequence_ids]
        base_columns = [f"{seq_id}_Base" for seq_id in sequence_ids]
        # Fixed string formatting to avoid f-string with backslashes
        header = "MSA_Position\t" + "\t".join(position_columns) + "\t" + "\t".join(base_columns) + "\n"
        output_file.write(header)

        # Process each reference sequence
        for reference_id, reference_seq in reference_sequences.items():
            msa_seq = msa_sequences[reference_id]
            reference_positions = {seq_id: -1 for seq_id in msa_sequences}

            # Iterate through each position in the MSA sequence
            for msa_position, reference_base in enumerate(msa_seq):
                row = [msa_position]
                positions = []
                bases = []

                for seq_id, seq in msa_sequences.items():
                    base = seq[msa_position]
                    if base != "-":
                        reference_positions[seq_id] += 1
                    positions.append(reference_positions[seq_id] if base != "-" else "-")
                    bases.append(base)

                row.extend(positions)
                row.extend(bases)

                output_line = "\t".join(map(str, row)) + "\n"
                output_file.write(output_line)

def main():
    parser = argparse.ArgumentParser(
        description="Map reference sequence positions to MSA positions and output valid bases along with individual sequence bases."
    )
    parser.add_argument(
        "--reference_fasta",
        required=True,
        help="Path to reference sequences in FASTA format"
    )
    parser.add_argument(
        "--msa_file",
        required=True,
        help="Path to MSA sequences in FASTA format"
    )
    parser.add_argument(
        "--output_file",
        required=True,
        help="Path for output TSV mapping file"
    )
    args = parser.parse_args()

    map_reference_to_msa(args.reference_fasta, args.msa_file, args.output_file)

if __name__ == "__main__":
    main()