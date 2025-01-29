import pysam
import argparse

def process_sam_to_alignment_tsv(sam_filepath, output_tsv_filepath):

    with pysam.AlignmentFile(sam_filepath, "r") as samfile, open(output_tsv_filepath, "w") as outfile:
        # Write header
        outfile.write("Read_Name\tRead_Position\tReference_Position\tRead_Base\tReference_Base\tType\tReference_Name\n")

        for read in samfile.fetch():
            read_name = read.query_name + ('+' if not read.is_reverse else '-')
            ref_name = samfile.get_reference_name(read.reference_id)

            query_sequence = read.query_sequence
            reference_start = read.reference_start

            for query_pos, ref_pos in read.get_aligned_pairs(matches_only=False):
                if query_pos is None:  # Deletion
                    read_base = "-"
                    ref_base = read.get_reference_sequence()[ref_pos - reference_start]
                    alignment_type = "deletion"
                elif ref_pos is None:  # Insertion
                    read_base = query_sequence[query_pos]
                    ref_base = "-"
                    alignment_type = "insertion"
                else:  # Match or mismatch
                    read_base = query_sequence[query_pos]
                    ref_base = read.get_reference_sequence()[ref_pos - reference_start]
                    alignment_type = "match" if read_base == ref_base else "mismatch"

                # Write data
                outfile.write(
                    f"{read_name}\t{query_pos if query_pos is not None else '-'}\t{ref_pos if ref_pos is not None else '-'}\t{read_base}\t{ref_base}\t{alignment_type}\t{ref_name}\n"
                )

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--sam", required=True, help="Path to input SAM file")
    parser.add_argument("--output", required=True, help="Path for output TSV file")
    args = parser.parse_args()

    process_sam_to_alignment_tsv(args.sam, args.output)