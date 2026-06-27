import tempfile
import unittest
from pathlib import Path

from src.pipeline.paradism_algo import load_msa, process_sam_to_dict, _process_sam_to_dict_parallel


class ParallelAssignmentTest(unittest.TestCase):
    def test_parallel_assignment_matches_serial_assignment(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            msa_path = temp_path / "ref.aln"
            sam_path = temp_path / "reads.sam"

            msa_path.write_text(
                ">GENE_A\n"
                "ACGA\n"
                ">GENE_B\n"
                "ACGT\n",
                encoding="utf-8",
            )
            sam_path.write_text(
                "@HD\tVN:1.6\tSO:unsorted\n"
                "@SQ\tSN:GENE_A\tLN:4\n"
                "@SQ\tSN:GENE_B\tLN:4\n"
                "read_a\t0\tGENE_A\t1\t255\t4M\t*\t0\t0\tACGA\tIIII\n"
                "read_b\t0\tGENE_B\t1\t255\t4M\t*\t0\t0\tACGT\tIIII\n"
                "read_none\t0\tGENE_A\t1\t255\t4M\t*\t0\t0\tACGC\tIIII\n"
                "read_multi\t0\tGENE_A\t1\t255\t4M\t*\t0\t0\tACGA\tIIII\n"
                "read_multi\t0\tGENE_B\t1\t255\t4M\t*\t0\t0\tACGT\tIIII\n",
                encoding="utf-8",
            )

            msa, seq_to_aln, gene_names = load_msa(str(msa_path))
            serial = process_sam_to_dict(str(sam_path), msa, seq_to_aln, gene_names)
            parallel = process_sam_to_dict(str(sam_path), msa, seq_to_aln, gene_names, workers=3)

        self.assertEqual(serial, parallel)
        self.assertEqual(
            {
                "read_a": "GENE_A",
                "read_b": "GENE_B",
                "read_none": "NONE",
                "read_multi": "NONE",
            },
            parallel,
        )

    def test_parallel_assignment_merges_evidence_across_chunks(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            msa_path = temp_path / "ref.aln"
            sam_path = temp_path / "reads.sam"

            msa_path.write_text(
                ">GENE_A\n"
                "ACGA\n"
                ">GENE_B\n"
                "ACGT\n",
                encoding="utf-8",
            )
            sam_path.write_text(
                "@HD\tVN:1.6\tSO:unsorted\n"
                "@SQ\tSN:GENE_A\tLN:4\n"
                "@SQ\tSN:GENE_B\tLN:4\n"
                "split_read\t0\tGENE_A\t1\t255\t4M\t*\t0\t0\tACGA\tIIII\n"
                "split_read\t0\tGENE_B\t1\t255\t4M\t*\t0\t0\tACGT\tIIII\n",
                encoding="utf-8",
            )

            msa, seq_to_aln, gene_names = load_msa(str(msa_path))
            parallel = _process_sam_to_dict_parallel(
                str(sam_path),
                msa,
                seq_to_aln,
                gene_names,
                min_anchors=1,
                workers=2,
                chunk_size=1,
            )

        self.assertEqual({"split_read": "NONE"}, parallel)


if __name__ == "__main__":
    unittest.main()
