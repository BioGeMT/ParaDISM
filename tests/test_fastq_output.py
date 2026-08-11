import gzip
import tempfile
import unittest
from pathlib import Path

from Bio import SeqIO

from src.pipeline.paradism_algo import write_fastq_outputs


class FastqOutputTest(unittest.TestCase):
    def test_writes_assigned_reads_from_gzipped_paired_fastqs(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            r1_path = temp_path / "reads_R1.fq.gz"
            r2_path = temp_path / "reads_R2.fq.gz"
            output_dir = temp_path / "outputs"

            with gzip.open(r1_path, "wt", encoding="utf-8") as r1_handle:
                r1_handle.write("@read_a/1\nACGT\n+\nIIII\n")
            with gzip.open(r2_path, "wt", encoding="utf-8") as r2_handle:
                r2_handle.write("@read_a/2\nTGCA\n+\nIIII\n")

            genes = write_fastq_outputs(
                {"read_a": "GENE_A"},
                str(r1_path),
                str(r2_path),
                str(output_dir),
            )
            records = [
                (record.id, str(record.seq))
                for record in SeqIO.parse(output_dir / "GENE_A.fq", "fastq")
            ]

            self.assertEqual(
                (["GENE_A"], [("read_a/1", "ACGT"), ("read_a/2", "TGCA")]),
                (genes, records),
            )


if __name__ == "__main__":
    unittest.main()
