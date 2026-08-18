import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch


PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from pipeline.executor import SimpleParaDISMExecutor
from paradism import build_parser


class IntermediateSamCompressionTest(unittest.TestCase):
    def test_cli_flag_enables_intermediate_compression(self):
        args = build_parser().parse_args(
            [
                "--read1",
                "reads.fq",
                "--reference",
                "reference.fa",
                "--compress-intermediate-sam",
            ]
        )

        self.assertTrue(args.compress_intermediate_sam)

    def test_successful_conversion_replaces_sam_atomically(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            sam_path = temp_path / "mapped_reads.sam"
            sam_path.write_text("sam", encoding="utf-8")
            executor = SimpleParaDISMExecutor(output_dir=temp_path / "output")

            def fake_run(command, **_kwargs):
                if command[1] == "view":
                    output_path = Path(command[command.index("-o") + 1])
                    output_path.write_bytes(b"bam")
                return subprocess.CompletedProcess(command, 0)

            with patch("pipeline.executor.subprocess.run", side_effect=fake_run):
                bam_path = executor._compress_intermediate_sam(sam_path, threads=3)

            self.assertEqual(temp_path / "mapped_reads.bam", bam_path)
            self.assertTrue(bam_path.exists())
            self.assertFalse(sam_path.exists())
            self.assertFalse(Path(f"{bam_path}.partial").exists())

    def test_failed_conversion_preserves_sam(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            sam_path = temp_path / "mapped_reads.sam"
            sam_path.write_text("sam", encoding="utf-8")
            executor = SimpleParaDISMExecutor(output_dir=temp_path / "output")

            with patch(
                "pipeline.executor.subprocess.run",
                side_effect=subprocess.CalledProcessError(1, ["samtools"]),
            ):
                with self.assertRaises(subprocess.CalledProcessError):
                    executor._compress_intermediate_sam(sam_path, threads=3)

            self.assertTrue(sam_path.exists())
            self.assertFalse((temp_path / "mapped_reads.bam").exists())
            self.assertFalse((temp_path / "mapped_reads.bam.partial").exists())


if __name__ == "__main__":
    unittest.main()
