import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


class GiabCompletionTest(unittest.TestCase):
    def test_partial_final_outputs_are_not_treated_as_complete(self):
        project_root = Path(__file__).resolve().parents[1]
        script = project_root / "benchmark" / "giab" / "run_giab.sh"

        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            reads_dir = temp_path / "reads"
            reads_dir.mkdir()
            (reads_dir / "HG002_R1.fq").write_text("", encoding="utf-8")
            (reads_dir / "HG002_R2.fq").write_text("", encoding="utf-8")
            output_dir = temp_path / "partial-run"
            (output_dir / "final_outputs").mkdir(parents=True)

            result = subprocess.run(
                [
                    "bash",
                    str(script),
                    "--reads-dir",
                    str(reads_dir),
                    "--output-dir",
                    str(output_dir),
                ],
                cwd=project_root,
                text=True,
                capture_output=True,
            )

        self.assertNotEqual(0, result.returncode)
        self.assertIn("incomplete", result.stdout + result.stderr)

    def test_post_processing_fails_when_final_bams_are_missing(self):
        project_root = Path(__file__).resolve().parents[1]
        script = (
            project_root
            / "benchmark"
            / "giab"
            / "call_variants_raw_g60.sh"
        )
        reference = project_root / "benchmark" / "references" / "pkd1_panel.fa"

        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            input_dir = temp_path / "partial-run"
            input_dir.mkdir()
            output_dir = temp_path / "variant-calling"
            environment = os.environ.copy()
            environment["PATH"] = (
                f"{Path(sys.executable).parent}:{environment['PATH']}"
            )

            result = subprocess.run(
                [
                    "bash",
                    str(script),
                    "--input-dir",
                    str(input_dir),
                    "--output-dir",
                    str(output_dir),
                    "--reference",
                    str(reference),
                ],
                cwd=project_root,
                env=environment,
                text=True,
                capture_output=True,
            )

        self.assertNotEqual(0, result.returncode)
        self.assertIn(
            "ParaDISM BAM directory not found",
            result.stdout + result.stderr,
        )

    def test_post_processing_wrapper_rejects_missing_completion_marker(self):
        project_root = Path(__file__).resolve().parents[1]
        script = (
            project_root
            / "benchmark"
            / "giab"
            / "run_variant_calling_to_vcf_out.sh"
        )
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            input_dir = temp_path / "partial-run"
            input_dir.mkdir()
            benchmark_bed = temp_path / "benchmark.bed"
            benchmark_bed.write_text("PKD1\t0\t1\n", encoding="utf-8")

            result = subprocess.run(
                [
                    "bash",
                    str(script),
                    "--run-dir",
                    str(input_dir),
                    "--out-dir",
                    str(temp_path / "variant-calling"),
                    "--benchmark-bed",
                    str(benchmark_bed),
                ],
                cwd=project_root,
                text=True,
                capture_output=True,
            )

        self.assertNotEqual(0, result.returncode)
        self.assertIn("completion marker not found", result.stdout + result.stderr)


if __name__ == "__main__":
    unittest.main()
