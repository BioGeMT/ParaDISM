import os
import subprocess
import tempfile
import textwrap
import unittest
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DOWNLOAD_SCRIPT = PROJECT_ROOT / "benchmark/giab/download_giab_hg002_reads.sh"


class GiabDownloadTest(unittest.TestCase):
    def _run_script_with_wget(self, wget_body: str):
        temporary_directory = tempfile.TemporaryDirectory()
        self.addCleanup(temporary_directory.cleanup)
        test_root = Path(temporary_directory.name)
        script_directory = test_root / "benchmark/giab"
        script_directory.mkdir(parents=True)
        script_path = script_directory / DOWNLOAD_SCRIPT.name
        script_path.write_text(DOWNLOAD_SCRIPT.read_text())

        fake_bin = test_root / "bin"
        fake_bin.mkdir()
        wget_path = fake_bin / "wget"
        wget_path.write_text(textwrap.dedent(wget_body).lstrip())
        wget_path.chmod(0o755)

        environment = os.environ.copy()
        environment["PATH"] = f"{fake_bin}:{environment['PATH']}"
        result = subprocess.run(
            ["bash", str(script_path)],
            capture_output=True,
            text=True,
            env=environment,
        )
        return result, script_directory / "giab_hg002_reads"

    def test_downloads_and_merges_only_completed_files(self):
        result, reads_directory = self._run_script_with_wget(
            """
            #!/usr/bin/env bash
            set -euo pipefail
            while [[ $# -gt 0 ]]; do
                if [[ "$1" == "-O" ]]; then
                    destination="$2"
                    shift 2
                else
                    shift
                fi
            done
            printf x >> "$destination"
            """
        )

        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(len(list(reads_directory.glob("D1_*.fastq.gz"))), 68)
        self.assertEqual((reads_directory / "HG002_R1.fq.gz").stat().st_size, 34)
        self.assertEqual((reads_directory / "HG002_R2.fq.gz").stat().st_size, 34)
        self.assertEqual(list(reads_directory.glob("*.partial")), [])

    def test_failed_download_is_not_promoted_to_final_name(self):
        result, reads_directory = self._run_script_with_wget(
            """
            #!/usr/bin/env bash
            set -euo pipefail
            while [[ $# -gt 0 ]]; do
                if [[ "$1" == "-O" ]]; then
                    destination="$2"
                    shift 2
                else
                    shift
                fi
            done
            printf partial > "$destination"
            exit 1
            """
        )

        first_shard = reads_directory / "D1_S1_L001_R1_001.fastq.gz"
        self.assertNotEqual(result.returncode, 0)
        self.assertFalse(first_shard.exists())
        self.assertTrue(Path(f"{first_shard}.partial").exists())


if __name__ == "__main__":
    unittest.main()
