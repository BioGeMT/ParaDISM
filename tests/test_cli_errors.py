import io
import sys
import unittest
from contextlib import redirect_stderr
from unittest.mock import patch

import paradism


class CliErrorTest(unittest.TestCase):
    def test_memory_error_reports_incomplete_run(self):
        stderr = io.StringIO()

        with (
            patch.object(
                sys,
                "argv",
                ["paradism.py", "--read1", "reads.fq", "--reference", "ref.fa"],
            ),
            patch("paradism.run_with_arguments", side_effect=MemoryError),
            redirect_stderr(stderr),
        ):
            with self.assertRaisesRegex(SystemExit, "1"):
                paradism.main()

        message = stderr.getvalue()
        self.assertIn("exhausted the available memory", message)
        self.assertIn("no completion marker was written", message)


if __name__ == "__main__":
    unittest.main()
