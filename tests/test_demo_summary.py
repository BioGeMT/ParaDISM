import unittest

from demo.summarize_demo import summarize_assignments


class DemoSummaryTest(unittest.TestCase):
    def test_counts_correct_unassigned_and_incorrect_assignments(self):
        truth_by_read = {
            "GENE_A_1": "GENE_A",
            "GENE_A_2": "GENE_A",
            "GENE_B_1": "GENE_B",
            "GENE_B_2": "GENE_B",
        }
        assignments = {
            "GENE_A_1": "GENE_A",
            "GENE_A_2": "NONE",
            "GENE_B_1": "GENE_A",
            "GENE_B_2": "GENE_B",
        }

        self.assertEqual(
            (
                ("GENE_A", 2, 1, 1, 0),
                ("GENE_B", 2, 1, 0, 1),
            ),
            summarize_assignments(truth_by_read, assignments),
        )


if __name__ == "__main__":
    unittest.main()
