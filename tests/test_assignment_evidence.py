import unittest
from collections import defaultdict

from src.pipeline.paradism_algo import _assign_from_collected_evidence


class AssignmentEvidenceTest(unittest.TestCase):
    def test_assignment_does_not_create_missing_anchor_entries(self):
        gene_names = ["GENE_A", "GENE_B"]
        anchor_cols = {
            "GENE_A": defaultdict(
                set,
                {
                    "read_a": {1},
                    "read_multi": {1},
                },
            ),
            "GENE_B": defaultdict(
                set,
                {
                    "read_b": {2},
                    "read_multi": {2},
                },
            ),
        }
        c2_false = {
            "GENE_A": set(),
            "GENE_B": {"read_b"},
        }
        all_qnames = {"read_a", "read_b", "read_multi", "read_none"}
        keys_before = {
            gene: set(anchor_cols[gene])
            for gene in gene_names
        }

        assignments = _assign_from_collected_evidence(
            anchor_cols,
            c2_false,
            all_qnames,
            gene_names,
            min_anchors=1,
        )

        self.assertEqual(
            {
                "read_a": "GENE_A",
                "read_b": "NONE",
                "read_multi": "NONE",
                "read_none": "NONE",
            },
            assignments,
        )
        self.assertEqual(
            keys_before,
            {gene: set(anchor_cols[gene]) for gene in gene_names},
        )


if __name__ == "__main__":
    unittest.main()
