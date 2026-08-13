import sys
import unittest
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
GIAB_DIR = PROJECT_ROOT / "benchmark" / "giab"
sys.path.insert(0, str(GIAB_DIR))

from make_benchmark_gene_coords_bed import map_chr16_to_gene_coords
from pkd1_panel_coordinates import EVALUATION_REGIONS, PANEL_COORDINATES


class GiabCoordinateTest(unittest.TestCase):
    def test_panel_mapping_matches_reference_lengths(self):
        reference = PROJECT_ROOT / "benchmark" / "references" / "pkd1_panel.fa"
        lengths = {}
        current = None
        with open(reference, encoding="utf-8") as handle:
            for line in handle:
                if line.startswith(">"):
                    current = line[1:].split()[0]
                    lengths[current] = 0
                else:
                    lengths[current] += len(line.strip())

        self.assertEqual(
            {name: coordinate.length for name, coordinate in PANEL_COORDINATES.items()},
            lengths,
        )

    def test_plus_and_minus_strand_endpoints(self):
        plus = PANEL_COORDINATES["PKD1P1"]
        self.assertEqual(
            1, plus.genomic_to_contig_position(plus.genomic_start0 + 1)
        )
        self.assertEqual(
            plus.length, plus.genomic_to_contig_position(plus.genomic_end0)
        )
        self.assertEqual(plus.genomic_start0 + 1, plus.contig_to_genomic_position(1))

        minus = PANEL_COORDINATES["PKD1"]
        self.assertEqual(
            minus.length,
            minus.genomic_to_contig_position(minus.genomic_start0 + 1),
        )
        self.assertEqual(1, minus.genomic_to_contig_position(minus.genomic_end0))
        self.assertEqual(minus.genomic_end0, minus.contig_to_genomic_position(1))
        self.assertEqual("AC", minus.orient_allele("GT"))

    def test_evaluation_regions_are_mapped_into_contig_coordinates(self):
        intervals = list(EVALUATION_REGIONS.values())
        mapped = map_chr16_to_gene_coords(intervals)

        self.assertEqual([(600, 47792)], mapped["PKD1"])
        self.assertEqual([(291, 30027)], mapped["PKD1P6"])


if __name__ == "__main__":
    unittest.main()
