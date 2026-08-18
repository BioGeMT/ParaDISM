"""Coordinate mapping between the ParaDISM PKD1 panel and GRCh38 chr16."""

from dataclasses import dataclass


@dataclass(frozen=True)
class PanelCoordinate:
    genomic_start0: int
    genomic_end0: int
    strand: str

    @property
    def length(self) -> int:
        return self.genomic_end0 - self.genomic_start0

    def genomic_to_contig_position(self, genomic_pos1: int) -> int | None:
        genomic_pos0 = genomic_pos1 - 1
        if not self.genomic_start0 <= genomic_pos0 < self.genomic_end0:
            return None
        if self.strand == "+":
            return genomic_pos0 - self.genomic_start0 + 1
        return self.genomic_end0 - genomic_pos0

    def contig_to_genomic_position(self, contig_pos1: int) -> int | None:
        if not 1 <= contig_pos1 <= self.length:
            return None
        if self.strand == "+":
            return self.genomic_start0 + contig_pos1
        return self.genomic_end0 - contig_pos1 + 1

    def genomic_to_contig_interval(
        self, genomic_start0: int, genomic_end0: int
    ) -> tuple[int, int] | None:
        overlap_start0 = max(genomic_start0, self.genomic_start0)
        overlap_end0 = min(genomic_end0, self.genomic_end0)
        if overlap_end0 <= overlap_start0:
            return None
        if self.strand == "+":
            return (
                overlap_start0 - self.genomic_start0,
                overlap_end0 - self.genomic_start0,
            )
        return (
            self.genomic_end0 - overlap_end0,
            self.genomic_end0 - overlap_start0,
        )

    def orient_allele(self, allele: str) -> str:
        if self.strand == "+":
            return allele
        return allele.translate(str.maketrans("ACGT", "TGCA"))[::-1]


# Exact end-to-end matches of benchmark/references/pkd1_panel.fa to GRCh38
# chromosome 16. Coordinates are 0-based and end-exclusive.
PANEL_COORDINATES: dict[str, PanelCoordinate] = {
    "PKD1": PanelCoordinate(2088107, 2136498, "-"),
    "PKD1P1": PanelCoordinate(16309740, 16334790, "+"),
    "PKD1P2": PanelCoordinate(16355623, 16378107, "+"),
    "PKD1P3": PanelCoordinate(14910950, 14936308, "+"),
    "PKD1P4": PanelCoordinate(18333799, 18353076, "-"),
    "PKD1P5": PanelCoordinate(18373920, 18402540, "-"),
    "PKD1P6": PanelCoordinate(15124641, 15155164, "-"),
}

# Genomic regions used for the manuscript evaluation. These are narrower than
# some panel contigs because the reference sequences include flanking bases.
EVALUATION_REGIONS: dict[str, tuple[int, int]] = {
    "PKD1": (2088706, 2135898),
    "PKD1P1": (16310132, 16334190),
    "PKD1P2": (16356222, 16377507),
    "PKD1P3": (14911549, 14935708),
    "PKD1P4": (18334398, 18352476),
    "PKD1P5": (18374519, 18402014),
    "PKD1P6": (15125137, 15154873),
}


def find_contig_position(genomic_pos1: int) -> tuple[str, int] | None:
    for contig, coordinate in PANEL_COORDINATES.items():
        contig_pos1 = coordinate.genomic_to_contig_position(genomic_pos1)
        if contig_pos1 is not None:
            return contig, contig_pos1
    return None
