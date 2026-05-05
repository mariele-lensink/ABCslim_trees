#!/usr/bin/env python3
import csv
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
GENOMEINFO_DIR = REPO_ROOT / "genomeinfo"
FIGURES_DIR = REPO_ROOT / "figures"
OUT_PATH = FIGURES_DIR / "genome_segments.svg"


SEGMENTS = [
    ("Chr1", GENOMEINFO_DIR / "gene100_chr1_10000000.csv", GENOMEINFO_DIR / "intergene100_chr1_10000000.csv"),
    ("Chr2", GENOMEINFO_DIR / "gene100_chr2_18000000.csv", GENOMEINFO_DIR / "intergene100_chr2_18000000.csv"),
    ("Chr3", GENOMEINFO_DIR / "gene100_chr3_22000000.csv", GENOMEINFO_DIR / "intergene100_chr3_22000000.csv"),
    ("Chr4", GENOMEINFO_DIR / "gene100_chr4_2000000.csv", GENOMEINFO_DIR / "intergene100_chr4_2000000.csv"),
    ("Chr5", GENOMEINFO_DIR / "gene100_chr5_1.csv", GENOMEINFO_DIR / "intergene100_chr5_1.csv"),
]


def read_intervals(path):
    intervals = []
    with open(path, newline="") as handle:
        for row in csv.reader(handle):
            if row:
                intervals.append((int(row[0]), int(row[1])))
    return intervals


def segment_length(gene_intervals, intergene_intervals):
    return max(stop for _, stop in gene_intervals + intergene_intervals)


def rect(x, y, width, height, fill, stroke, stroke_width):
    return (
        f'<rect x="{x:.2f}" y="{y:.2f}" width="{width:.2f}" height="{height:.2f}" '
        f'fill="{fill}" stroke="{stroke}" stroke-width="{stroke_width}"/>'
    )


def text(x, y, value, size=14, anchor="start", weight="normal"):
    safe = (
        str(value)
        .replace("&", "&amp;")
        .replace("<", "&lt;")
        .replace(">", "&gt;")
    )
    return (
        f'<text x="{x:.2f}" y="{y:.2f}" font-size="{size}" text-anchor="{anchor}" '
        f'font-family="Helvetica, Arial, sans-serif" font-weight="{weight}" fill="#111827">{safe}</text>'
    )


def main():
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)

    loaded = []
    max_length = 0
    for label, gene_path, intergene_path in SEGMENTS:
        gene_intervals = read_intervals(gene_path)
        intergene_intervals = read_intervals(intergene_path)
        length = segment_length(gene_intervals, intergene_intervals)
        max_length = max(max_length, length)
        loaded.append((label, gene_intervals, intergene_intervals, length))

    width = 1000
    left_margin = 180
    right_margin = 40
    top_margin = 60
    row_gap = 34
    bar_height = 34
    row_height = bar_height + row_gap
    plot_width = width - left_margin - right_margin
    height = top_margin + len(loaded) * row_height + 50
    blue = "#3B82F6"

    svg = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        rect(0, 0, width, height, "#FFFFFF", "#FFFFFF", 0),
        text(left_margin, 28, "Genic and Intergenic Layout for the Five Simulation Segments", size=20, weight="bold"),
        text(left_margin, 48, "Blue = genic intervals, white = intergenic intervals; all rows share the same length scale.", size=11),
    ]

    for idx, (label, gene_intervals, _, length) in enumerate(loaded):
        y = top_margin + idx * row_height
        svg.append(text(left_margin - 18, y + 14, label, size=14, anchor="end", weight="bold"))
        svg.append(text(left_margin - 18, y + 29, f"{length:,} bp", size=11, anchor="end"))
        svg.append(rect(left_margin, y, plot_width * (length / max_length), bar_height, "#FFFFFF", "#1F2937", 1.2))

        for start, stop in gene_intervals:
            x = left_margin + ((start - 1) / max_length) * plot_width
            gene_width = ((stop - start + 1) / max_length) * plot_width
            svg.append(rect(x, y, gene_width, bar_height, blue, blue, 0))

    axis_y = top_margin + len(loaded) * row_height - row_gap / 2
    svg.append(f'<line x1="{left_margin:.2f}" y1="{axis_y:.2f}" x2="{left_margin + plot_width:.2f}" y2="{axis_y:.2f}" stroke="#1F2937" stroke-width="1"/>')
    svg.append(text(left_margin, axis_y + 20, "0", size=11))
    svg.append(text(left_margin + plot_width, axis_y + 20, f"{max_length:,} bp", size=11, anchor="end"))
    svg.append("</svg>")

    OUT_PATH.write_text("\n".join(svg))
    print(OUT_PATH)


if __name__ == "__main__":
    main()
