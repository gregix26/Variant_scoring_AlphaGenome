"""
VCF Merger Script
-----------------
Traverses: main_dir -> sub1 -> sub2(s) -> *.vcf files

Usage:
    python merge_vcfs.py --root /path/to/main_directory
"""

import argparse
from pathlib import Path

# ── Standard VCF header ───────────────────────────────────────────────────────
VCF_HEADER = """##fileformat=VCFv4.2
##reference=GRCh38
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"""

# ── CLI ───────────────────────────────────────────────────────────────────────
parser = argparse.ArgumentParser(description='Merge VCF files per sub1, output named after sub1')
parser.add_argument('--root', required=True, help='Path to main directory')
args = parser.parse_args()

root = Path(args.root).resolve()
if not root.exists():
    raise FileNotFoundError(f"Root directory not found: {root}")

print(f"Root directory: {root}\n")

total_merged   = 0
total_files    = 0
total_variants = 0

# ── Walk: root -> sub1 -> sub2 -> *.vcf ──────────────────────────────────────
for sub1 in sorted(root.iterdir()):
    if not sub1.is_dir():
        continue

    # Collect all .vcf files from every sub2 inside this sub1
    vcf_files = sorted(sub1.glob("*/*.vcf"))

    if not vcf_files:
        print(f"[SKIP] {sub1.name} — no VCF files found in sub-directories")
        continue

    output_path = sub1 / f"{sub1.name}.vcf"

    print(f"[MERGE] {sub1.name}/")
    print(f"        Output → {output_path.relative_to(root)}")
    print(f"        Found {len(vcf_files)} VCF file(s):")

    variant_rows = []

    for vcf_path in vcf_files:
        print(f"        + {vcf_path.relative_to(sub1)}")
        file_variants = 0

        try:
            with open(vcf_path, 'r') as fh:
                for line in fh:
                    line = line.rstrip('\n')
                    if line.startswith('#') or not line.strip():
                        continue  # skip all header/meta/empty lines
                    variant_rows.append(line)
                    file_variants += 1

            print(f"          → {file_variants} variant(s)")
            total_files += 1

        except Exception as e:
            print(f"          ✗ Error reading {vcf_path.name}: {e}")

    # Write merged output
    try:
        with open(output_path, 'w') as out:
            out.write(VCF_HEADER + '\n')
            for row in variant_rows:
                out.write(row + '\n')

        print(f"        ✓ Written: {len(variant_rows)} total variants\n")
        total_merged   += 1
        total_variants += len(variant_rows)

    except Exception as e:
        print(f"        ✗ Failed to write output: {e}\n")

