"""
AlphaGenome Variant Scorer — Command Line Version
Usage:
    python score_variants.py --input variants.vcf --api_key YOUR_API_KEY [--output results.csv]

Input vcf must have columns: CHROM, POS, REF, ALT
"""
import argparse
import os
import time
import pandas as pd
from tqdm import tqdm
from io import StringIO
from alphagenome.data import genome
from alphagenome.models import dna_client, variant_scorers

# ── CLI arguments ──────────────────────────────────────────────────────────────
parser = argparse.ArgumentParser(description='Score variants with AlphaGenome')
parser.add_argument('--input',    required=True,  help='Path to input vcf files')
parser.add_argument('--output',   default='variant_scores.csv', help='Output CSV path')
parser.add_argument('--api_key',  default=None,   help='AlphaGenome API key (or set ALPHAGENOME_API_KEY env var)')
parser.add_argument('--sleep',    default=0.5, type=float,
                    help='Seconds to sleep between API calls (default: 0.5)')
# Scorer selection flags
parser.add_argument('--score_rna_seq',           action=argparse.BooleanOptionalAction, default=False)
parser.add_argument('--score_cage',              action=argparse.BooleanOptionalAction, default=False)
parser.add_argument('--score_procap',            action=argparse.BooleanOptionalAction, default=True)
parser.add_argument('--score_atac',              action=argparse.BooleanOptionalAction, default=True)
parser.add_argument('--score_dnase',             action=argparse.BooleanOptionalAction, default=True)
parser.add_argument('--score_chip_histone',      action=argparse.BooleanOptionalAction, default=True)
parser.add_argument('--score_chip_tf',           action=argparse.BooleanOptionalAction, default=True)
parser.add_argument('--score_polyadenylation',   action=argparse.BooleanOptionalAction, default=True)
parser.add_argument('--score_splice_sites',      action=argparse.BooleanOptionalAction, default=True)
parser.add_argument('--score_splice_site_usage', action=argparse.BooleanOptionalAction, default=True)
parser.add_argument('--score_splice_junctions',  action=argparse.BooleanOptionalAction, default=True)

args = parser.parse_args()

# ── API key ────────────────────────────────────────────────────────────────────
api_key = args.api_key or os.environ.get('ALPHAGENOME_API_KEY')
if not api_key:
    raise ValueError(
        "Provide an API key via --api_key or the ALPHAGENOME_API_KEY environment variable."
    )

# ── Load model ─────────────────────────────────────────────────────────────────
print("Loading AlphaGenome model...")
dna_model = dna_client.create(api_key)

organism_map = {'human': dna_client.Organism.HOMO_SAPIENS}
organism = organism_map['human']

sequence_length = dna_client.SUPPORTED_SEQUENCE_LENGTHS['SEQUENCE_LENGTH_500KB']

# ── Load variants ──────────────────────────────────────────────────────────────
print(f"Reading variants from: {args.input}")
with open(args.input, 'r') as f:
    lines = [l for l in f if not l.startswith('##')]
df = pd.read_csv(StringIO(''.join(lines)), sep='\t')
df = df.rename(columns={'#CHROM': 'CHROM'})

required_cols = {'CHROM', 'POS', 'REF', 'ALT'}
missing = required_cols - set(df.columns)
if missing:
    raise ValueError(f"Input VCF is missing columns: {missing}")

print(f"Found {len(df)} variants to score.\n")

# ── Scorer selection ───────────────────────────────────────────────────────────
scorer_selections = {
    'rna_seq':            args.score_rna_seq,
    'cage':               args.score_cage,
    'procap':             args.score_procap,
    'atac':               args.score_atac,
    'dnase':              args.score_dnase,
    'chip_histone':       args.score_chip_histone,
    'chip_tf':            args.score_chip_tf,
    'polyadenylation':    args.score_polyadenylation,
    'splice_sites':       args.score_splice_sites,
    'splice_site_usage':  args.score_splice_site_usage,
    'splice_junctions':   args.score_splice_junctions,
}

all_scorers = variant_scorers.RECOMMENDED_VARIANT_SCORERS
selected_scorers = [
    all_scorers[key]
    for key in all_scorers
    if scorer_selections.get(key.lower(), False)
]

print(f"Active scorers: {[k for k, v in scorer_selections.items() if v]}\n")

# ── Score variants ─────────────────────────────────────────────────────────────
# Collect raw rows as a list of dicts — build DataFrame once at the end
all_rows = []

for i, df_row in tqdm(df.iterrows(), total=len(df)):
    variant_id = f"{df_row['CHROM']}:{df_row['POS']}:{df_row['REF']}>{df_row['ALT']}"
    try:
        variant = genome.Variant(
            chromosome=str(df_row['CHROM']),
            position=int(df_row['POS']),
            reference_bases=df_row['REF'],
            alternate_bases=df_row['ALT'],
            name=df_row['ID'],
        )

        interval = variant.reference_interval.resize(sequence_length)

        variant_scores = dna_model.score_variant(
            interval=interval,
            variant=variant,
            variant_scorers=selected_scorers,
            organism=organism,
        )

        # tidy_scores returns a DataFrame — convert each row to a dict
        df_scores = variant_scorers.tidy_scores(variant_scores)

        # Tag every row with the variant info so you know which variant it came from
        df_scores['variant_id']  = variant_id
        df_scores['CHROM']       = df_row['CHROM']
        df_scores['POS']         = df_row['POS']
        df_scores['REF']         = df_row['REF']
        df_scores['ALT']         = df_row['ALT']
        df_scores['ID']          = df_row['ID']

        all_rows.append(df_scores)
        print(f"  ✓ {variant_id} ({len(df_scores)} scores)")

    except Exception as e:
        print(f"  ✗ Failed {variant_id}: {e}")
        # Add a single error row so failures are traceable
        all_rows.append(pd.DataFrame([{
            'variant_id': variant_id,
            'CHROM': df_row['CHROM'],
            'POS':   df_row['POS'],
            'REF':   df_row['REF'],
            'ALT':   df_row['ALT'],
            'ID':    df_row['ID'],
            'error': str(e),
        }]))

    time.sleep(args.sleep)

# ── Build final DataFrame ──────────────────────────────────────────────────────
print(f"\nBuilding final DataFrame...")
df_all = pd.concat(all_rows, ignore_index=True)

# Reorder so variant info columns come first
id_cols = ['variant_id', 'CHROM', 'POS', 'ID', 'REF', 'ALT']
other_cols = [c for c in df_all.columns if c not in id_cols]
df_all = df_all[id_cols + other_cols]

print(df_all)          # preview in terminal
print(f"\nShape: {df_all.shape[0]} rows x {df_all.shape[1]} columns")

# ── Save to CSV ────────────────────────────────────────────────────────────────
df_all.to_csv(args.output, index=False)
print(f"\nSaved to: {args.output}")
