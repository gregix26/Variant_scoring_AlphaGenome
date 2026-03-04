"""
Filter AlphaGenome scores for brain ontologies and extract top N variants per modality.
Usage:
    python filter_top10.py --input variant_scores.csv --output_top top10.csv --top_n 20
"""

import argparse
import pandas as pd

parser = argparse.ArgumentParser(description='Filter brain ontologies and extract top N per modality')
parser.add_argument('--input',      required=True,              help='Path to variant_scores.csv')
parser.add_argument('--output_filtered', default='brain_scores.csv',      help='Brain-filtered scores CSV')
parser.add_argument('--output_top', default='top10_per_modality.csv',     help='Top N per modality CSV')
parser.add_argument('--top_n',      default=10, type=int,       help='Top N variants per modality (default: 10)')
args = parser.parse_args()

# ── Brain ontology terms ───────────────────────────────────────────────────────
BRAIN_ONTOLOGY_TERMS = [
    'UBERON:0000955',   #brain
    'UBERON:0001851',   # cortex
    'UBERON:0002037',   # cerebellum
    'UBERON:0002421',   # hippocampus
    'UBERON:0001876',   # amygdala
    'UBERON:0001897',   # hypothalamus
    'UBERON:0001898',   # thalamus
    'CL:0002590', # brain 
    'CL:2000043', # brain
    'CL:2000044',  # brain
    'CL:0000100', # neuron 
    'CL:0000047',  # neuron 
    'CL:0000540', # neuron 
    'CL:0000679', # neuron 
    'CL:0002603', # astrocyte
    'CL:0002605', # astrocyte
    'CL:0000127', # astrocyte
    'CL:0002606' # astrocyte
    'CL:0002453' # OPC 
]


# ── Load ───────────────────────────────────────────────────────────────────────
print(f"Loading: {args.input}")
df = pd.read_csv(args.input)
print(f"Total rows    : {len(df)}")
print(f"Columns       : {df.columns.tolist()}\n")

# ── Filter by ontology ─────────────────────────────────────────────────────────
# Try ontology column first (most precise)
ontology_col = next(
    (c for c in ['ontology_curie'] if c in df.columns),
    None
)

if ontology_col:
    print(f"Filtering by ontology column: '{ontology_col}'")
    ontology_mask = df[ontology_col].isin(BRAIN_ONTOLOGY_TERMS)
else:
    print("No ontology ID column found.")
    ontology_mask = pd.Series(False, index=df.index)

# Combine: match either ontology ID or keyword
brain_mask = ontology_mask
df_brain = df[brain_mask].copy()

print(f"Rows after brain filter : {len(df_brain)}")
print(f"Rows excluded           : {len(df) - len(df_brain)}\n")

if len(df_brain) == 0:
    print("WARNING: No brain rows found. Check your column names:")
    print(df.columns.tolist())
    exit(1)

# ── Report what was found ──────────────────────────────────────────────────────
if ontology_col:
    print(f"Ontology terms found:")
    print(df_brain[ontology_col].value_counts().to_string())
    print()

# ── Save filtered ──────────────────────────────────────────────────────────────
df_brain.to_csv(args.output_filtered, index=False)
print(f"Brain-filtered scores saved to: {args.output_filtered}\n")

# ── Extract top N per modality ─────────────────────────────────────────────────
score_col = next(
    (c for c in ['quantile_score'] if c in df_brain.columns),
    None
)
modality_col = next(
    (c for c in ['Assay title'] if c in df_brain.columns),
    None
)

if not score_col:
    print(f"WARNING: No score column found. Available: {df_brain.columns.tolist()}")
    exit(1)
if not modality_col:
    print(f"WARNING: No modality column found. Available: {df_brain.columns.tolist()}")
    exit(1)

print(f"Score column    : '{score_col}'")
print(f"Modality column : '{modality_col}'")
print(f"Extracting top {args.top_n} variants per modality by absolute score...\n")

df_top = (
    df_brain
    .dropna(subset=[score_col])
    .assign(abs_score=lambda x: x[score_col].abs())   # capture both up and down effects
    .sort_values('abs_score', ascending=False)
    .groupby(modality_col)
    .head(args.top_n)
    .drop(columns='abs_score')
    .reset_index(drop=True)
)

print(f"Top {args.top_n} per modality — variant counts:")
print(df_top[modality_col].value_counts().to_string())

df_top.to_csv(args.output_top, index=False)
print(f"\nTop {args.top_n} per modality saved to: {args.output_top}")
