"""
AlphaGenome Variant Visualiser — Fully Annotated Version
"""

import argparse   
import os         
import io         
import time       
import pandas as pd  
import matplotlib
matplotlib.use('Agg') # 'Agg' = Anti-Grain Geometry — renders to file instead of a screen window
import matplotlib.pyplot as plt  
from matplotlib.backends.backend_pdf import PdfPages 
from io import StringIO 
from tqdm import tqdm 

from alphagenome.data import gene_annotation, genome, transcript
# gene_annotation — functions to filter the GTF gene table (e.g. protein-coding only)
# genome         — Variant class: holds chrom/pos/ref/alt
# transcript     — TranscriptExtractor: finds transcripts overlapping a region

from alphagenome.models import dna_client
# dna_client — connects to the AlphaGenome API, holds Organism enum and OutputType enum

from alphagenome.visualization import plot_components
# plot_components — all the plotting classes: OverlaidTracks, Sashimi, TranscriptAnnotation etc.

# ── CLI arguments ──────────────────────────────────────────────────────────────
parser = argparse.ArgumentParser(description='Visualise AlphaGenome variant effects')
parser.add_argument('--input',       required=True,           help='Top variants CSV')
parser.add_argument('--api_key',     default=None,            help='AlphaGenome API key')
parser.add_argument('--output_dir',  default='variant_plots', help='Directory to save PDFs')
parser.add_argument('--sleep',       default=1.0, type=float, help='Sleep between API calls')
parser.add_argument('--sequence_length', default='500KB',
                    choices=['2KB', '16KB', '100KB', '500KB', '1MB'],
                    help='Sequence window around variant')
parser.add_argument('--plot_width',  default=43008, type=int, help='Width of plot window in bp')

parser.add_argument('--plot_rna_seq',           action=argparse.BooleanOptionalAction, default=True)
parser.add_argument('--plot_cage',              action=argparse.BooleanOptionalAction, default=True)
parser.add_argument('--plot_atac',              action=argparse.BooleanOptionalAction, default=True)
parser.add_argument('--plot_dnase',             action=argparse.BooleanOptionalAction, default=True)
parser.add_argument('--plot_chip_histone',      action=argparse.BooleanOptionalAction, default=True)
parser.add_argument('--plot_chip_tf',           action=argparse.BooleanOptionalAction, default=True)
parser.add_argument('--plot_splice_sites',      action=argparse.BooleanOptionalAction, default=True)
parser.add_argument('--plot_splice_site_usage', action=argparse.BooleanOptionalAction, default=True)
parser.add_argument('--plot_splice_junctions',  action=argparse.BooleanOptionalAction, default=True)
parser.add_argument('--plot_contact_maps',      action=argparse.BooleanOptionalAction, default=True)

# store_true = flag presence sets it True, absence leaves it False (no --no- version needed)
parser.add_argument('--positive_strand_only', action='store_true', default=False)
parser.add_argument('--negative_strand_only', action='store_true', default=False)

args = parser.parse_args()


# ── API key ────────────────────────────────────────────────────────────────────
api_key = args.api_key or os.environ.get('ALPHAGENOME_API_KEY')
if not api_key:
    raise ValueError("Provide --api_key or set ALPHAGENOME_API_KEY env var.")

os.makedirs(args.output_dir, exist_ok=True)

# ── Load model ─────────────────────────────────────────────────────────────────
print("Loading AlphaGenome model...")
dna_model = dna_client.create(api_key)

# Set organism — HOMO_SAPIENS uses hg38 reference genome
organism = dna_client.Organism.HOMO_SAPIENS

sequence_length = dna_client.SUPPORTED_SEQUENCE_LENGTHS[f'SEQUENCE_LENGTH_{args.sequence_length}']

# ── Load GTF gene annotation ───────────────────────────────────────────────────
# GTF = Gene Transfer Format — contains gene/transcript coordinates for hg38
# This feather file is a fast binary version of the GTF hosted by Google
HG38_GTF_FEATHER = (
    'https://storage.googleapis.com/alphagenome/reference/gencode/'
    'hg38/gencode.v46.annotation.gtf.gz.feather'
)

print("Loading gene annotation (hg38)...")
gtf = pd.read_feather(HG38_GTF_FEATHER)  # load entire annotation table into DataFrame

# Keep only protein-coding genes (removes lncRNA, pseudogenes etc.)
# Then keep only transcripts with support level '1' (highest confidence)
gtf_transcript = gene_annotation.filter_transcript_support_level(
    gene_annotation.filter_protein_coding(gtf), ['1']
)

# TranscriptExtractor: given a genomic interval, returns all transcripts overlapping it
# Used for drawing the gene track at the top of the plot
transcript_extractor = transcript.TranscriptExtractor(gtf_transcript)

# A second extractor that only keeps the longest transcript per gene
# Avoids cluttered plots with many isoforms
longest_transcript_extractor = transcript.TranscriptExtractor(
    gene_annotation.filter_to_longest_transcript(gtf_transcript)
)

# ── Requested outputs ──────────────────────────────────────────────────────────
requested_outputs = [*dna_client.OutputType]

# ── Load variants CSV ──────────────────────────────────────────────────────────
print(f"Loading variants from: {args.input}")
df = pd.read_csv(args.input) 
print(f"Found {len(df)} variants to visualise.\n")

ontology_col = next(
    (c for c in ['ontology_curie'] if c in df.columns),
    None
)

# ── Prediction cache ───────────────────────────────────────────────────────────
# If the same variant + interval + ontology is requested twice, return cached result
# Saves API calls when the same variant appears in multiple modality rows
_prediction_cache = {}

def predict_variant_cached(interval, variant, ontology_terms):
    # Build a hashable key from all inputs that affect the result
    cache_key = (str(interval), str(variant), tuple(ontology_terms))

    # Return cached result if we've seen this exact combination before
    if cache_key in _prediction_cache:
        return _prediction_cache[cache_key]

    # Otherwise call the API — returns REF and ALT track predictions
    result = dna_model.predict_variant(
        interval=interval,
        variant=variant,
        organism=organism,
        requested_outputs=requested_outputs,
        ontology_terms=ontology_terms,
    )
    _prediction_cache[cache_key] = result  # store for reuse
    return result

# ── Colour scheme for REF vs ALT tracks ───────────────────────────────────────
# These colours are passed to OverlaidTracks so REF and ALT are visually distinct
ref_alt_colors = {'REF': 'dimgrey', 'ALT': 'red'}

# ── Main loop — one PDF per variant ───────────────────────────────────────────
for i, df_row in tqdm(df.iterrows(), total=len(df)):
    chrom = str(df_row.get('CHROM', df_row.get('chromosome', '')))
    pos   = int(df_row.get('POS',   df_row.get('position', 0)))
    ref   = str(df_row.get('REF',   df_row.get('reference_bases', '')))
    alt   = str(df_row.get('ALT',   df_row.get('alternate_bases', '')))
    rsid  = str(df_row.get('ID',    df_row.get('name', 'unknown')))

    # Human-readable ID used in titles and filenames
    variant_id = f"{chrom}_{pos}_{ref}_{alt}"
    safe_id    = variant_id.replace(':', '_').replace('>', '_')

    # Use per-row ontology if available, otherwise fall back to generic brain term
    if ontology_col and pd.notna(df_row.get(ontology_col)):
        ontology_terms = [df_row[ontology_col]]
    else:
        ontology_terms = ['UBERON:0000955']   # brain (generic fallback)

    # Full path for this variant's PDF output file
    pdf_path = os.path.join(args.output_dir, f"{safe_id}.pdf")
    print(f"\n[{i+1}/{len(df)}] {variant_id} | ontology: {ontology_terms}")

    try:
        # Create a Variant object — AlphaGenome's representation of a SNP/indel
        variant = genome.Variant(
            chromosome=chrom,
            position=pos,
            reference_bases=ref,
            alternate_bases=alt,
            name=rsid,
        )

        # reference_interval = a genomic interval centred on the variant position
        # .resize() expands it to the requested window (e.g. 500KB)
        interval = variant.reference_interval.resize(sequence_length)

        # Call predict_variant (or return cached result)
        # Returns an object with .reference and .alternate — both full track outputs
        output = predict_variant_cached(interval, variant, ontology_terms)

        # Unpack REF and ALT track data
        ref_out, alt_out = output.reference, output.alternate

        # Optionally restrict to one DNA strand
        # filter_to_strand removes tracks from the opposite strand
        if args.positive_strand_only:
            ref_out = ref_out.filter_to_strand(strand='+')
            alt_out = alt_out.filter_to_strand(strand='+')
        elif args.negative_strand_only:
            ref_out = ref_out.filter_to_strand(strand='-')
            alt_out = alt_out.filter_to_strand(strand='-')

        # ── Build list of plot components ──────────────────────────────────────
        # Each component is one "panel" in the stacked genome browser plot
        components = []

        # First panel: transcript annotation track (gene arrows, exon boxes)
        transcripts = longest_transcript_extractor.extract(interval)
        components.append(plot_components.TranscriptAnnotation(transcripts))

        # Map each flag name → (ref_data, alt_data, output_type_name, y-axis label template)
        # ylabel_template uses {biosample_name}, {strand} etc. as placeholders
        # that AlphaGenome fills in from track metadata
        plot_map = {
            'plot_rna_seq':        (ref_out.rna_seq,        alt_out.rna_seq,        'RNA_SEQ',        '{output_type}: {biosample_name} ({strand})\n{name}'),
            'plot_cage':           (ref_out.cage,           alt_out.cage,           'CAGE',           '{output_type}: {biosample_name} ({strand})\n{name}'),
            'plot_atac':           (ref_out.atac,           alt_out.atac,           'ATAC',           '{output_type}: {biosample_name} ({strand})\n{name}'),
            'plot_dnase':          (ref_out.dnase,          alt_out.dnase,          'DNASE',          '{output_type}: {biosample_name} ({strand})\n{name}'),
            'plot_chip_histone':   (ref_out.chip_histone,   alt_out.chip_histone,   'CHIP_HISTONE',   '{output_type}: {biosample_name} ({strand})\n{histone_mark}'),
            'plot_chip_tf':        (ref_out.chip_tf,        alt_out.chip_tf,        'CHIP_TF',        '{output_type}: {biosample_name} ({strand})\n{transcription_factor}'),
            'plot_splice_sites':   (ref_out.splice_sites,   alt_out.splice_sites,   'SPLICE_SITES',   '{output_type}: {name} ({strand})'),
            'plot_splice_site_usage': (ref_out.splice_site_usage, alt_out.splice_site_usage, 'SPLICE_SITE_USAGE', '{output_type}: {biosample_name} ({strand})\n{name}'),
            'plot_contact_maps':   (ref_out.contact_maps,   alt_out.contact_maps,   'CONTACT_MAPS',   '{output_type}: {biosample_name} ({strand})'),
            'plot_splice_junctions': (ref_out.splice_junctions, alt_out.splice_junctions, 'SPLICE_JUNCTIONS', '{output_type}: {biosample_name} ({strand})\n{name}'),
        }

        for flag, (ref_data, alt_data, output_type, ylabel_template) in plot_map.items():
            if not getattr(args, flag):
                continue  # user turned this track type off — skip

            # None means AlphaGenome didn't return this output type
            if ref_data is None or alt_data is None:
                continue

            # shape[-1] is the number of tracks — 0 means no data for this ontology
            if ref_data.values.shape[-1] == 0:
                print(f"  ⚠ No tracks for {output_type} — skipping")
                continue

            # Replace {output_type} placeholder in label template
            ylabel = ylabel_template.replace('{output_type}', output_type)

            if output_type == 'CONTACT_MAPS':
                # Contact maps show the DIFFERENCE (ALT - REF) as a heatmap
                components.append(plot_components.ContactMapsDiff(
                    tdata=alt_data - ref_data,
                    ylabel_template=ylabel,
                ))
            elif output_type == 'SPLICE_JUNCTIONS':
                # Sashimi plots: arc diagram showing splice junction usage
                # REF and ALT shown as separate panels (can't overlay arcs)
                components.append(plot_components.Sashimi(ref_data, ylabel_template='REF: ' + ylabel))
                components.append(plot_components.Sashimi(alt_data, ylabel_template='ALT: ' + ylabel))
            else:
                # OverlaidTracks: REF and ALT drawn on same axes in different colours
                components.append(plot_components.OverlaidTracks(
                    tdata={'REF': ref_data, 'ALT': alt_data},
                    colors=ref_alt_colors,
                    ylabel_template=ylabel,
                ))

        # ── Render the stacked plot ────────────────────────────────────────────
        plot_interval = interval.resize(args.plot_width)

        fig = plot_components.plot(
            components=components,        # list of track panels built above
            interval=plot_interval,       # genomic region to display on x-axis
            annotations=[
                plot_components.VariantAnnotation([variant]),  # vertical line at variant position
            ],
        )

        # ── Save to PDF ────────────────────────────────────────────────────────
        # PdfPages context manager: everything saved inside becomes pages of one PDF
        with PdfPages(pdf_path) as pdf:

            # Page 1: title page — plain figure with text only
            title_fig, ax = plt.subplots(figsize=(12, 2))
            ax.axis('off')   # hide axes/ticks — we just want text on white background
            ax.text(0.5, 0.6, f'Variant: {variant_id}',
                    ha='center', va='center', fontsize=16, fontweight='bold',
                    transform=ax.transAxes)   # transform=ax.transAxes means 0-1 coordinates
            ax.text(0.5, 0.3, f'rsID: {rsid}  |  Ontology: {ontology_terms[0]}',
                    ha='center', va='center', fontsize=11, color='gray',
                    transform=ax.transAxes)
            pdf.savefig(title_fig, bbox_inches='tight')  # save title as page 1
            plt.close(title_fig)   # free memory — important in long loops

            # Page 2: the actual genome browser plot
            pdf.savefig(fig, bbox_inches='tight')
            plt.close(fig)   # free memory

        print(f"  ✓ Saved → {pdf_path}")

    except Exception as e:
        # Catch any failure so one bad variant doesn't stop the whole run
        print(f"  ✗ Failed {variant_id}: {e}")

    # Polite pause between API calls to avoid rate limiting
    time.sleep(args.sleep)

print(f"\nDone! PDFs saved to: {args.output_dir}/")
