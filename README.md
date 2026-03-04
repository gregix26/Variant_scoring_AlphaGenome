## Variant_scoring_AlphaGenome

# Install AlphaGenome
Follow installation instruction from AlphaGenome: https://www.alphagenomedocs.com/installation.html
Request API key. 
It is recommended to create a separate environment with Python=3.11

# Prepare your variants for batch processing
Use merge_vcfs.py to merge separate VCF files of variants into one big VCF file that the main script can loop through. ALphaGenome works with variant position, alternative and reference bases.

# Run AlphaGenome
Use run_alpgagenome.py. Specify modality outputs (all tracks or subset). Output file is csv with all score. Takes about 10 mins to run 50 variants with all tracks enabled.

# Visualization
