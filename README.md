# Gene → mRNA → Protein: Translation, AA Composition, and Mutation Effect (Python)

This project reproduces a classic gene-expression pipeline:
1) **Build mRNA** from DNA using exon positions,
2) **Translate** mRNA to a protein (AUG → stop),
3) **Analyze** amino-acid composition (overall, by category, within category, specific AA),
4) **Mutate** the DNA at position 30049 (0-based, T→A), re-translate, and compare protein lengths.

## Inputs
Place these in `inputs/`:
- `gene_sequence.txt` — full DNA sequence (A/T/G/C)
- `exon_positions.txt` — exon start/stop (0-based, end exclusive)
- `code.txt` — mRNA codon → 1-letter AA (stop codons mapped to `X`)

## How to run
- **Notebook:** open `notebook.ipynb` and run all cells, or
- **Script:** `python main.py`

The script:
- Builds `mRNA` from exons (replaces T→U),
- Translates to `protein` (starts at first AUG, stops at UAA/UAG/UGA),
- Shows AA stats interactively (ALL / Per category / Within category / Specific AA),
- Saves plots to `figures/`,
- Applies the **T→A** mutation at index `30049`, re-translates, and prints both protein lengths.

## Notes
- Exon coordinates are assumed **0-based, end-exclusive**.
- Stop codons in `code.txt` must be mapped to `X`.
- Figures are saved to `figures/`.

## Keywords
Bioinformatics · Gene Expression · Translation · Amino-Acid Composition · Mutation Analysis · Python
