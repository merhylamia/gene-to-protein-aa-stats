import os
from collections import Counter
import matplotlib.pyplot as plt

# ---------- Paths ----------
IN = "inputs"
GENE_FP = os.path.join(IN, "gene_sequence.txt")
EXON_FP = os.path.join(IN, "exon_positions.txt")
CODE_FP = os.path.join(IN, "code.txt")
FIG_DIR = "figures"
os.makedirs(FIG_DIR, exist_ok=True)

# ---------- Helpers ----------
def read_gene(path):
    with open(path, "r") as f:
        return "".join([line.strip() for line in f])

def read_exons(path):
    exons = []
    with open(path, "r") as f:
        for line in f:
            a, b = line.split()[:2]
            exons.append((int(a), int(b)))   # assume 0-based, end-exclusive
    return exons

def read_codon_table(path):
    table = {}
    with open(path, "r") as f:
        for line in f:
            # first two tokens: CODON  AA1
            toks = line.split()
            if len(toks) >= 2:
                table[toks[0]] = toks[1]
    return table

def build_mrna(dna_seq, exons):
    # concat exons, replace T->U
    pieces = [dna_seq[s:e] for (s, e) in exons]
    return "".join(pieces).replace("T", "U")

def translate(mrna, codon_table):
    start = "AUG"
    stops = {"UAA", "UAG", "UGA"}
    # find first AUG
    i0 = mrna.find(start)
    if i0 == -1: return ""
    aa = []
    for i in range(i0, len(mrna), 3):
        codon = mrna[i:i+3]
        if len(codon) < 3: break
        if codon in stops: break
        aa1 = codon_table.get(codon, "")
        if aa1 == "X": break
        if aa1:
            aa.append(aa1)
    return "".join(aa)

# ---------- Load inputs ----------
dna = read_gene(GENE_FP)
exons = read_exons(EXON_FP)
codon = read_codon_table(CODE_FP)

# ---------- Part 1: Translation ----------
mrna = build_mrna(dna, exons)
protein = translate(mrna, codon)
print(f"Protein length (original): {len(protein)}")

# ---------- Part 2: AA stats ----------
AA_ALL = list("ACDEFGHIKLMNPQRSTVWY")
CATEGORIES = {
    "Positively charged": list("RHK"),
    "Negatively charged": list("DE"),
    "Polar": list("NCQSTY"),
    "Non-polar": list("AILMFPWV G".replace(" ", "")),
}

def count_all(seq):
    cnt = Counter(seq)
    freqs = [cnt.get(a, 0) for a in AA_ALL]
    total = sum(freqs) or 1
    perc = [f * 100.0 / total for f in freqs]
    return freqs, perc

def count_per_category(seq):
    cnt = Counter(seq)
    out = {}
    total = len(seq) or 1
    for cat, aas in CATEGORIES.items():
        s = sum(cnt.get(a, 0) for a in aas)
        out[cat] = (s, s * 100.0 / total)
    return out

def count_within_category(seq, cat_name):
    cnt = Counter(seq)
    aas = CATEGORIES[cat_name]
    total = sum(cnt.get(a, 0) for a in aas) or 1
    rows = []
    for a in aas:
        f = cnt.get(a, 0)
        p = f * 100.0 / total
        rows.append((a, f, p))
    rows.sort(key=lambda x: x[1], reverse=True)
    return rows

def plot_bar_line(labels, freq, perc, title, fig_path):
    x = list(range(len(labels)))
    fig, ax1 = plt.subplots()
    ax2 = ax1.twinx()
    ax1.bar(x, freq, alpha=0.7)
    ax2.plot(x, perc, "o--")
    ax1.set_xticks(x)
    ax1.set_xticklabels(labels, rotation=0)
    ax1.set_ylabel("Frequency")
    ax2.set_ylabel("Percentage")
    ax1.set_title(title)
    fig.tight_layout()
    plt.savefig(fig_path, dpi=300)
    plt.close(fig)

print("\nMenu:\n 1) All\n 2) Per category\n 3) Within category\n 4) Specific AA")
choice = input("> ").strip()

if choice in {"1", "All", "all"}:
    freqs, perc = count_all(protein)
    # sort by freq desc
    rows = sorted(zip(AA_ALL, freqs, perc), key=lambda x: x[1], reverse=True)
    for a, f, p in rows:
        print(f"{a}\t{f}\t{p:.2f}%")
    labels = [r[0] for r in rows]
    fvals = [r[1] for r in rows]
    pvals = [r[2] for r in rows]
    plot_bar_line(labels, fvals, pvals, "AA Composition (All)", os.path.join(FIG_DIR, "aa_all.png"))

elif choice in {"2", "Per category", "per"}:
    cats = count_per_category(protein)
    rows = sorted(cats.items(), key=lambda kv: kv[1][0], reverse=True)
    for cat, (f, p) in rows:
        print(f"{cat}\t{f}\t{p:.2f}%")
    labels = [r[0] for r in rows]
    fvals = [r[1][0] for r in rows]
    pvals = [r[1][1] for r in rows]
    plot_bar_line(labels, fvals, pvals, "AA Composition (By Category)", os.path.join(FIG_DIR, "aa_by_category.png"))

elif choice in {"3", "Within category", "within"}:
    print("Categories:", ", ".join(CATEGORIES.keys()))
    cat = input("Pick a category exactly as shown: ").strip()
    rows = count_within_category(protein, cat)
    for a, f, p in rows:
        print(f"{a}\t{f}\t{p:.2f}%")
    labels = [r[0] for r in rows]
    fvals = [r[1] for r in rows]
    pvals = [r[2] for r in rows]
    plot_bar_line(labels, fvals, pvals, f"AA Within Category: {cat}", os.path.join(FIG_DIR, f"aa_within_{cat.replace(' ','_')}.png"))

elif choice in {"4", "Specific AA", "specific"}:
    a = input("One-letter AA code (A,C,D,...,Y): ").strip().upper()
    if a not in AA_ALL:
        print("Invalid AA code.")
    else:
        cnt = Counter(protein)
        f = cnt.get(a, 0)
        p = f * 100.0 / (len(protein) or 1)
        print(f"{a}\t{f}\t{p:.2f}%")

# ---------- Part 3: Mutation ----------
print("\n--- Mutation (T→A at DNA index 30049, 0-based) ---")
dna_mut = list(dna)
if 0 <= 30049 < len(dna_mut):
    dna_mut[30049] = "A"
dna_mut = "".join(dna_mut)
mrna_mut = build_mrna(dna_mut, exons)
protein_mut = translate(mrna_mut, codon)
print(f"Protein length (mutated):   {len(protein_mut)}")
print(f"Length difference:          {len(protein) - len(protein_mut)} aa (positive means truncation)")
