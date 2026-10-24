# 🧬🎵 Codon Melody: DNA & FASTA Sequence Sonification

![Python Version](https://img.shields.io/badge/python-3.7%2B-blue)
![Dependencies](https://img.shields.io/badge/dependencies-none%20(stdlib%20only)-brightgreen)
![License](https://img.shields.io/badge/license-MIT-green)

**Codon Melody** is a lightweight, zero-dependency Python toolkit that translates DNA and RNA coding sequences into playable 16-bit WAV audio melodies. By mapping the 64 genetic codons onto a 3-octave C Major scale, this project lets you *hear* the primary structure of genes—and audibly detect single-nucleotide point mutations.

---

## ✨ Key Features

- **Zero External Dependencies:** Built entirely with Python's standard library (`wave`, `struct`, `math`, `itertools`, `argparse`). No `pip install` required.
- **Sickle-Cell Mutation Demo (`gene_to_melody.py`):** Translates the 444-bp human beta-globin (*HBB*) coding sequence (148 codons) into a WAV file, introduces the classic sickle-cell point mutation at base position 20 (`A` $\rightarrow$ `T`), and compares the two sequences across codons, amino acids, and musical pitches.
- **Custom FASTA CLI Tool (`fasta_to_music.py`):** Converts any single- or multi-record `.fasta` file into `.wav` audio from the terminal, complete with customizable tempo, automatic RNA (`U` $\rightarrow$ `T`) conversion, and graceful handling of ambiguous bases (`N` $\rightarrow$ musical rest).
- **Warm Audio Synthesis:** Synthesizes each note using a fundamental sine wave blended with a subtle second harmonic and a 15 ms linear attack/release envelope to prevent audio clicking.

---

## 🎼 How the Musical Mapping Works

1. **Lexicographical Codon Sorting:** All 64 standard codons (`AAA`, `AAC`, `AAG`, `AAT`, ..., `TTT`) are sorted alphabetically from index `0` to `63`.
2. **3-Octave C Major Scale Assignment:** Each codon is mapped cyclically (`index % 21`) onto a **21-note C Major scale** spanning three octaves (`C3` at **130.81 Hz** through `B5` at **987.77 Hz**).

### Why a 21-Note Scale Matters for Mutation Detection
In alphabetical order, the wild-type 7th codon of *HBB* (`GAG`, coding for Glutamate) sits at index `34`, while the sickle-cell mutant codon (`GTG`, coding for Valine) sits at index `46`—a separation of exactly 12 positions.
- If mapped onto a standard 12-note scale (`index % 12`), both `GAG` and `GTG` would alias onto the exact same musical pitch.
- By cycling across the **21-note C Major scale** (`index % 21`), `GAG` (`34 % 21 = 13`) maps to **B4 (493.88 Hz)** and `GTG` (`46 % 21 = 4`) maps to **G3 (196.00 Hz)**, creating a dramatic drop of more than an octave on the 7th note of the melody.

---

## 📂 Repository Structure

```text
codon-melody/
├── README.md                  # Project documentation
├── LICENSE                    # MIT License
├── .gitignore                 # Git ignore rules
├── gene_to_melody.py          # Main HBB wild-type vs. A20T mutation comparison script
├── fasta_to_music.py          # Command-line converter for custom FASTA files
└── examples/
    └── hbb_sample.fasta       # Sample FASTA file (wild-type & mutant HBB)
```

---

## 🚀 Getting Started

### Prerequisites
- **Python 3.7+** (No external packages required)

### Clone the Repository
```bash
git clone [https://github.com/](https://github.com/)<your-username>/codon-melody.git
cd codon-melody
```

---

## 💻 Usage

### 1. Run the *HBB* Sickle-Cell Mutation Demo (`gene_to_melody.py`)

Run the main script to print the full 64-codon mapping table, generate both the wild-type and mutant WAV files, and verify that exactly **1 codon**, **1 amino acid**, and **1 musical note** differ between the two sequences:

```bash
python gene_to_melody.py
```

**Generated Files:**
- `hbb_wildtype.wav` (~26.6 seconds, 148 notes)
- `hbb_mutant_A20T.wav` (~26.6 seconds, 148 notes)

**Sample Console Output:**
```text
Saved wild-type melody to: hbb_wildtype.wav
Saved mutant melody to:    hbb_mutant_A20T.wav

=== COMPARISON SUMMARY (Wild-Type vs. A20T Mutant) ===
Codon 7 change:       GAG -> GTG
Amino acid 7 change:  E (Glu) -> V (Val)
Note 7 change:        B4 (493.88 Hz) -> G3 (196.0 Hz)
------------------------------------------------------
Codons that differ:      1
Amino acids that differ: 1
Notes that differ:       1
```

---

### 2. Convert Any Custom FASTA File (`fasta_to_music.py`)

Use `fasta_to_music.py` to convert your own `.fasta`, `.fa`, or `.fna` files into playable WAV melodies:

```bash
python fasta_to_music.py my_sequence.fasta
```

#### Command-Line Options
| Flag | Description | Default |
| :--- | :--- | :--- |
| `fasta_file` | Path to the input FASTA file (single or multi-sequence) | *Required* |
| `-o`, `--output` | Custom output WAV filename | `<fasta_filename>.wav` |
| `-d`, `--duration` | Duration of each codon note in seconds | `0.18` |
| `--show-mapping` | Print the 64-codon to note mapping table before synthesis | `False` |

#### Examples
```bash
# Convert the included sample FASTA (creates two WAV files, one per record)
python fasta_to_music.py examples/hbb_sample.fasta

# Specify a custom output filename and faster tempo (0.12s per codon)
python fasta_to_music.py my_sequence.fasta -o custom_melody.wav -d 0.12

# Print the full 64-codon scale table while converting
python fasta_to_music.py my_sequence.fasta --show-mapping
```

---

## 🎹 Complete 64-Codon Mapping Table

| Codon (AA) $\rightarrow$ Note | Codon (AA) $\rightarrow$ Note | Codon (AA) $\rightarrow$ Note | Codon (AA) $\rightarrow$ Note |
| :--- | :--- | :--- | :--- |
| `AAA` (K) $\rightarrow$ **C3** (130.81 Hz) | `AAC` (N) $\rightarrow$ **D3** (146.83 Hz) | `AAG` (K) $\rightarrow$ **E3** (164.81 Hz) | `AAT` (N) $\rightarrow$ **F3** (174.61 Hz) |
| `ACA` (T) $\rightarrow$ **G3** (196.00 Hz) | `ACC` (T) $\rightarrow$ **A3** (220.00 Hz) | `ACG` (T) $\rightarrow$ **B3** (246.94 Hz) | `ACT` (T) $\rightarrow$ **C4** (261.63 Hz) |
| `AGA` (R) $\rightarrow$ **D4** (293.66 Hz) | `AGC` (S) $\rightarrow$ **E4** (329.63 Hz) | `AGG` (R) $\rightarrow$ **F4** (349.23 Hz) | `AGT` (S) $\rightarrow$ **G4** (392.00 Hz) |
| `ATA` (I) $\rightarrow$ **A4** (440.00 Hz) | `ATC` (I) $\rightarrow$ **B4** (493.88 Hz) | `ATG` (M) $\rightarrow$ **C5** (523.25 Hz) | `ATT` (I) $\rightarrow$ **D5** (587.33 Hz) |
| `CAA` (Q) $\rightarrow$ **E5** (659.26 Hz) | `CAC` (H) $\rightarrow$ **F5** (698.46 Hz) | `CAG` (Q) $\rightarrow$ **G5** (783.99 Hz) | `CAT` (H) $\rightarrow$ **A5** (880.00 Hz) |
| `CCA` (P) $\rightarrow$ **B5** (987.77 Hz) | `CCC` (P) $\rightarrow$ **C3** (130.81 Hz) | `CCG` (P) $\rightarrow$ **D3** (146.83 Hz) | `CCT` (P) $\rightarrow$ **E3** (164.81 Hz) |
| `CGA` (R) $\rightarrow$ **F3** (174.61 Hz) | `CGC` (R) $\rightarrow$ **G3** (196.00 Hz) | `CGG` (R) $\rightarrow$ **A3** (220.00 Hz) | `CGT` (R) $\rightarrow$ **B3** (246.94 Hz) |
| `CTA` (L) $\rightarrow$ **C4** (261.63 Hz) | `CTC` (L) $\rightarrow$ **D4** (293.66 Hz) | `CTG` (L) $\rightarrow$ **E4** (329.63 Hz) | `CTT` (L) $\rightarrow$ **F4** (349.23 Hz) |
| `GAA` (E) $\rightarrow$ **G4** (392.00 Hz) | `GAC` (D) $\rightarrow$ **A4** (440.00 Hz) | `GAG` (E) $\rightarrow$ **B4** (493.88 Hz) | `GAT` (D) $\rightarrow$ **C5** (523.25 Hz) |
| `GCA` (A) $\rightarrow$ **D5** (587.33 Hz) | `GCC` (A) $\rightarrow$ **E5** (659.26 Hz) | `GCG` (A) $\rightarrow$ **F5** (698.46 Hz) | `GCT` (A) $\rightarrow$ **G5** (783.99 Hz) |
| `GGA` (G) $\rightarrow$ **A5** (880.00 Hz) | `GGC` (G) $\rightarrow$ **B5** (987.77 Hz) | `GGG` (G) $\rightarrow$ **C3** (130.81 Hz) | `GGT` (G) $\rightarrow$ **D3** (146.83 Hz) |
| `GTA` (V) $\rightarrow$ **E3** (164.81 Hz) | `GTC` (V) $\rightarrow$ **F3** (174.61 Hz) | `GTG` (V) $\rightarrow$ **G3** (196.00 Hz) | `GTT` (V) $\rightarrow$ **A3** (220.00 Hz) |
| `TAA` (*) $\rightarrow$ **B3** (246.94 Hz) | `TAC` (Y) $\rightarrow$ **C4** (261.63 Hz) | `TAG` (*) $\rightarrow$ **D4** (293.66 Hz) | `TAT` (Y) $\rightarrow$ **E4** (329.63 Hz) |
| `TCA` (S) $\rightarrow$ **F4** (349.23 Hz) | `TCC` (S) $\rightarrow$ **G4** (392.00 Hz) | `TCG` (S) $\rightarrow$ **A4** (440.00 Hz) | `TCT` (S) $\rightarrow$ **B4** (493.88 Hz) |
| `TGA` (*) $\rightarrow$ **C5** (523.25 Hz) | `TGC` (C) $\rightarrow$ **D5** (587.33 Hz) | `TGG` (W) $\rightarrow$ **E5** (659.26 Hz) | `TGT` (C) $\rightarrow$ **F5** (698.46 Hz) |
| `TTA` (L) $\rightarrow$ **G5** (783.99 Hz) | `TTC` (F) $\rightarrow$ **A5** (880.00 Hz) | `TTG` (L) $\rightarrow$ **B5** (987.77 Hz) | `TTT` (F) $\rightarrow$ **C3** (130.81 Hz) |

---

## 📄 License

This project is licensed under the [MIT License](LICENSE).
