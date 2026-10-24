import itertools
import math
import struct
import wave

# 1. Input coding sequence (Human Beta-Globin / HBB, 444 bases = 148 codons)
WT_SEQ = (
    "ATGGTGCATCTGACTCCTGAGGAGAAGTCTGCCGTTACTGCCCTGTGGGGCAAGGTGAAC"
    "GTGGATGAAGTTGGTGGTGAGGCCCTGGGCAGGCTGCTGGTGGTCTACCCTTGGACCCAG"
    "AGGTTCTTTGAGTCCTTTGGGGATCTGTCCACTCCTGATGCTGTTATGGGCAACCCTAAG"
    "GTGAAGGCTCATGGCAAGAAAGTGCTCGGTGCCTTTAGTGATGGCCTGGCTCACCTGGAC"
    "AACCTCAAGGGCACCTTTGCCACACTGAGTGAGCTGCACTGTGACAAGCTGCACGTGGAT"
    "CCTGAGAACTTCAGGCTCCTGGGCAACGTGCTGGTCTGTGTGCTGGCCCATCACTTTGGC"
    "AAAGAATTCACCCCACCAGTGCAGGCTGCCTATCAGAAAGTGGTGGCTGGTGTGGCTAAT"
    "GCCCTGGCCCACAAGTATCACTAA"
)

# 2. Build a 3-octave C Major scale (21 notes: C3 to B5)
NOTE_NAMES = ["C", "D", "E", "F", "G", "A", "B"]
SEMITONE_OFFSETS = [0, 2, 4, 5, 7, 9, 11]  # Major scale intervals in semitones

SCALE_NOTES = []
for octave in (3, 4, 5):
    for name, offset in zip(NOTE_NAMES, SEMITONE_OFFSETS):
        midi_num = 12 * (octave + 1) + offset  # C4 = MIDI 60, A4 = MIDI 69 (440 Hz)
        freq = 440.0 * (2.0 ** ((midi_num - 69) / 12.0))
        SCALE_NOTES.append((f"{name}{octave}", round(freq, 2)))

# 3. Sort all 64 codons alphabetically (AAA .. TTT) and map onto the 21-note scale
BASES = ["A", "C", "G", "T"]
ALL_CODONS = ["".join(c) for c in itertools.product(BASES, repeat=3)]
CODON_TO_NOTE = {
    codon: SCALE_NOTES[i % len(SCALE_NOTES)]
    for i, codon in enumerate(ALL_CODONS)
}

# Standard genetic code for amino acid translation
CODON_TO_AA = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}


def split_codons(seq):
    """Splits a DNA sequence string into a list of 3-base codons."""
    return [seq[i:i + 3] for i in range(0, len(seq), 3)]


def write_melody_wav(filename, notes, note_duration=0.18, sample_rate=44100):
    """Synthesizes a list of (note_name, freq_hz) tuples into a 16-bit mono WAV file."""
    n_samples_per_note = int(sample_rate * note_duration)
    fade_samples = int(sample_rate * 0.015)  # 15 ms fade-in/out to prevent clicks
    frames = bytearray()

    for _, freq in notes:
        for s in range(n_samples_per_note):
            t = s / sample_rate
            # Fundamental sine wave + 2nd harmonic for a warm, clear tone
            wave_val = (
                0.8 * math.sin(2.0 * math.pi * freq * t)
                + 0.2 * math.sin(4.0 * math.pi * freq * t)
            )

            # Linear attack and release envelope
            if s < fade_samples:
                env = s / fade_samples
            elif s > n_samples_per_note - fade_samples:
                env = (n_samples_per_note - s) / fade_samples
            else:
                env = 1.0

            sample = int(32767 * 0.45 * env * wave_val)
            frames.extend(struct.pack("<h", sample))

    with wave.open(filename, "wb") as wav_file:
        wav_file.setnchannels(1)       # Mono
        wav_file.setsampwidth(2)       # 16-bit PCM (2 bytes per sample)
        wav_file.setframerate(sample_rate)
        wav_file.writeframes(frames)


def main():
    # Print the explicit 64-codon mapping
    print("=== CODON TO NOTE MAPPING (64 Codons Sorted -> 3-Octave C Major Scale C3-B5) ===")
    for i in range(0, 64, 4):
        row = [
            f"{c} ({CODON_TO_AA[c]}) -> {CODON_TO_NOTE[c][0]:<3} ({CODON_TO_NOTE[c][1]:>6.2f} Hz)"
            for c in ALL_CODONS[i:i + 4]
        ]
        print("  |  ".join(row))
    print()

    # Create mutant sequence: 1-based position 20 (0-based index 19), A -> T
    pos_1based = 20
    idx_0based = pos_1based - 1
    assert WT_SEQ[idx_0based] == "A", (
        f"Expected 'A' at position {pos_1based}, found '{WT_SEQ[idx_0based]}'"
    )
    mut_seq = WT_SEQ[:idx_0based] + "T" + WT_SEQ[idx_0based + 1:]

    # Translate both sequences into codons, amino acids, and notes
    wt_codons = split_codons(WT_SEQ)
    mut_codons = split_codons(mut_seq)

    wt_aas = [CODON_TO_AA[c] for c in wt_codons]
    mut_aas = [CODON_TO_AA[c] for c in mut_codons]

    wt_notes = [CODON_TO_NOTE[c] for c in wt_codons]
    mut_notes = [CODON_TO_NOTE[c] for c in mut_codons]

    # Save both WAV files
    wt_wav = "hbb_wildtype.wav"
    mut_wav = "hbb_mutant_A20T.wav"
    write_melody_wav(wt_wav, wt_notes)
    write_melody_wav(mut_wav, mut_notes)
    print(f"Saved wild-type melody to: {wt_wav}")
    print(f"Saved mutant melody to:    {mut_wav}\n")

    # Compare the two versions
    codons_diff = sum(1 for a, b in zip(wt_codons, mut_codons) if a != b)
    aas_diff = sum(1 for a, b in zip(wt_aas, mut_aas) if a != b)
    notes_diff = sum(1 for a, b in zip(wt_notes, mut_notes) if a != b)

    print("=== COMPARISON SUMMARY (Wild-Type vs. A20T Mutant) ===")
    print(f"Codon 7 change:       {wt_codons[6]} -> {mut_codons[6]}")
    print(f"Amino acid 7 change:  {wt_aas[6]} (Glu) -> {mut_aas[6]} (Val)")
    print(f"Note 7 change:        {wt_notes[6][0]} ({wt_notes[6][1]} Hz) -> {mut_notes[6][0]} ({mut_notes[6][1]} Hz)")
    print("-" * 54)
    print(f"Codons that differ:      {codons_diff}")
    print(f"Amino acids that differ: {aas_diff}")
    print(f"Notes that differ:       {notes_diff}")


if __name__ == "__main__":
    main()
