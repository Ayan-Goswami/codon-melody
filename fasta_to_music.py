#!/usr/bin/env python3
"""
fasta_to_music.py
Reads a DNA/RNA FASTA file from the terminal and converts each sequence
into a playable WAV melody by mapping codons onto a 3-octave C Major scale.
"""

import argparse
import itertools
import math
import os
import re
import struct
import sys
import wave

# 1. Build a 3-octave C Major scale (21 notes: C3 to B5)
NOTE_NAMES = ["C", "D", "E", "F", "G", "A", "B"]
SEMITONE_OFFSETS = [0, 2, 4, 5, 7, 9, 11]

SCALE_NOTES = []
for octave in (3, 4, 5):
    for name, offset in zip(NOTE_NAMES, SEMITONE_OFFSETS):
        midi_num = 12 * (octave + 1) + offset  # C4 = MIDI 60, A4 = MIDI 69 (440 Hz)
        freq = 440.0 * (2.0 ** ((midi_num - 69) / 12.0))
        SCALE_NOTES.append((f"{name}{octave}", round(freq, 2)))

# 2. Sort all 64 codons alphabetically (AAA .. TTT) and map onto the 21-note scale
BASES = ["A", "C", "G", "T"]
ALL_CODONS = ["".join(c) for c in itertools.product(BASES, repeat=3)]
CODON_TO_NOTE = {
    codon: SCALE_NOTES[i % len(SCALE_NOTES)]
    for i, codon in enumerate(ALL_CODONS)
}


def parse_fasta(filepath):
    """
    Parses a FASTA file and yields (header_id, cleaned_sequence) tuples.
    Supports single or multi-record FASTA files, lowercase letters, and RNA (U -> T).
    """
    if not os.path.isfile(filepath):
        print(f"Error: File '{filepath}' not found.", file=sys.stderr)
        sys.exit(1)

    records = []
    current_header = "sequence_1"
    current_seq = []

    with open(filepath, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith(";"):
                continue
            if line.startswith(">"):
                if current_seq:
                    records.append((current_header, "".join(current_seq)))
                    current_seq = []
                # Take the first word after '>' as a clean record ID
                raw_header = line[1:].strip()
                current_header = raw_header.split()[0] if raw_header else f"seq_{len(records)+1}"
            else:
                # Keep alphabetic characters, uppercase them, and convert RNA U -> T
                cleaned = re.sub(r"[^A-Za-z]", "", line).upper().replace("U", "T")
                current_seq.append(cleaned)

    if current_seq:
        records.append((current_header, "".join(current_seq)))

    return records


def sequence_to_notes(seq, header):
    """Splits a sequence into 3-base codons and maps each codon to (note_name, freq_hz)."""
    remainder = len(seq) % 3
    if remainder != 0:
        print(
            f"  [Note] '{header}' length ({len(seq)} bp) is not a multiple of 3. "
            f"Trimming last {remainder} base(s)."
        )
        seq = seq[:-remainder]

    notes = []
    ambiguous_count = 0
    for i in range(0, len(seq), 3):
        codon = seq[i:i + 3]
        if codon in CODON_TO_NOTE:
            notes.append(CODON_TO_NOTE[codon])
        else:
            # Play a silent rest for ambiguous codons (e.g., containing 'N')
            notes.append(("REST", 0.0))
            ambiguous_count += 1

    if ambiguous_count > 0:
        print(f"  [Note] '{header}' contained {ambiguous_count} ambiguous codon(s) mapped to silence.")

    return notes


def write_melody_wav(filename, notes, note_duration=0.18, sample_rate=44100):
    """Synthesizes a list of (note_name, freq_hz) tuples into a 16-bit mono WAV file."""
    n_samples_per_note = int(sample_rate * note_duration)
    fade_samples = min(int(sample_rate * 0.015), n_samples_per_note // 2)
    frames = bytearray()

    for _, freq in notes:
        if freq <= 0.0:
            # Silent rest
            frames.extend(b"\x00\x00" * n_samples_per_note)
            continue

        for s in range(n_samples_per_note):
            t = s / sample_rate
            wave_val = (
                0.8 * math.sin(2.0 * math.pi * freq * t)
                + 0.2 * math.sin(4.0 * math.pi * freq * t)
            )

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
        wav_file.setsampwidth(2)       # 16-bit PCM
        wav_file.setframerate(sample_rate)
        wav_file.writeframes(frames)


def sanitize_filename(name):
    """Removes characters that are unsafe in filenames."""
    return re.sub(r"[^A-Za-z0-9_.-]", "_", name)


def print_mapping():
    """Prints the 64-codon to musical note mapping table."""
    print("=== CODON TO NOTE MAPPING (64 Codons Sorted -> 3-Octave C Major Scale C3-B5) ===")
    for i in range(0, 64, 4):
        row = [
            f"{c} -> {CODON_TO_NOTE[c][0]:<3} ({CODON_TO_NOTE[c][1]:>6.2f} Hz)"
            for c in ALL_CODONS[i:i + 4]
        ]
        print("  |  ".join(row))
    print()


def main():
    parser = argparse.ArgumentParser(
        description="Convert DNA/RNA sequences in a FASTA file into WAV melodies."
    )
    parser.add_argument(
        "fasta_file",
        help="Path to the input FASTA file (e.g., my_sequence.fasta)",
    )
    parser.add_argument(
        "-o", "--output",
        help="Custom output WAV filename (default: <fasta_filename>.wav)",
        default=None,
    )
    parser.add_argument(
        "-d", "--duration",
        type=float,
        default=0.18,
        help="Duration of each codon note in seconds (default: 0.18)",
    )
    parser.add_argument(
        "--show-mapping",
        action="store_true",
        help="Print the 64-codon to note mapping table before generating audio",
    )

    args = parser.parse_args()

    if args.show_mapping:
        print_mapping()

    records = parse_fasta(args.fasta_file)
    if not records:
        print(f"Error: No valid FASTA sequences found in '{args.fasta_file}'.", file=sys.stderr)
        sys.exit(1)

    base_stem = os.path.splitext(os.path.basename(args.fasta_file))[0]

    print(f"Loaded {len(records)} sequence(s) from '{args.fasta_file}':")

    for idx, (header, seq) in enumerate(records, start=1):
        notes = sequence_to_notes(seq, header)
        if not notes:
            print(f"  [Skipping] '{header}' has fewer than 3 bases.")
            continue

        # Determine output WAV filename
        if len(records) == 1:
            out_wav = args.output if args.output else f"{base_stem}.wav"
        else:
            safe_id = sanitize_filename(header)
            if args.output:
                out_stem, out_ext = os.path.splitext(args.output)
                out_wav = f"{out_stem}_{safe_id}{out_ext or '.wav'}"
            else:
                out_wav = f"{base_stem}_{safe_id}.wav"

        write_melody_wav(out_wav, notes, note_duration=args.duration)
        total_time = len(notes) * args.duration
        print(
            f"  [{idx}/{len(records)}] '{header}': {len(seq)} bp ({len(notes)} codons, "
            f"{total_time:.1f}s) -> Saved to '{out_wav}'"
        )


if __name__ == "__main__":
    main()
