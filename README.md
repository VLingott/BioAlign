# BioAlign - DNA Sequence Alignment Tool

BioAlign is a user-friendly tool for DNA sequence alignment and visualization. It uses Clustal Omega for alignment and creates nicely formatted Word documents with customizable sequence highlighting.

## Features

- Automatic DNA sequence alignment using Clustal Omega
- Triplet notation formatting for better readability
- Search and highlight specific DNA sequences in the alignment
- Option for spaced or exact sequence matching
- Different highlight colors for each sequence (optional)
- Caching of alignment results for unchanged sequences

## Installation

No installation required! The release zip file contains everything you need:

1. Unzip the file to any location on your computer
2. Run `start.bat` to launch the application

The package includes:
- Embedded Python 3.13.2 runtime
- Clustal Omega 1.2.2 executable
- All required Python dependencies

## Usage

### Preparing Your Sequences

Create or edit the `sequences.json` file in the application folder. This file should contain your DNA sequences in the following format:

```json
{
  "Sequence1": "ATGCCTGACCTAGTCGATCGATGCTA",
  "Sequence2": "ATGCGTGACCTAGTTGATCGATGCTA",
  "Sequence3": "ATGCCTGACCAAGTCGATCTATGCTA"
}
```

Where:
- Each key is the sequence name
- Each value is the DNA sequence
- You can add as many sequences as needed

### Running the Application

1. Double-click the `start.bat` file
2. The program will:
   - Load your sequences from sequences.json
   - Perform sequence alignment (or use cached results if unchanged)
   - Prompt you for search options

### Search Options

When prompted:

1. Enter a DNA sequence to search for (e.g., "CTG") or leave empty to disable highlighting
2. Choose search mode:
   - `exact`: Matches only exact sequences without spaces
   - `spaced`: Matches sequences allowing for spaces between nucleotides
3. Choose whether to use separate colors for each sequence (yes/no)

### Output

The program generates:
- `sequences.fasta`: The input file for Clustal Omega
- `sequences.aln`: The alignment result from Clustal Omega
- `sequences.docx`: The final Word document with formatted alignment and highlighting

## Example

For the provided example sequences, searching for "CTG" with spaced mode enabled will highlight this pattern in all sequences, allowing you to easily compare variations.

## Notes

- The tool caches alignment results to avoid redundant calculations
- The Word document uses Courier New font for consistent spacing
- Highlighting uses yellow by default, or green/blue/pink when using separate colors

## Requirements

This package is self-contained and works on Windows systems without additional installations.

## Acknowledgements

- [Clustal Omega](http://www.clustal.org/omega/) for sequence alignment
- [Biopython](https://biopython.org/) for biological sequence handling
- [python-docx](https://python-docx.readthedocs.io/) for Word document generation
