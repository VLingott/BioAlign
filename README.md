# BioAlign - DNA Sequence Alignment Tool

BioAlign is a user-friendly tool for DNA sequence alignment and visualization. It uses Clustal Omega for alignment and creates nicely formatted Word documents with customizable sequence highlighting.

## Table of Contents
- [Features](#features)
- [Installation](#installation)
- [Usage](#usage)
  - [Preparing Your Sequences](#preparing-your-sequences)
  - [Running the Application](#running-the-application)
  - [Search Options](#search-options)
  - [Output](#output)
- [Example](#example)
- [Notes](#notes)
- [Developer Setup](#developer-setup)
  - [Installing Python Dependencies](#installing-python-dependencies)
  - [Installing Clustal Omega](#installing-clustal-omega)
  - [Required Files Structure](#required-files-structure)
- [Requirements](#requirements)
- [Troubleshooting](#troubleshooting)
- [Contributing](#contributing)
- [License](#license)
- [Acknowledgements](#acknowledgements)

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

## Developer Setup

If you only want to work with the source code (main.py), you'll need to set up the following dependencies:

### Installing Python Dependencies

```bash
pip install biopython python-docx
```

### Installing Clustal Omega

1. **Windows**:
   - Download Clustal Omega from http://www.clustal.org/omega/
   - Extract the files to a directory named `clustal-omega-1.2.2-win64` in the same folder as main.py
   - Ensure that `clustalo.exe` is directly inside this directory

2. **macOS**:
   - Install via Homebrew: `brew install clustal-omega`
   - Modify the path in the prepare_seq function where Clustal Omega is called to use:
     ```python
     command_return = subprocess.run(f"clustalo --infile {input_file_name} --outfile {output_file_name} --outfmt clustal --force", shell=True)
     ```

3. **Linux**:
   - Install via package manager: `sudo apt install clustal-omega` (Ubuntu/Debian) or equivalent
   - Modify the path in the prepare_seq function where Clustal Omega is called to use:
     ```python
     command_return = subprocess.run(f"clustalo --infile {input_file_name} --outfile {output_file_name} --outfmt clustal --force", shell=True)
     ```

### Required Files Structure

Create the following files in your working directory:
- `main.py` - The main program
- `sequences.json` - Your DNA sequences in JSON format
- `start.bat` (optional) - For easy launching on Windows

## Requirements

For end users, the release package is self-contained and works on Windows systems without additional installations.

For developers working with just the source code:
- Python 3.6 or higher
- Biopython library
- python-docx library
- Clustal Omega (installed as described above)
- Access to write files in the working directory

## Troubleshooting

### Common Issues

1. **"No module named 'Bio'"**
   - Make sure you've installed Biopython: `pip install biopython`

2. **"No module named 'docx'"**
   - Make sure you've installed python-docx: `pip install python-docx`

3. **"An error occurred while running clustal-omega"**
   - Verify that Clustal Omega is installed correctly in the expected location
   - Check that you have executable permissions for the Clustal Omega binary
   - Try running Clustal Omega manually to verify it works outside the application

4. **"ModuleNotFoundError: No module named 'sequences'"**
   - Ensure that `sequences.json` exists in the same directory as main.py
   - Check that the JSON file is properly formatted

5. **Word document has no highlighting**
   - Make sure you entered a valid sequence pattern that exists in your DNA sequences
   - Verify you selected the appropriate search mode for your pattern

## Contributing

Contributions to improve BioAlign are welcome! Here's how you can contribute:

1. Fork the repository
2. Create a feature branch: `git checkout -b new-feature`
3. Make your changes
4. Test your changes thoroughly
5. Submit a pull request with a clear description of the improvements

Please ensure your code maintains compatibility with the existing functionality and follows the project's coding style.

## License

This project is licensed under the GNU General Public License v3.0 - see the LICENSE file for details.

## Acknowledgements

- [Clustal Omega](http://www.clustal.org/omega/) for sequence alignment
- [Biopython](https://biopython.org/) for biological sequence handling
- [python-docx](https://python-docx.readthedocs.io/) for Word document generation
