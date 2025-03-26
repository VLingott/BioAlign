import json
import hashlib
import subprocess
import os.path
import sys

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from docx import Document
from docx.shared import Pt, Inches
from docx.enum.text import WD_COLOR_INDEX
from docx.oxml.ns import qn


def prepare_seq(seqs: dict, output_file_name: str):
    if "." in output_file_name:
        input_file_name = "".join(output_file_name.split(".")[:-1]) + ".fasta"
    else:
        input_file_name = output_file_name + ".fasta"
        output_file_name = output_file_name + ".aln"
    with open(input_file_name, "w") as f:
        for key, value in seqs.items():
            SeqIO.write(SeqRecord(Seq(value.replace(" ", "").replace("\n", "").strip()), id=key), f, "fasta")
    command_return = subprocess.run(f".\\clustal-omega-1.2.2-win64\\clustalo.exe --infile {input_file_name} --outfile {output_file_name} --outfmt clustal --force")
    if command_return.returncode != 0:
        print(f"An error occurred while running clustal-omega.\nCommand: {command_return.args}\nReturncode: {command_return.returncode}")
        sys.exit(-1)


def prepare_formatted_seq(aln_file_name: str) -> str:
    with open(aln_file_name, "r") as f:
        aln_lines = f.readlines()

    aln_seq = "".join(aln_lines[:3])

    i = 0
    c = ''
    while True:
        c = aln_lines[3][i]
        if c == ' ':
            break
        i += 1
    while True:
        c = aln_lines[3][i]
        if c != ' ':
            break
        i += 1

    for line in aln_lines[3:]:
        if line.strip() == "":
            aln_seq += line
            continue
        ID = line[:i]
        seq = line[i:]
        seq_n = ""
        for k, c in enumerate(seq):
            seq_n += c
            if (k + 1) % 3 == 0:
                seq_n += " "
        aln_seq += ID + seq_n.rstrip() + "\n"
    return aln_seq


def create_word_document_and_mark(
        filename,
        raw_text,
        font_name="Courier New",
        font_size=9,
        margin_inches=0.5,
        marks=[],
    ):
    document = Document()

    # Set narrow margins
    sections = document.sections
    for section in sections:
        section.top_margin = Inches(margin_inches)
        section.bottom_margin = Inches(margin_inches)
        section.left_margin = Inches(margin_inches)
        section.right_margin = Inches(margin_inches)

    # Add a paragraph style for monospace font
    style = document.styles.add_style("Monospace", 1)  # 1 is for paragraph style
    style.font.name = font_name
    style.font.size = Pt(font_size)
    style._element.rPr[0].set(qn("w:ascii"), font_name)

    paragraph = document.add_paragraph(raw_text, style="Monospace")
    run = paragraph.runs[0]

    # Clear existing runs and apply marks
    if len(marks) != 0:
        paragraph._element.clear()
        i = 0
        while i + 1 <= len(marks):
            t = type(marks[i][0])
            if t != int:
                print(t)

            # Add the non-marked part before the current mark
            if i == 0:
                unmarked_text = raw_text[:marks[i][0]]
            else:
                unmarked_text = raw_text[marks[i - 1][1]:marks[i][0]]
            run = paragraph.add_run(unmarked_text)
            run.font.name = font_name
            run.font.size = Pt(font_size)

            # Add the marked part
            run = paragraph.add_run(raw_text[marks[i][0]:marks[i][1]])
            run.font.name = font_name
            run.font.size = Pt(font_size)
            run.font.highlight_color = marks[i][2]

            i += 1

        if marks[-1][1] < len(raw_text):
            # Add unmarked part at the end
            run = paragraph.add_run(raw_text[marks[-1][1]:])
            run.font.name = font_name
            run.font.size = Pt(font_size)

    document.save(filename)


def mark_sequences(sequence: str, search: str, spaced: bool, separate_colors: bool) -> [tuple[int, int]]:
    def mark_sequence(s: int):
        def skip_lines(lines: int):
            nonlocal sequence_i
            for _ in range(lines):
                if sequence_i + 1 >= len(sequence):
                    return False
                while sequence[sequence_i] != '\n':
                    sequence_i += 1
                    if sequence_i >= len(sequence):
                        return False
                sequence_i += 1
            return True

        sequence_i = 0
        search_i = 0
        start_i = -1
        color = ([WD_COLOR_INDEX.BRIGHT_GREEN, WD_COLOR_INDEX.YELLOW, WD_COLOR_INDEX.PINK])[s % num_sequences] if separate_colors else WD_COLOR_INDEX.YELLOW

        # skip header
        skip_lines(3 + s)

        marked_over_line = False

        while sequence_i < len(sequence):
            # skip label
            while sequence[sequence_i] != ' ':
                sequence_i += 1

            # skip whitespace
            while sequence[sequence_i] == ' ':
                sequence_i += 1

            # find next match
            while sequence[sequence_i] != '\n':
                if sequence[sequence_i] == ' ':
                    if not spaced:
                        start_i = -1
                        search_i = 0
                elif sequence[sequence_i] == search[search_i]:  # current sequence chr == current search term chr
                    if start_i == -1:
                        start_i = sequence_i
                    if search_i + 1 == len(search):  # append search if the whole search term is found
                        matches.append((start_i, sequence_i + 1, color))
                        marked_over_line = False
                        start_i = -1
                        search_i = 0
                    else:
                        search_i += 1
                else:
                    if marked_over_line:
                        matches.pop()
                        marked_over_line = False
                    start_i = -1
                    search_i = 0

                sequence_i += 1

                if sequence[sequence_i] == '\n' and start_i != -1:
                    if spaced:
                        matches.append((start_i, sequence_i + 1, color))
                        marked_over_line = True
                    else:
                        search_i = 0
                    start_i = -1

            skip_lines(num_sequences + 2)

            if marked_over_line and sequence_i + 1 >= len(sequence):
                matches.pop()

    def num_sequences():
        sequence_i = 0

        # skip header
        for _ in range(3):
            while sequence[sequence_i] != '\n':
                sequence_i += 1
            sequence_i += 1

        # count sequences
        sequences = 0
        while sequence[sequence_i] != ' ':
            while sequence[sequence_i] != '\n':
                sequence_i += 1
            sequence_i += 1
            sequences += 1

        return sequences

    matches: [tuple[int, int]] = []
    num_sequences = num_sequences()
    for s in range(num_sequences):
        mark_sequence(s)

    return sorted(matches, key=lambda x: x[0])


if __name__ == "__main__":
    # change sequence in sequences.json file. create if it does not exists.
    # in json file: change this if you want to have other DNA sequences. add as many as you like in the following syntax: {"Name for sequence": "sequence", "Name for sequence": "sequence", ...}
    sequences = json.load(open("sequences.json"))

    # skipping generation of DNA alignment if sequence is unchanged
    path_hash = ".sequencehash"
    last_hash = None
    if os.path.exists(path_hash) and os.path.isfile(path_hash):
        with open(path_hash, "r", encoding="utf-8") as file:
            last_hash = file.read()
    new_hash = hashlib.md5(str(sequences).encode("utf-8"), usedforsecurity=False).hexdigest()
    if last_hash != new_hash or not os.path.exists("sequences.aln") or not os.path.isfile("sequences.aln"):
        with open(path_hash, "w", encoding="utf-8") as file:
            file.write(new_hash)
        print("Computing DNA sequence alignment...", end='')
        prepare_seq(sequences, "sequences.aln")
        print("\rComputed DNA sequence alignment.  \n")
    else:
        print("Reusing unchanged DNA alignment file.\n")

    # format DNA sequence to use triplet notation
    text = prepare_formatted_seq("sequences.aln")

    # get user settings
    word_filename = "sequences.docx"
    matches = []
    search_word = input("Input DNA sequence to search for (e.g. \"ACC\" or \"\" to disable marking): ").strip()
    if len(search_word) != 0:
        skip_spaces = True if "space" in input("Search mode (exact/spaced): ") else False
        separate_marking_colors = True if "y" in input("Separate marking colors for each sequence (yes/no): ") else False

        # compute markings for DNA sequence
        matches = mark_sequences(text, search_word, skip_spaces, separate_marking_colors)

    matches_num = len(matches)
    print(f"{matches_num if matches_num != 0 else "No"} matches of \"{search_word}\" found.")

    # create word document with markings
    create_word_document_and_mark(word_filename, text, marks=matches)

    print(f"\nWord document \"{word_filename}\" created successfully.")
