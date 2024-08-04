import argparse
import sys
import os
import subprocess
from Bio import AlignIO, Phylo
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from Bio import SeqIO
from Bio.Seq import translate
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor, DistanceCalculator
from Bio.Data import CodonTable

def read_protein_alignment(file_path):
    """Reads protein alignment file and determines its format."""
    with open(file_path, "r") as file:
        first_line = file.readline().strip()
        if first_line.startswith("CLUSTAL"):
            format = "clustal"
        elif first_line.startswith(">"):
            format = "fasta"
        elif first_line.startswith("Gblocks"):
            format = "gblocks"
        else:
            format = "clustal"
    return AlignIO.read(file_path, format)

def align_protein_sequences(input_file, output_file, aligner):
    """Aligns protein sequences using the specified aligner."""
    if aligner == "clustalo":
        cmd = ["clustalo", "-i", input_file, "-o", output_file, "--force", "--auto"]
    elif aligner == "mafft":
        cmd = ["mafft", "--auto", input_file]
    elif aligner == "muscle":
        cmd = ["muscle", "-in", input_file, "-out", output_file]

    if aligner == "mafft":
        with open(output_file, "w") as outfile:
            result = subprocess.run(cmd, stdout=outfile, stderr=subprocess.PIPE)
    else:
        result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    
    if result.returncode != 0:
        raise RuntimeError(f"Error running {aligner}: {result.stderr.decode()}")

def read_nucleotide_sequences(file_paths):
    """Reads nucleotide sequences from multiple FASTA files."""
    nucleotide_sequences = {}
    for file_path in file_paths:
        with open(file_path, "r") as file:
            seq_records = SeqIO.parse(file, "fasta")
            for record in seq_records:
                if record.id in nucleotide_sequences:
                    nucleotide_sequences[record.id].seq += record.seq
                else:
                    nucleotide_sequences[record.id] = record
    return nucleotide_sequences

def translate_nucleotide_to_protein(nucleotide_sequences, codon_table):
    """Translates nucleotide sequences into protein sequences using the specified codon table."""
    protein_sequences = {}
    for id, record in nucleotide_sequences.items():
        protein_seq = record.seq.translate(table=codon_table)
        protein_sequences[id] = SeqRecord(protein_seq, id=id, description="")
    return protein_sequences

def align_nucleotide_sequences(protein_alignment, nucleotide_sequences, codon_table=1):
    """Aligns nucleotide sequences based on the given protein alignment."""
    aligned_nucleotides = []
    for protein_record in protein_alignment:
        nucleotide_seq = nucleotide_sequences[protein_record.id]
        aligned_seq = ""
        nuc_index = 0
        for aa in protein_record.seq:
            if aa == "-":
                aligned_seq += "---"
            else:
                aligned_seq += str(nucleotide_seq.seq[nuc_index:nuc_index+3])
                nuc_index += 3
        aligned_nucleotides.append(SeqRecord(Seq(aligned_seq), id=protein_record.id))
    return MultipleSeqAlignment(aligned_nucleotides)

def remove_gaps(alignment):
    """Removes gaps from the alignment."""
    for record in alignment:
        record.seq = Seq(str(record.seq).replace("-", ""))
    return alignment

def exclude_mismatches(alignment, protein_alignment):
    """Excludes mismatches from the alignment."""
    filtered_records = []
    for nuc_record, prot_record in zip(alignment, protein_alignment):
        if not any(nuc == "-" and prot != "-" for nuc, prot in zip(nuc_record.seq, prot_record.seq)):
            filtered_records.append(nuc_record)
    return MultipleSeqAlignment(filtered_records)

def output_blocks(alignment, block_size=10):
    """Outputs only blocks of aligned sequences without gaps."""
    blocks = []
    for record in alignment:
        block_seq = ""
        for i in range(0, len(record.seq), block_size):
            block = record.seq[i:i+block_size]
            if "---" not in block:
                block_seq += str(block)
        blocks.append(SeqRecord(Seq(block_seq), id=record.id))
    return MultipleSeqAlignment(blocks)

def write_output(alignment, output_file="output.aln", format="clustal", html=False):
    """Writes the alignment output in the specified format."""
    if format == "clustal":
        AlignIO.write(alignment, output_file, "clustal")
    elif format == "fasta":
        SeqIO.write(alignment, output_file, "fasta")
    elif format == "paml":
        with open(output_file, "w") as file:
            for record in alignment:
                file.write(f"{record.id} {record.seq}\n")
    elif format == "codon":
        with open(output_file, "w") as file:
            for record in alignment:
                codon_seq = " ".join([str(record.seq[i:i+3]) for i in range(0, len(record.seq), 3)])
                file.write(f"{record.id} {codon_seq}\n")
    elif format == "phylip":
        AlignIO.write(alignment, output_file, "phylip")
    if html:
        with open("output.html", "w") as html_file:
            html_file.write("<html><body><pre>")
            with open(output_file, "r") as output:
                html_file.write(output.read())
            html_file.write("</pre></body></html>")

def generate_newick_tree(alignment, output_file):
    """Generates a Newick tree from the codon alignment."""
    calculator = DistanceCalculator('identity')
    constructor = DistanceTreeConstructor(calculator, 'nj')
    tree = constructor.build_tree(alignment)
    Phylo.write(tree, output_file, "newick")

def main():
    parser = argparse.ArgumentParser(
        description="Align nucleotide sequences based on protein alignment.\n\n"
                    "Usage examples:\n"
                    "1. User has both nucleotide and protein alignment files:\n"
                    "   python AlNuclfromProtAl.py prot_aln.clustal nuc1.fasta nuc2.fasta -output clustal -outputfile output.aln\n\n"
                    "2. User has nucleotide file and protein file (not aligned):\n"
                    "   python AlNuclfromProtAl.py prot_seq.fasta nuc1.fasta nuc2.fasta -align_prot -aligner mafft -output clustal -outputfile output.aln\n\n"
                    "3. User has only nucleotide file:\n"
                    "   python AlNuclfromProtAl.py nuc1.fasta nuc2.fasta -translate -aligner muscle -codontable 1 -output clustal -outputfile output.aln\n\n"
                    "To generate a Newick tree using codon alignment, add the -tree option to any of the above commands:\n"
                    "   python AlNuclfromProtAl.py ... -tree\n",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument("pep_aln", nargs='?', help="Protein alignment file in CLUSTAL or FASTA format")
    parser.add_argument("nuc_fasta", nargs='+', help="DNA sequences in FASTA format (single or multiple files)")
    parser.add_argument("-output", choices=["clustal", "paml", "fasta", "codon", "phylip"], default="clustal", help="Output format")
    parser.add_argument("-outputfile", default="output.aln", help="Output file name")
    parser.add_argument("-nogap", action="store_true", help="Remove gaps from the output")
    parser.add_argument("-html", action="store_true", help="Output in HTML format")
    parser.add_argument("-nostderr", action="store_true", help="Suppress standard error messages")
    parser.add_argument("-blockonly", action="store_true", help="Output only blocks of aligned sequences")
    parser.add_argument("-nomismatch", action="store_true", help="Exclude mismatches")
    parser.add_argument("-codontable", type=int, default=1, help=(
        "Codon table to use:\n"
        "  1  Universal code (default)\n"
        "  2  Vertebrate mitochondrial code\n"
        "  3  Yeast mitochondrial code\n"
        "  4  Mold, Protozoan, and Coelenterate Mitochondrial code\n"
        "     and Mycoplasma/Spiroplasma code\n"
        "  5  Invertebrate mitochondrial\n"
        "  6  Ciliate, Dasycladacean and Hexamita nuclear code\n"
        "  9  Echinoderm and Flatworm mitochondrial code\n"
        " 10  Euplotid nuclear code\n"
        " 11  Bacterial, archaeal and plant plastid code\n"
        " 12  Alternative yeast nuclear code\n"
        " 13  Ascidian mitochondrial code\n"
        " 14  Alternative flatworm mitochondrial code\n"
        " 15  Blepharisma nuclear code\n"
        " 16  Chlorophycean mitochondrial code\n"
        " 21  Trematode mitochondrial code\n"
        " 22  Scenedesmus obliquus mitochondrial code\n"
        " 23  Thraustochytrium mitochondrial code\n"
    ))
    parser.add_argument("-tree", action="store_true", help="Generate Newick tree using codon alignment")
    parser.add_argument("-align_prot", action="store_true", help="Align protein sequences first if not already aligned")
    parser.add_argument("-aligner", choices=["clustalo", "mafft", "muscle"], default="clustalo", help="Aligner to use if protein sequences need to be aligned")
    parser.add_argument("-translate", action="store_true", help="Translate nucleotide sequences to protein if protein file is not provided")

    args = parser.parse_args()

    if args.nostderr:
        sys.stderr = open(os.devnull, 'w')

    if args.translate and args.pep_aln is None:
        nucleotide_sequences = read_nucleotide_sequences(args.nuc_fasta)
        protein_sequences = translate_nucleotide_to_protein(nucleotide_sequences, args.codontable)
        temp_protein_file = "temp_protein.fasta"
        SeqIO.write(protein_sequences.values(), temp_protein_file, "fasta")
        temp_aligned_prot_file = "temp_aligned_prot.aln"
        align_protein_sequences(temp_protein_file, temp_aligned_prot_file, args.aligner)
        protein_alignment = read_protein_alignment(temp_aligned_prot_file)
        os.remove(temp_protein_file)
    elif args.align_prot:
        temp_aligned_prot_file = "temp_aligned_prot.aln"
        align_protein_sequences(args.pep_aln, temp_aligned_prot_file, args.aligner)
        protein_alignment = read_protein_alignment(temp_aligned_prot_file)
    else:
        protein_alignment = read_protein_alignment(args.pep_aln)

    nucleotide_sequences = read_nucleotide_sequences(args.nuc_fasta)
    aligned_nucleotides = align_nucleotide_sequences(protein_alignment, nucleotide_sequences, args.codontable)
    
    if args.nogap:
        aligned_nucleotides = remove_gaps(aligned_nucleotides)
    
    if args.nomismatch:
        aligned_nucleotides = exclude_mismatches(aligned_nucleotides, protein_alignment)
    
    if args.blockonly:
        aligned_nucleotides = output_blocks(aligned_nucleotides)
    
    write_output(aligned_nucleotides, output_file=args.outputfile, format=args.output, html=args.html)
    
    if args.tree:
        tree_file = args.outputfile.rsplit(".", 1)[0] + ".nwk"
        generate_newick_tree(aligned_nucleotides, tree_file)

    if args.align_prot or args.translate:
        os.remove(temp_aligned_prot_file)

if __name__ == "__main__":
    main()
