# Codon Alignment Suite

`CodonAlignmentSuite.py` is a versatile script for aligning nucleotide sequences based on protein alignments. It includes options for translating nucleotide sequences to protein sequences using different codon tables, aligning protein sequences with various aligners, and generating phylogenetic trees in Newick format.

## Features

-   **Protein Sequence Alignment**: Supports alignment using ClustalOmega, MAFFT, and MUSCLE.
-   **Nucleotide Sequence Translation**: Translates nucleotide sequences into protein sequences using specified codon tables.
-   **Phylogenetic Tree Generation**: Generates Newick trees based on codon alignments.
-   **Flexible Output Formats**: Outputs aligned sequences in various formats including CLUSTAL, PAML, FASTA, CODON, and PHYLIP.

## Prerequisites

Before using this script, you need to have the following software installed in your path:

-   **Python 3.6 or higher**
-   **Biopython**: Install using `pip install biopython`
-   **ClustalOmega**: Install from ClustalOmega official site
-   **MAFFT**: Install from MAFFT official site
-   **MUSCLE**: Install from MUSCLE official site

## Usage

Here are examples of how to use the script in different scenarios:

### Example 1: User has only nucleotide file
python CodonAlignmentSuite.py zur.ffn -translate -aligner mafft -codontable 11 -output clustal -outputfile output.aln -tree
### Example 2: User has both nucleotide and protein alignment files
python CodonAlignmentSuite.py prot_aln.clustal nuc1.fasta nuc2.fasta -output clustal -outputfile output.aln
### Example 3: User has nucleotide file and protein file (not aligned)
python CodonAlignmentSuite.py prot_seq.fasta nuc1.fasta nuc2.fasta -align_prot -aligner mafft -output clustal -outputfile output.aln

To generate a Newick tree using codon alignment, add the `-tree` option to any of the above commands.

## Command-line Options

-   `pep_aln` (optional): Protein alignment file in CLUSTAL or FASTA format
-   `nuc_fasta`: DNA sequences in FASTA format (single or multiple files)
-   `-output`: Output format (clustal, paml, fasta, codon, phylip)
-   `-outputfile`: Output file name
-   `-nogap`: Remove gaps from the output
-   `-html`: Output in HTML format
-   `-nostderr`: Suppress standard error messages
-   `-blockonly`: Output only blocks of aligned sequences
-   `-nomismatch`: Exclude mismatches
-   `-codontable`: Codon table to use
    -   1: Universal code (default)
    -   2: Vertebrate mitochondrial code
    -   3: Yeast mitochondrial code
    -   4: Mold, Protozoan, and Coelenterate Mitochondrial code and Mycoplasma/Spiroplasma code
    -   5: Invertebrate mitochondrial
    -   6: Ciliate, Dasycladacean and Hexamita nuclear code
    -   9: Echinoderm and Flatworm mitochondrial code
    -   10: Euplotid nuclear code
    -   11: Bacterial, archaeal and plant plastid code
    -   12: Alternative yeast nuclear code
    -   13: Ascidian mitochondrial code
    -   14: Alternative flatworm mitochondrial code
    -   15: Blepharisma nuclear code
    -   16: Chlorophycean mitochondrial code
    -   21: Trematode mitochondrial code
    -   22: Scenedesmus obliquus mitochondrial code
    -   23: Thraustochytrium mitochondrial code
-   `-tree`: Generate Newick tree using codon alignment
-   `-align_prot`: Align protein sequences first if not already aligned
-   `-aligner`: Aligner to use if protein sequences need to be aligned (clustalo, mafft, muscle)
-   `-translate`: Translate nucleotide sequences to protein if protein file is not provided

## Example Input and Output
### Command Used: python CodonAlignmentSuite.py zur.ffn -translate -aligner mafft -codontable 11 -output clustal -outputfile output.aln -tree
### Input File: `zur.ffn`
### Output Files
-   **Aligned Sequences**: `output.aln`
-   **Newick Tree**: `output.nwk`

## Citation

If you are using the `CodonAlignmentSuite.py` script, please cite it as follows: Sharma, V. (2023). CodonAlignmentSuite.py [Python script]. Retrieved from [https://github.com/vsmicrogenomics/codon-alignment-suite](https://github.com/vsmicrogenomics/codon-alignment-suite)

## Acknowledgements

This script relies on several excellent tools and libraries:

-   Biopython: [Biopython](https://biopython.org/) is an open-source collection of Python tools for biological computation.
-   ClustalOmega: [ClustalOmega](http://www.clustal.org/omega/) is a general-purpose multiple sequence alignment program for proteins.
-   MAFFT: [MAFFT](https://mafft.cbrc.jp/alignment/server/index.html) is a multiple sequence alignment program for unix-like operating systems.
-   MUSCLE: [MUSCLE](http://www.drive5.com/muscle/) is a multiple sequence alignment program with reduced time and space complexity.

Please make sure to cite and acknowledge these tools if you use them in your research.
