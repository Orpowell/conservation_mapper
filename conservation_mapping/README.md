# Conservation Mapper

Conservation Mapper is PyMol Plug-in for the rudimentary mapping of amino acid variation onto a protein structure in PyMol.
The general method for doing this is as follows:
1. Generate a multiple sequence alignment.
2. Load a protein structure of a sequence in the MSA.
3. Use this script to map amino acid variability onto the protein structure from the MSA.

Residues are coloured 4 colours from white to increasingly darker shades of red, where white residues 
have no variability at the position (I.E total conservation) in the alignment and darkest red 
residues have high variability (3+ amino acids in the alignment). Note: deletions, represented by "-"
are ignored by the tool and don't count towards variability scores.


## Dependencies

- Python3 (included with PyMol)
- Pymol (only tested with Open source version)
- Biotite (Automatically installed by the script)

## Usage

Once you have selected to run the script in PyMol the following command should become available:

	conserve selection, reference, alignment

|parameters|description|
|---|---|
|selection | Pymol object (protein structure). |
|reference | Name of sequence to use as reference from the alignment. This must match to the alignment. The sequence must be identical to the sequence of the selected pymol object|
|alignment | Path to the multiple sequence alignment in a FASTA format.|

## Installation and Usage Notes

Importing and running python libraries in PyMol can be tricky. In our experience, conservation mapper did not work when using pre-built binaries of PyMol (Open Source version).
However, if you build PyMol from source yourself the script works perfectly! I would advise only attempting to build PyMol on a linux system such as Ubuntu.
The instructions for building PyMol from source can be found [here](https://github.com/schrodinger/pymol-open-source)

