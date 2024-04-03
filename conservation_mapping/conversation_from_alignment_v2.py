from pymol import cmd
import logging
import os
import sys

# Check if biotite is installed in pymol python3 and install if not found
# N.B this requires open-source pymol built from the github
try:
    import biotite.sequence.io.fasta as fasta

except ModuleNotFoundError:
    import pip

    pip.main(["install", "biotite"])
    import biotite.sequence.io.fasta as fasta

nu_color_scheme = {
            1: "0x10f8a38",  # green for conserved
            2: "0x60f090",  # Light green
            3: "0x8c7fac",  # light purple
            4: "0x692D94",  # purple
        }

nu_nu_color_scheme = {
            1: "0xffffff",  # white
            2: "0xff8A8A",  # light pink
            3: "0xD9381E",  # red
            4: "0x750000",  # dark red
        }

class ConservationMapper:
    def __init__(self, msa, reference, selection) -> None:
        self.selection = selection
        self.msa = msa
        self.reference = reference
        self.alignment_data = None
        self.conservation = []
        self.color_dict = nu_nu_color_scheme


    def load_alignment_data(self):
        # Check alignment file exists
        if not os.path.isfile(self.msa):
            logging.error("Alignment file not found")
            sys.exit(1)

        self.alignment_data = fasta.FastaFile.read(self.msa)

    def validate_reference(self):
        # Check reference in alignment
        if self.reference not in self.alignment_data.keys():
            logging.error("Reference not found in alignment")
            sys.exit(1)

    def compare_pdb_msa_sequences(self):
        # Exctract reference sequence from MSA
        reference_sequence = (
            self.alignment_data[self.reference].replace("*", "").replace("-", "")
        )

        # Extract pdb sequence
        selection_sequence = "".join(cmd.get_fastastr(self.selection).split("\n")[1:])

        print(selection_sequence)

        # Check sequences are identical
        if selection_sequence != reference_sequence:
            logging.error("PDB and alignment sequences must be identical")
            sys.exit(1)

    def generate_conservation_profile(self):
        ref_positions = [
            n
            for n, aa in enumerate(self.alignment_data[self.reference])
            if aa not in ["-", "*"]
        ]

        sequences_cleaned = [
            "".join([sequence[i] for i in ref_positions])
            for sequence in self.alignment_data.values()
        ]

        reference_position_residues = [
            set(residues) for residues in zip(*sequences_cleaned)
        ]

        gapless_residues = [
            residues.discard("-") or residues
            for residues in reference_position_residues
        ]

        # generate binary conservation array
        # 0 = A single residue at the position (conserved)
        # 1 = Mutliple residues at the positon (non-conserved)
        self.conservation = [len(position) for position in gapless_residues]

    def map_conserved_residues(self):
        # Clone structure and color red (default assumption initially is that all residues are unconserved)
        cp_structure = f"{self.selection}_conservation_profile"
        cmd.create(cp_structure, self.selection)
        # cmd.hide(self.selection)
        cmd.color("grey60", cp_structure)

        # color residues according to colour scheme
        # Residues with 4+ are coloured the same
        for n, score in enumerate(self.conservation):
            if score < 4:
                cmd.color(self.color_dict[score], f"resi {n+1} and {cp_structure}")

            else:
                cmd.color(self.color_dict[4], f"resi {n+1} and {cp_structure}")

    def generate_statistics(self):
        total_residues = len(self.conservation)
        conserved_postions = self.conservation.count(1)
        unconserved_postions = total_residues - conserved_postions

        print(
            f"Numer of conserved residues: {conserved_postions} ({conserved_postions/total_residues} %)"
        )
        print(
            f"Numer of unconserved residues: {unconserved_postions} ({unconserved_postions/total_residues} %)"
        )

    def run(self):
        self.load_alignment_data()
        self.validate_reference()
        self.compare_pdb_msa_sequences()
        self.generate_conservation_profile()
        self.map_conserved_residues()
        self.generate_statistics()
        sys.exit(0)


def map_conservation(selection, ref, alignment):
    mapper = ConservationMapper(alignment, ref, selection=selection)
    mapper.run()


cmd.extend("conserve", map_conservation)
