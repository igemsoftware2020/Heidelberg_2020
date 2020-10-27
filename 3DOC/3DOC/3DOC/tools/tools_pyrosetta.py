import pyrosetta
import os

def check_pumby_ppr_length(sequence_list):
    '''
    Function checks for Pumby modules published in Adamala K.P., Martin-Alarcon D.A., Boyden E.S. (2016). Programmable RNA-binding protein composed of repeats of a single modular unit. Proc Natl Acad Sci U S A. doi:10.1073/pnas.1519368113 and PPR modules published in Coquille S., Filipovska A., Chia T., Rajappa L., Lingford J.P., Razif M.F.M., Thore S., Rackham O. (2014). An artificial PPR scaffold for programmable RNA recognition. Nat Commun. doi:10.1038/ncomms6729, whether the modules are ordered correctly and have the right amount as specified in the articles.
    
    :param sequence_list: list of sequences inputted by the user
    '''
    pumby_positions = []
    ppr_positions = []
    pumby_start = 0
    pumby_end = 0
    pumby_end_bool = False
    pumby_start_bool = False

    for sequence_choice in range(0,len(sequence_list)):
        sequence = sequence_list[sequence_choice]
        shift_pumby_start = 0
        if sequence[0:2] == "RG":
            pumby_start = sequence_choice
            pumby_start_bool = True
            shift_pumby_start = 2
        if sequence[0+shift_pumby_start:15+shift_pumby_start] == "ELHQHTEQLVQDQYG" and sequence[16+shift_pumby_start:19+shift_pumby_start] == "YVI" and sequence[20+shift_pumby_start:36+shift_pumby_start] == "HVLEHGRPEDKSKIVA":
            pumby_positions.append(sequence_choice)
            if sequence[36:38] == "GR":
                pumby_end = sequence_choice
                pumby_end_bool = True
        if sequence[:3] == "VTY" and sequence[4:33] == "TLISGLGKAGRLEEALELFEEMKEKGIVP" and sequence[34:35] == "V":
            ppr_positions.append(sequence_choice)

    if len(pumby_positions) > 1:
        pumby_connected = True
        for element in range(0, len(pumby_positions)-1):
            if pumby_positions[element + 1] != pumby_positions[element] + 1:
                pumby_connected = False

        if pumby_positions[0] != pumby_start or pumby_start_bool == False:
            print("You did not position the Pumby start module at the beginning of your Pumby sequences. Please check [Adamala K.P., Martin-Alarcon D.A., Boyden E.S. (2016). Programmable RNA-binding protein composed of repeats of a single modular unit. Proc Natl Acad Sci U S A. doi:10.1073/pnas.1519368113].")
        if pumby_positions[-1] != pumby_end or pumby_end_bool == False:
            print("You did not position the Pumby end module at the end of your Pumby sequences. Please check [Adamala K.P., Martin-Alarcon D.A., Boyden E.S. (2016). Programmable RNA-binding protein composed of repeats of a single modular unit. Proc Natl Acad Sci U S A. doi:10.1073/pnas.1519368113].")
        if not pumby_connected:
            print("You did not position the Pumby modules direct neighborhood. Please check [Adamala K.P., Martin-Alarcon D.A., Boyden E.S. (2016). Programmable RNA-binding protein composed of repeats of a single modular unit. Proc Natl Acad Sci U S A. doi:10.1073/pnas.1519368113].")
        if len(pumby_positions) != 6 and len(pumby_positions) != 10 and len(pumby_positions) != 12 and len(pumby_positions) != 18 and len(pumby_positions) != 24:
            print("You have chosen a deviating number of Pumby modules. Stable RNA linkers may only be produced with 6/ 10/ 12/ 18 or 24 modules in total. Please check [Adamala K.P., Martin-Alarcon D.A., Boyden E.S. (2016). Programmable RNA-binding protein composed of repeats of a single modular unit. Proc Natl Acad Sci U S A. doi:10.1073/pnas.1519368113]")

    if len(ppr_positions) > 1:
        if len(ppr_positions) != 8:
            print("You have chosen a deviating number of PPR modules. Stable RNA linkers may only be produced with 8 modules in total. Please check [Coquille S., Filipovska A., Chia T., Rajappa L., Lingford J.P., Razif M.F.M., Thore S., Rackham O. (2014). An artificial PPR scaffold for programmable RNA recognition. Nat Commun. doi:10.1038/ncomms6729]")
        if not pumby_connected:
            print("You did not position the PPR modules direct neighborhood. Please check [Coquille S., Filipovska A., Chia T., Rajappa L., Lingford J.P., Razif M.F.M., Thore S., Rackham O. (2014). An artificial PPR scaffold for programmable RNA recognition. Nat Commun. doi:10.1038/ncomms6729]")


def modification_pose(pose, blast_fragment_start, blast_fragment_end, mutated_amino_acids, additional_amino_acids_begin, additional_amino_acids_end):
    '''
    Function modifies a pose based on inputs mutated_amino_acids, additional_amino_acids_begin, additional_amino_acids_end.
    
    :param pose: pose, which is getting modified
    :param blast_fragment_start: residue of pose, where BLAST fragment starts
    :param blast_fragment_end: residue of pose, where BLAST fragment ends
    :param mutated_amino_acids: amino acids in pose differing from user-inputted sequence
    :param additional_amino_acids_begin: amino acids at the beginning of pose differing from user-inputted sequence
    :param additional_amino_acids_end: amino acids at the end of pose differing from userinputted sequence
    :return pose: pose with mutated amino acids
    '''
    
    # Delete sequence parts before and after the blast fragment

    pyrosetta.rosetta.protocols.grafting.delete_region(pose, blast_fragment_end+1+len(additional_amino_acids_end.keys()), len(pose.sequence()))
    pyrosetta.rosetta.protocols.grafting.delete_region(pose, 1, blast_fragment_start-len(additional_amino_acids_begin.keys()))

    # Mutate the changed amino_acids

    mutated_amino_acids = {**mutated_amino_acids, **additional_amino_acids_begin, **additional_amino_acids_end}

    for aa in mutated_amino_acids.keys():
        pyrosetta.toolbox.mutants.mutate_residue(pose, aa+1, mutated_amino_acids[aa])

    # Repacking for new pose

    scorefxn = pyrosetta.rosetta.protocols.loops.get_fa_scorefxn()

    if mutated_amino_acids > 0:
        task = pyrosetta.standard_packer_task(pose)
        task.restrict_to_repacking()
        task.temporarily_fix_everything()
        for aa in mutated_amino_acids.keys():
            task.temporarily_set_pack_residue(aa, True)

        pack_mover = pyrosetta.rosetta.protocols.minimization_packing.PackRotamersMover(scorefxn, task)
        pack_mover.apply(pose)

    return pose

def check_pose_energy(pose):
    '''
    Function determines energy of a pose via the ref2015 score function.
    
    :param pose: pyrosetta pose
    :return energy: energy of the pose
    '''

    #OPTIONAL
    '''
    pmm = pyrosetta.rosetta.protocols.moves.PyMOLMover()
    pmm.apply(pose)
    pmm.send_energy(pose, label = True)
    '''

    scorefxn = pyrosetta.get_fa_scorefxn()
    energy = scorefxn(pose)

    return energy


def make_rna_pose(sequence):
    '''
    Function converts a RNA sequence in a PyRosetta comaptible format and creates as well as returns a pose for the RNA sequence.

    :param sequence: RNA sequence
    :return pose: pose of RNA sequence
    '''
    
    sequence = sequence.upper()
    sequence_rna = ""
    sequence_base = ""

    for base in sequence:
        if base == "A":
            sequence_base = "A[RAD]"
        if base == "C":
            sequence_base = "C[RCY]"
        if base == "G":
            sequence_base = "G[RGU]"
        if base == "U":
            sequence_base = "U[URA]"

        sequence_rna = sequence_rna + sequence_base

    pose = pyrosetta.io.pose_from_sequence(sequence_rna, "centroid")

    return pose

