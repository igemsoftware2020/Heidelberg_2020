from Bio.Blast import NCBIWWW
from Bio import SearchIO

def blast_aa_sequence(sequence, program):
    """BLAST with blastp 
    
    The function takes a fasta file as input and performs a qblast (blastp) search against sequences in the database "rcsb". It outputs the  and an blasted sequence.
    
    :param fasta_string: sequence
    :param program: BLAST database
    :return id: id of the entry chosen of the database searched
    :return hit: hit of the BLAST search
    :return mutated_amino_acids: Amino acids in query differing from hit
    :return additional_amino_acids_begin: Amino acids at the end of query sequence, but not found in hit
    :return additional_amino_acids_end: Amino acids at the beginning of query sequence, but not found in hit
    """

    result_handle = NCBIWWW.qblast("blastp", program, sequence)

    with open("blast_output.xml", "w") as out_handle:
      out_handle.write(result_handle.read())
    result_handle.close()

    blastp_result = SearchIO.read("blast_output.xml", "blast-xml")

    if len(blastp_result) > 1:

        blast_hsp = blastp_result[0][0]

        if blast_hsp.bitscore > 10 and blast_hsp.evalue < 0.05:  # ANPASSEN
            
            print("Bitscore")
            print(blast_hsp.bitscore)
            print("E-value")
            print(blast_hsp.evalue)

            blast_frag = blastp_result[0][0][0]# first hit, first hsp, first fragment

            id = str(blast_frag.hit_id)
            id = id.split("|")[1]
            hit = str(blast_frag.hit.seq)
            query = str(blast_frag.query.seq)

            mutated_amino_acids = {}
            additional_amino_acids_begin = {}
            additional_amino_acids_end = {}

            begin = False
            end_counter = 0

            for aa in range(len(hit)-1, -1, -1):
                if query[aa] != hit[aa]:
                    additional_amino_acids_end[aa] = query[aa]
                    end_counter += 1
                else:
                    break

            for aa in range(0, len(hit)):
                if aa < len(hit) - end_counter:
                    if begin == False:
                        if query[aa] != hit[aa]:
                            additional_amino_acids_begin[aa] = query[aa]
                            if query[aa + 1] == hit[aa + 1]:
                                begin = True
                        else:
                            begin = True
                    if begin == True:
                        if query[aa] != hit[aa]:
                            mutated_amino_acids[aa] = query[aa]

            percentage_changed_amino_acids = (len(mutated_amino_acids) + len(additional_amino_acids_begin) + len(additional_amino_acids_end)) / len(hit)
            print("Percentage changed amino acids")
            print(percentage_changed_amino_acids)
            percentage_changed_amino_acids_begin = len(additional_amino_acids_begin) / len(hit)
            percentage_changed_amino_acids_end = len(additional_amino_acids_end) / len(hit)
            if percentage_changed_amino_acids < 0.2 and percentage_changed_amino_acids_begin < 0.1 and percentage_changed_amino_acids_end < 0.1:
                return id, hit, mutated_amino_acids, additional_amino_acids_begin, additional_amino_acids_end
            else:
                return 0, 0, 0, 0, 0
        else:
            return 0, 0, 0, 0, 0
    else:
        return 0, 0, 0, 0, 0

def find_fragment_start_end(pose, hit):
    
    pose_sequence = pose.sequence()
    
    blast_fragment_start = pose_sequence.find(hit)
    blast_fragment_end = blast_fragment_start + len(hit)
    
    return blast_fragment_start, blast_fragment_end

'''
def find_fragment_start_end_old(pose, hit):

    pose_sequence = pose.sequence()

    for aa_pose in range(0, len(pose_sequence)):
        if pose_sequence[aa_pose] == hit[0]:
                blast_fragment_start = aa_pose
                counter_hit = 1
                print("aa_pose")
                print(hit[0])
                print(pose_sequence[blast_fragment_start])
                for aa_hit in range(blast_fragment_start, blast_fragment_start + len(hit)):
                    print("aa_hit")
                    #print(hit[0])
                    #print(pose_sequence[aa_hit])
                    try:
                        if hit[aa_hit-blast_fragment_start] == pose_sequence[aa_hit]:
                            counter_hit += 1
                    except:
                        continue
                if counter_hit == len(hit):
                    break

    blast_fragment_end = blast_fragment_start + len(pose_sequence) - 1

    return blast_fragment_start, blast_fragment_end
'''
