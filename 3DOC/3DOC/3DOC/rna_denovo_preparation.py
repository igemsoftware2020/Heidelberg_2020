import RNA
import argparse

def secondary_structure_rna(rna_motif):
    '''
    Function takes a RNA-motif and outputs a Secondary Structure of the RNA-motif based on the RNA fold function in VIennaRNA2.0.
    
    :param rna_motif: RNA-motif
    :return secstruct_rna: secondary structure of RNA
    '''

    secstruct_rna = RNA.fold(rna_motif)[0]
	
    return secstruct_rna


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--rna_motif", help = "RNA-motif bound by protein", type = str, required=True)

    args = parser.parse_args()
    rna_motif = args.rna_motif
	
    secstruct_rna = secondary_structure_rna(rna_motif)
    
    print(secstruct_rna)
