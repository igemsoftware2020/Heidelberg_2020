import argparse

def translate_dna_to_rna(dna):
    '''
    Function translates a DNA into a RNA based on Table 1 in Article: Kunkler, Charlotte N and Hulewicz, Jacob P and Hickman, Sarah C and Wang, Matthew C and McCown, Phillip J and Brown, Jessica A (2019) {Stability of an RNA•DNA–DNA triple helix depends on base triplet composition and length of the RNA third strand}.Nucleic Acids Research 47, 7213-7222.'
   
    :param dna: DNA sequence
    :param rna: RNA sequence
    '''
    
    rna = ""
    
    for nt in dna:
        if nt == "A" or nt == "a":
            rna += "U"
        elif nt == "T" or nt == "t":
            rna += "A"
        elif nt == "C" or nt == "c":
            rna += "U"
        elif nt == "G" or nt == "g":
            rna += "C"
        else:
            print("Error: input other than DNA sequence.")
            rna = ""
            break
            
    return rna
       

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Translate DNA in Triple-Helix RNA based on Table 1 in Article: Kunkler, Charlotte N and Hulewicz, Jacob P and Hickman, Sarah C and Wang, Matthew C and McCown, Phillip J and Brown, Jessica A (2019) {Stability of an RNA•DNA–DNA triple helix depends on base triplet composition and length of the RNA third strand}.Nucleic Acids Research 47, 7213-7222.')
    parser.add_argument('--dna_sequence', type=str, help='DNA sequence to translate in RNA for Triple-Helix Design.')
    
    args = parser.parse_args()
    dna = args.dna_sequence

    #variable initialization
    rna = ""
    error = False

    rna = translate_dna_to_rna(dna)
    
    if rna != "":
        print(rna)
  