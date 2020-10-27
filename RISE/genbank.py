from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
 
from Bio.SeqIO import InsdcIO
def save(name,record):
    outputfile = open(f"{name}.gb", 'w')
    SeqIO.write(record, outputfile, 'genbank')
    outputfile.close()

def entry_to_genbank(dfrow,filename):


    sequence_string = Seq(dfrow.sequence.upper())
    
    record = SeqRecord(sequence_string,
                       id=str(dfrow.part_name), # random accession number
                       name=dfrow.part_name,
                       description=dfrow.short_desc)
    
    record.annotations['length'] = len(record.seq)   
    record.annotations['molecule_type'] = 'DNA'
    record.annotations['topology'] = 'linear'
    #record.annotations['data_file_divison'] =
    record.annotations['date'] = dfrow.creation_date 
    record.annotations['organism'] = dfrow.source
    target = None
    for val in dfrow.seq_edit_cache.split(';'):
        if 'seqFeatures' in val:
            target = val
            break
    if target is None:
        save(f"{filename}_{''.join(dfrow.part_name.split('.'))}",record)
        return

    array = target.split('new Array')[1].split('(')[1].split(')')[0].split('],')
    for feat in array:
        
        sp = feat.split(',')
        labels = {'label': sp[3]}
        if int(sp[1]) > int(sp[2]):
            continue
        
        feature = SeqFeature(FeatureLocation(start=int(sp[1]), end=int(sp[2])), type='misc_feature',qualifiers=labels)
        record.features.append(feature)
    save(f"{filename}_{''.join(dfrow.part_name.split('.'))}",record)
    

