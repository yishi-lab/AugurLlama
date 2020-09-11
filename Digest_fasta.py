from enum import Enum
import sys
from Bio import SeqIO
import argparse

class Enzyme(Enum):
    TRYPSIN=1
    CHYMO=2
    GLUC=3

def digest(name:str,sequence:str,miscleavage_allowed:int,min_peptide_length_allowed:int,enzyme):
    peptides=set()

    #find all cleavage position
    cleavage_pos=[0]
    for pos in range(len(sequence)-1):
        aa=sequence[pos]
        aa_=sequence[pos+1]
        if enzyme == Enzyme.CHYMO:
            if (aa == 'F' or aa == 'Y' or aa == 'W' or aa == 'L'):
                cleavage_pos.append(pos+1)
        elif enzyme == Enzyme.TRYPSIN:
            if (aa == 'K' or aa == 'R') and aa_ != 'P':
                cleavage_pos.append(pos+1)
        elif enzyme == Enzyme.GLUC:
            if( aa == 'D' or aa == 'E') and aa_ != 'P':
                cleavage_pos.append(pos+1)

    cleavage_pos.append(len(sequence))

    #get cleaved peptides
    for miscleavage in range(miscleavage_allowed+1):
        for i in range(len(cleavage_pos)-1-miscleavage):
            _start = cleavage_pos[i]
            _end = cleavage_pos[i+1+miscleavage]
            if _end - _start >= min_peptide_length_allowed:
                _pep = sequence[_start:_end]
                peptides.add(_pep)

    #assemble a peptides->protein index for protein assembling in the future
    peptides_protein_index=name+"\n"+"#".join(peptides)+"\n"


    return [peptides,peptides_protein_index]


def process_file(file,enzyme,miscleavage,minlength=5):
    peptides_count={}
    digested_file = open(enzyme.name+'_predigested.fasta','w+')
    peptides_protein_file = open(enzyme.name+'_peptides_protein.fasta','w+')
    for record in SeqIO.parse(file, 'fasta'):
        name = '>'+str(record.id)
        sequence = str(record.seq)
        peptides,peptides_protein_index = digest(name,sequence,miscleavage,minlength,enzyme)
        peptides_protein_file.write(peptides_protein_index)
        for peptide in peptides:
            if peptide in peptides_count:
                peptides_count[peptide] +=1
            else:
                peptides_count[peptide] = 1
    l=1
    for peptide,count in peptides_count.items():
        digested_file.write(">"+str(l)+"-"+str(count)+"\n")
        digested_file.write(peptide+"\n")
        l+=1

    digested_file.close()
    peptides_protein_file.close()


def main():
    parser = argparse.ArgumentParser(description="Description of the protein sequence digestion")
    parser.add_argument('-file',help="fasta file",required=True)
    parser.add_argument('-enzyme',type=int, choices=[1,2,3], required=True,help="enzyme type: [1\tTrypsin] [2\tChymo] [3\tGluC] ")
    parser.add_argument('-miscleavage', type=int, choices=[0, 1, 2, 3],required=True, help="Allowed miscleavage")
    parser.add_argument('-L',type=int,help='minimum peptide length allowed')

    args = parser.parse_args()


    infile = args.file
    type = args.enzyme
    miscleavage = args.miscleavage
    minlength = args.L

    enzyme = Enzyme(type)

    process_file(infile,enzyme,miscleavage,minlength)



main()