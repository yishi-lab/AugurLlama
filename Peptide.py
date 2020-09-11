from enum import Enum
import numpy as np
class Belonging(Enum):
    CDR1=1
    CDR2=2
    CDR3=3
    CONSTANT=4

class PeptideFragments(object):
    def __init__(self,fragments:list):
        self.fragments = fragments
        self.fragment_ions = np.sum(fragments)
    def __len__(self):
        return len(self.fragments)

class Domain(object):
    def __init__(self,region:list):
        self.start = region[0]
        self.end = region[1]
    def __len__(self):
        return self.end- self.start

class Peptide(object):

    def __init__(self,sequence:str,mz:float,charge:int,retention_time:float,PEP:float, qvalue:float,index:int):
        '''
        Initialzer for peptide object
        :param sequence:
        :param mz:
        :param charge:
        :param retention_time:
        :param PEP:
        '''
        self.sequence = sequence
        #self.sequence =
        self.mz = mz
        self.charge = charge
        self.retention_time = retention_time
        self.PEP = PEP
        self.qvalue = qvalue
        self.proteins_ids = set()
        self.index = index
        self.fragments = None
        self.additional = list()

    def __eq__(self, other):
        if self.sequence == other.sequence:
            return True
        return False

    def __hash__(self):
        return self.sequence.__hash__()

    def __repr__(self):
        return self.__str__()

    def __str__(self):
        _str_peptide  = self.sequence +"," + str(round(self.mz,4))+","+str(self.charge) +","+str(round(self.retention_time,2))+","
        _str_peptide += str(self.PEP) +"," + str(self.qvalue) +","+str(len(self.proteins_ids))+"\n"

        for _additional in self.additional:
            _str_peptide += str(_additional)

        return _str_peptide

    def add_protein_id(self,protein_id):
        '''
        Add relevant protein to this peptide
        :param protein_id:
        :return:
        '''
        self.proteins_ids.add(protein_id)

    def annotate_fragment_ions(self,fragments:PeptideFragments):
        '''
        Attach fragments ions information on the peptide
        :param fragments:
        :return:
        '''
        self.fragments = fragments


class PeptideMatching(object):
    def __init__(self,peptide:Peptide,span:list,total_fragment_ions:int):
        self.peptide = peptide
        self.span = Domain(span)
        self.belonging = Belonging.CONSTANT
        self.total_fragment_ions = total_fragment_ions

    def set_belonging(self,belonging:Belonging):
        self.belonging = belonging


    def set_coverage(self,coverage:float):
        self.cdr_coverage = coverage

    def set_fragmentation_coverage(self,coverage:float):
        self.cdr_fragmentation_coverage = coverage

    def set_cdr_fragmentation_ions(self,ions:int):
        self.cdr_fragment_ions = ions

    def __eq__(self, other):
        return self.peptide == other.peptide

    def __hash__(self):
        return self.peptide.__hash__()

    def __repr__(self):
        return self.__str__()
    def __str__(self):
        __str_matching = self.peptide.sequence+"|" +str(round(self.cdr_coverage,3))+"|"
        if self.belonging != Belonging.CONSTANT:
            __str_matching += str(round(self.cdr_fragmentation_coverage,3)) +"|"
        __str_matching += str(self.total_fragment_ions) +"|"
        if self.belonging != Belonging.CONSTANT:
            __str_matching += str(self.cdr_fragment_ions)

        return __str_matching

    def __len__(self):
        return len(self.span)

    def __lt__(self, other):
        if self.span.start == other.span.start:
            return len(self) < len(other)
        return self.span.start < other.span.start