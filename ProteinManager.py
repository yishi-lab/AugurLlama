from Base import *
from Peptide import  Domain,Peptide,PeptideMatching,Belonging,PeptideFragments
from Config import Param
import time

W_RPattern = re.compile(r'W\wR')
Y_CPattern = re.compile(r'(Y\wC)|(YY\w)')
WG_GPattern = re.compile(r'(WG\wG)|(W\wQG)|(\wGQG)|(\w\wQG)|(WG\w\w)|(W\w\wG)|(W\wQ\w)')


class CDRAnnotator(object):

    @staticmethod
    def find_cdr1(sequence):
        # STARTING POS. OF CDR1:
        left_area = sequence[20:26]  # look from pos. 20 - 26 of seq (0-based)
        left_cdr = -1
        la_i = left_area.find('SC')
        if la_i < 0:
            # didn't find 'SC', look for 'C'
            la_i = left_area.find('C')
        else:
            la_i += 1  # 'C' is our marker, so advance past 'S'
        if la_i >= 0:
            left_cdr = la_i + 20 + 5  # CDR1 starts at 'C' + 5 (add 20 to put it back in the full sequence)

        # ENDING POS. OF CDR1:
        right_area = sequence[32:40]  # look from pos. 32 - 40 of seq (0-based)
        ra_i = -1
        right_cdr = -1
        W_R = W_RPattern.search(right_area)
        if W_R != None:
            # if we found 'WXR', find its index
            ra_i = right_area.find(W_R[0])

        else:
            ra_i = right_area.find('W')  # didn't find 'WXR', look for 'W'
        if ra_i >= 0:
            right_cdr = ra_i + 32 - 1 + 1  # CDR1 ends at 'W' - 1 (add 32 to put it back in the full sequence)

        # check if st/end found and if not follow rules:
        if left_cdr == -1 and right_cdr == -1:
            left_cdr = 28
            right_cdr = 36
        elif left_cdr == -1:
            left_cdr = right_cdr - 8
        elif right_cdr == -1:
            right_cdr = left_cdr + 8

        return [left_cdr, right_cdr]

    @staticmethod
    def find_cdr2(sequence):
        # STARTING POS. OF CDR2:
        left_area = sequence[32:40]  # look from pos. 32 - 40 of seq (0-based)
        la_i = -1
        left_cdr = -1
        W_R = W_RPattern.search(left_area)
        if W_R != None:
            # if we found 'WXR', find its index
            la_i = W_R.start(0)
        else:
            la_i = left_area.find('W')  # didn't find 'WXR', look for 'W'
        if la_i >= 0:
            left_cdr = la_i + 32 + 14  # CDR2 starts at 'W' + 14 (add 32 to put it back in the full sequence)

        # ENDING POS. OF CDR2:
        right_area = sequence[63:72]  # look from pos. 63 - 72 of seq (0-based)
        right_cdr = -1
        ra_i = right_area.find('RF')
        if ra_i >= 0:
            right_cdr = ra_i + 63 - 8 + 1  # CDR2 ends at 'R' - 8 (add 63 to put it back in the full sequence)

        # check if st/end found and if not follow rules:
        if left_cdr == -1 and right_cdr == -1:
            left_cdr = 51
            right_cdr = 60
        elif left_cdr == -1:
            left_cdr = right_cdr - 9
        elif right_cdr == -1:
            right_cdr = left_cdr + 9

        return [left_cdr, right_cdr]

    @staticmethod
    def find_cdr3(sequence):
        sequence = sequence[:-39]

        left_area = sequence[90:105]
        la_i = -1
        left_cdr = -1
        Y_C = Y_CPattern.search(left_area)
        if Y_C != None:
            # if we found 'YXR', find its index
            la_i = Y_C.start(0) + 2
        else:
            la_i = left_area.find('C')  # didn't find 'YXC', look for 'C'

        if la_i >= 0:
            left_cdr = la_i + 90 + 3
        n = len(sequence) - 1
        n1 = n - 14
        subtract_amount = 1
        right_area = sequence[n1:n - 4]
        ra_i = -1
        right_cdr = -1
        WG_G = WG_GPattern.search(right_area)
        if WG_G != None:
            ra_i = WG_G.start(0)

        if ra_i >= 0:
            right_cdr = ra_i + n1 - subtract_amount + 1  # CDR3 ends at 'W' - 1 (or 'Q' - 3) (add n-14 to put it back in the full sequence)
        # check
        if left_cdr == -1 and right_cdr == -1:
            left_cdr = n - 21
            right_cdr = n - 10
        elif left_cdr == -1:
            left_cdr = right_cdr - 11
        elif right_cdr == -1:
            if left_cdr + 11 <= n:
                right_cdr = left_cdr + 11
            else:
                right_cdr = n
        if left_cdr > right_cdr:
            left_cdr = n - 1
            right_cdr = n
        return [left_cdr, right_cdr]


class Protein(object):
    def __init__(self, id: str, sequence: str, id_file:str):
        self.id = id
        self.sequence = sequence
        self.source = id_file
        self.peptides_matchings = list()
        self.cdr1_coverage = 0
        self.cdr2_coverage = 0
        self.cdr3_coverage = 0
        self.cdr_coverage = 0
        self.nonessential_peps=[]
        self.coverage = 0
        self.annotate_cdr()

    def annotate_cdr(self):

        self.domain_cdr1 = Domain(CDRAnnotator.find_cdr1(self.sequence))
        self.domain_cdr2 = Domain(CDRAnnotator.find_cdr2(self.sequence))
        self.domain_cdr3 = Domain(CDRAnnotator.find_cdr3(self.sequence))

    def merge(self, other):
        _peptides_matching_set = set(self.peptides_matchings)
        for _peptide_matching in other.peptides_matchings:
            if not _peptide_matching in _peptides_matching_set:
                self.peptides_matchings.append(_peptide_matching)

    def add_peptide(self, peptide: Peptide, peptide_fragmentation: PeptideFragments):
        _added = False
        _peptide_sequence = peptide.sequence
        _span_start = self.sequence.find(_peptide_sequence)
        if _span_start < 0:
            return None
        _span = [_span_start, _span_start + len(_peptide_sequence)]
        _ions = peptide_fragmentation.fragment_ions
        _peptide_matching = PeptideMatching(peptide, _span, _ions)
        if self.determine_peptides_belongings_coverage(_peptide_matching, peptide_fragmentation.fragments):
            # peptide belong to CDR region
            # fullfill the peptide coverage requirement, add to protein
            self.peptides_matchings.append(_peptide_matching)
            _added = True


        elif _peptide_matching.belonging == Belonging.CONSTANT:
            # peptide belongs to constant region only
            #self.peptides_matchings.append(_peptide_matching)
            self.nonessential_peps.append(_span)
            _added = True


        return (_added,_peptide_matching.belonging)

    def determine_peptides_belongings_coverage(self, peptide_matching: PeptideMatching, peptide_fragmentation: list):
        _peptide_span = peptide_matching.span

        '''
        If can be aligned to CDR1
        '''
        overlap_span = find_overlap(self.domain_cdr1, _peptide_span)
        overlap = overlap_span[1] - overlap_span[0]

        if overlap >= 0:
            coverage = overlap / len(self.domain_cdr1)
            if coverage >= Param.CDR1_PEPTIDE_COVERAGE:
                cdr_fragments = find_cdr_fragments(_peptide_span.start, overlap_span, peptide_fragmentation)
                cdr_fragment_coverage = cdr_fragments / len(self.domain_cdr1)
                peptide_matching.set_belonging(Belonging.CDR1)
                peptide_matching.set_coverage(coverage)
                peptide_matching.set_cdr_fragmentation_ions(cdr_fragments)
                peptide_matching.set_fragmentation_coverage(cdr_fragment_coverage)
                return True
            return False

        '''
        If can be aligned to CDR2
        '''
        overlap_span = find_overlap(self.domain_cdr2, _peptide_span)
        overlap = overlap_span[1] - overlap_span[0]

        if overlap > 0:
            coverage = overlap / len(self.domain_cdr2)
            if coverage >= Param.CDR2_PEPTIDE_COVERAGE:
                cdr_fragments = find_cdr_fragments(_peptide_span.start, overlap_span, peptide_fragmentation)
                cdr_fragment_coverage = cdr_fragments / len(self.domain_cdr2)
                peptide_matching.set_belonging(Belonging.CDR2)
                peptide_matching.set_coverage(coverage)
                peptide_matching.set_cdr_fragmentation_ions(cdr_fragments)
                peptide_matching.set_fragmentation_coverage(cdr_fragment_coverage)
                return True
            return False

        '''
        If can be aligned to CDR3
        '''
        overlap_span = find_overlap(self.domain_cdr3, _peptide_span)
        overlap = overlap_span[1] - overlap_span[0]

        if overlap > 0:
            coverage = overlap / len(self.domain_cdr3)
            peptide_matching.set_belonging(Belonging.CDR3)
            peptide_matching.set_coverage(coverage)
            if coverage >= Param.CDR3_PEPTIDE_COVERAGE:
                cdr_fragments = find_cdr_fragments(_peptide_span.start, overlap_span, peptide_fragmentation)
                cdr_fragment_coverage = cdr_fragments / len(self.domain_cdr3)
                '''
                updated 10/16/2019
                deal with peptides containing C-terminus
                discard peptides covering entire cdr3(big peptides population)
                For peptides containing only YYC, apply one fragmentation filter
                For peptides containing only QVTVS, apply another fragmentation filter.
                '''


                '''
                Identify big tryptic peptide covering entire CDR3 
                '''
                if Param.ENZYME == "Trypsin" and coverage == 1.0 and Y_CPattern.search(peptide_matching.peptide.sequence)!= None and WG_GPattern.search(peptide_matching.peptide.sequence)!= None:
                    return False

                '''
                Get the enzyme name specified
                '''
                _enzyme = Param.ENZYME

                '''
                Identify peptides containing C-terminus
                '''
                if Param.ENZYME == "Trypsin" and "QVTVS" in peptide_matching.peptide.sequence:
                    #change the enzyme name
                    _enzyme = "Trypsin_C_Term"

                if cdr_fragment_coverage >= get_fragmentation_coverage_threashold(len(self.domain_cdr3),_enzyme):
                    peptide_matching.set_cdr_fragmentation_ions(cdr_fragments)
                    peptide_matching.set_fragmentation_coverage(cdr_fragment_coverage)
                    return True
            return False

        '''
        Aligned to FR
        '''

        return False

    def assign_cdr3_type(self,type):
        self.cdr3_type = type
    def assign_cdr3_type_2(self,type):
        self.cdr3_type2 = type

    def remove_redundant_peptides(self):
        '''
        Remove shorter peptide that can be included by a longer peptide
        '''
        _peptides_matchings = sorted(self.peptides_matchings)
        _new_peptides_matchings = list()

        current_span_start = 0

        for _peptide_matching in _peptides_matchings:
            if len(_new_peptides_matchings) == 0:
                _new_peptides_matchings.append(_peptide_matching)
                current_span_start = _peptide_matching.span.start
                continue
            if current_span_start != _peptide_matching.span.start:
                if _peptide_matching.span.end<=_new_peptides_matchings[-1].span.end:
                    continue
                current_span_start = _peptide_matching.span.start
                _new_peptides_matchings.append(_peptide_matching)
            else:
                _new_peptides_matchings.pop()
                _new_peptides_matchings.append(_peptide_matching)


        self.peptides_matchings = _new_peptides_matchings





    def calculate_coverage(self):
        '''
        Calculate the coverage of the protein
        '''
        sequencemap = [0 for i in range(len(self.sequence))]

        self.remove_redundant_peptides()

        for _peptide_matching in self.peptides_matchings:
            # label as mapped
            for i in range(_peptide_matching.span.start, _peptide_matching.span.end):
                sequencemap[i] = 1


        _cdr1_start = self.domain_cdr1.start
        _cdr1_end = self.domain_cdr1.end
        cdr1map = sequencemap[_cdr1_start:_cdr1_end]
        self.cdr1_coverage = np.sum(cdr1map) / len(cdr1map)

        _cdr2_start = self.domain_cdr2.start
        _cdr2_end = self.domain_cdr2.end
        cdr2map = sequencemap[_cdr2_start:_cdr2_end]
        self.cdr2_coverage = np.sum(cdr2map) / len(cdr2map)

        _cdr3_start = self.domain_cdr3.start
        _cdr3_end = self.domain_cdr3.end
        cdr3map = sequencemap[_cdr3_start:_cdr3_end]
        self.cdr3_coverage = np.mean(cdr3map)

        self.cdr_coverage = (self.cdr1_coverage + self.cdr2_coverage + self.cdr3_coverage) / 3

        #Non essential peptides,only increase whole sequence coverage
        for _pep_span in self.nonessential_peps:
            for i in range(_pep_span[0],_pep_span[1]):
                sequencemap[i] = 1
        self.coverage = np.mean(sequencemap)

    def add_conditionwise_abundance(self, cdr3_conditionwise_abundance):

        self.cdr3_conditionwise_abundance = cdr3_conditionwise_abundance


    def get_peptides_to_quantify(self):
        cdr3_peptides=set()

        for _peptide_matching in self.peptides_matchings:
            if _peptide_matching.belonging == Belonging.CDR3:
                cdr3_peptides.add(_peptide_matching.peptide)
        return cdr3_peptides

    def simple_str(self):
        _str_protein = self.id + "," + self.sequence + "," + str(self.coverage) + "," + str(self.cdr_coverage) + ","
        _cdr1_sequence = self.sequence[self.domain_cdr1.start:self.domain_cdr1.end]
        _cdr2_sequence = self.sequence[self.domain_cdr2.start:self.domain_cdr2.end]
        _cdr3_sequence = self.sequence[self.domain_cdr3.start:self.domain_cdr3.end]

        _str_protein += _cdr1_sequence + "," + str(self.cdr1_coverage) + ","
        _str_protein += "#".join([str(_peptide_matching) for _peptide_matching in self.peptides_matchings if
                                  _peptide_matching.belonging == Belonging.CDR1]) + ","

        _str_protein += _cdr2_sequence + "," + str(self.cdr2_coverage) + ","
        _str_protein += "#".join([str(_peptide_matching) for _peptide_matching in self.peptides_matchings if
                                  _peptide_matching.belonging == Belonging.CDR2]) + ","

        _str_protein += _cdr3_sequence + "," + str(self.cdr3_coverage) + ","
        _str_protein += "#".join([str(_peptide_matching) for _peptide_matching in self.peptides_matchings if
                                  _peptide_matching.belonging == Belonging.CDR3]) + ","

        _str_protein += self.source

        _str_protein += "\n"

        return _str_protein

    def __str__(self):
        _str_protein = self.id + "," + self.sequence + "," +str(self.coverage)+"," +str(self.cdr_coverage) + ","
        _cdr1_sequence = self.sequence[self.domain_cdr1.start:self.domain_cdr1.end]
        _cdr2_sequence = self.sequence[self.domain_cdr2.start:self.domain_cdr2.end]
        _cdr3_sequence = self.sequence[self.domain_cdr3.start:self.domain_cdr3.end]

        _str_protein += _cdr1_sequence + "," + str(self.cdr1_coverage) + ","
        _str_protein += "#".join([str(_peptide_matching) for _peptide_matching in self.peptides_matchings if
                                  _peptide_matching.belonging == Belonging.CDR1]) + ","


        _str_protein += _cdr2_sequence + "," + str(self.cdr2_coverage) + ","
        _str_protein += "#".join([str(_peptide_matching) for _peptide_matching in self.peptides_matchings if
                                  _peptide_matching.belonging == Belonging.CDR2]) + ","


        _str_protein += _cdr3_sequence + "," + str(self.cdr3_coverage) + ","
        _str_protein += "#".join([str(_peptide_matching) for _peptide_matching in self.peptides_matchings if
                                  _peptide_matching.belonging == Belonging.CDR3]) + ","

        _str_protein += ",".join(self.cdr3_conditionwise_abundance.astype(str)) + ","

        _str_protein += str(self.cdr3_type) +"," +str(self.cdr3_type2)+","



        #_str_protein += ",".join(self.weighted_conditionwise_abundance.astype(str)) + ","
        #_str_protein += ",".join(self.best_frag_conditionwise_abundance.astype(str))+","

        _str_protein += self.source



        _str_protein+= "\n"

        return _str_protein

















