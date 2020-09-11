from Base import *
from Experiment import Experiemnt
from Peptide import Peptide,PeptideMatching,PeptideFragments

from Config import Param



class PeptideManager(object):
    def __init__(self):
        '''
        Initializer for peptide manager
        It contains the peptides collection and a cleavage manager
        '''
        self.peptides = dict()

    def load_peptides_from_pd(self,experiment:Experiemnt):
        '''
        Load peptide identification result from proteome discover
        Apply some filters on identification results
        Fill survived peptides into peptides collection with key information reserved
        :param fpath:
        :return:
        '''
        data = pd.read_excel(experiment.peptide_identification_file)

        # ignore post-translational modifications
        # turn sequence into upper case
        # remove 'Rejected' PSMs
        # remove duplicate sequences

        data['Sequence'] = data['Annotated Sequence'].apply(lambda seq: seq.upper())

        data['PSM Ambiguity'] = data[data['PSM Ambiguity'] != 'Rejected']

        data = data.sort_values(['Sequence', 'Charge', 'Isolation Interference [%]']).drop_duplicates(['Annotated Sequence','Charge'], 'first')

        for index,row in data.iterrows():
            _sequence = row['Sequence']
            _mz = round(float(row['m/z [Da]']),4)
            _charge = int(row['Charge'])
            _retention_time = round(float(row['RT [min]']),4)
            _PEP = float(row['Percolator PEP'])
            _qvalue = float(row['Percolator q-Value'])
            if _sequence in self.peptides.keys():
                self.peptides[_sequence].additional.append(Peptide(sequence = _sequence,mz = _mz,charge = _charge, retention_time = _retention_time, PEP = _PEP,qvalue=_qvalue,index=experiment.index))
            else:
                self.peptides[_sequence]=Peptide(sequence = _sequence,mz = _mz,charge = _charge, retention_time = _retention_time, PEP = _PEP,qvalue=_qvalue,index=experiment.index)

    def load_peptide_fragments_annotation(self,experiment:Experiemnt):
        dpath = experiment.fragment_spectra_folder
        fragments_dict=dict()
        for file in os.listdir(dpath):
            try:
                if not ".txt" in file:
                    continue
                fullpath = os.path.join(dpath,file)
                _seq = file.split("_")[0]
                if not _seq in self.peptides.keys():
                    continue

                _fragments = np.zeros(len(_seq),int)
                data = pd.read_csv(fullpath, sep="\t")
                data = data[~data['matches'].isnull()]
                for hit in data['matches']:
                    if 'NH3' in hit or 'H2O' in hit:
                        continue
                        # in case there are multiple matches sharing same m/z
                        # example: y (2) (2+), b (10) (3+)

                        # we will take the first one
                        # example: y (2) (2+)
                    ion = hit.split(',')[0]

                        # remove information about charge/ loss of water/amonia
                        # example:    ["y","(2)"]
                        #      [0]ion type   [1]fragmented postion
                    ion = ion.strip().split(' ')[:2]

                    if ion[0] == 'b':
                        pos = int(ion[1][1:-1])
                        _fragments[pos - 1] = 1
                    elif ion[0] == 'y':
                        pos = int(ion[1][1:-1])
                        _fragments[len(_seq) - pos] = 1
                    else:
                        continue
                if _seq in fragments_dict.keys():
                    merged_fragments = _fragments | fragments_dict[_seq]
                    fragments_dict[_seq] = merged_fragments
                else:
                    fragments_dict[_seq] = _fragments
            except:
                raise IOError("Error loading fragment file:{}".format(os.path.join(dpath,file)))
            for _seq,_fragments in fragments_dict.items():
                self.peptides[_seq].annotate_fragment_ions(PeptideFragments(_fragments))


    # def load_peptide_fragments_annotation(self,experiment:Experiemnt):
    #     '''
    #     Provided the directory path containing spectra ion match information generated from proteome discover
    #     parse peptides fragments information
    #     :param dpath:
    #     :return:
    #     '''
    #     dpath = experiment.fragment_spectra_folder
    #     for peptide in self.peptides.values():
    #         _seq = peptide.sequence
    #         _fragments=[ 0 for i in range(len(_seq))]
    #
    #         #
    #         # Merge fragments distribution for the identical peptide
    #         #
    #         files = [os.path.join(dpath,file) for file in os.listdir(dpath) if _seq == file.split("_")[0]]
    #         for file in files:
    #             data = pd.read_csv(file,sep="\t")
    #             data = data[~data['matches'].isnull()]
    #             for hit in data['matches']:
    #                 #
    #                 # in case there are multiple matches sharing same m/z
    #                 # example: y (2) (2+), b (10) (3+)
    #
    #                 # we will take the first one
    #                 # example: y (2) (2+)
    #                 ion = hit.split(',')[0]
    #
    #                 # remove information about charge/ loss of water/amonia
    #                 # example:    ["y","(2)"]
    #                 #      [0]ion type   [1]fragmented postion
    #                 ion = ion.strip().split(' ')[:2]
    #
    #                 if ion[0] == 'b':
    #                     pos = int(ion[1][1:-1])
    #                     _fragments[pos - 1] = 1
    #                 elif ion[0] == 'y':
    #                     pos = int(ion[1][1:-1])
    #                     _fragments[len(_seq) - pos] = 1
    #                 else:
    #                     continue
    #         self.peptides[_seq].annotate_fragment_ions(_fragments)

    def get_peptide_fragmentation(self,peptide_sequence:str):
        if peptide_sequence not in self.peptides.keys():
            warnings.warn("{} not included,skipped")
            return None
        return self.peptides[peptide_sequence].fragments



    def find_parent_protein_id_for_peptides(self,db_file_path:str):
        '''
        Provided with entire protein database
        In-silico digest each protein and
        find out proteins that are relevant to peptides that are identified
        :param db_file_path:
        :return:
        '''


        index_file_path = CleavageManager.get_peptides_index_file(db_file_path)
        db_entries = 0
        for record in SeqIO.parse(index_file_path,'fasta'):
            db_entries=+1
            _id = str(record.name)
            _peptides = set(str(record.seq).split("#"))

            # find whether there is any peptide that is identified comes from current protein
            _peptides_set =  self.peptides.keys() & _peptides

            for peptide in _peptides_set:
                self.peptides[peptide].add_protein_id(_id)

        '''
        Remove peptides that can be mapped 20% sequences in database 
        '''
        '''
        for _peptide in list(self.peptides.keys()):
            if len(self.peptides[_peptide].proteins_ids) > 5/db_entries:
                self.peptides.pop(_peptide)
                
        '''


    def output_peptides_information(self):
        output = Param.TASK_NAME + "_peptide_information.csv"

        with open(output,'w+') as out:
            out.write("peptide,m/z,charge,retention time,PEP,q-value,num(Proteins)\n")
            for peptide in self.peptides.values():
                out.write(str(peptide))









class CleavageManager(object):

    @staticmethod
    def verify_index(db_file):

        dirname = os.path.dirname(db_file)
        basename = os.path.basename(db_file)

        index_file = os.path.join(dirname, "_".join(
            [os.path.splitext(basename)[0], Param.ENZYME, str(Param.MISS_CLEAVAGE_ALLOWED)]) + ".fasta")

        if not verify_file_path(os.path.join(dirname,index_file)):
            print("Peptide index for {} does not exist, build index now".format(db_file))
            CleavageManager.build_index(db_file,index_file)
            print("Finished building index. Location:{}".format(index_file))

    @staticmethod
    def get_peptides_index_file(db_file):
        dirname = os.path.dirname(db_file)
        basename = os.path.basename(db_file)

        index_file = os.path.join(dirname,"_".join([os.path.splitext(basename)[0], Param.ENZYME, str(Param.MISS_CLEAVAGE_ALLOWED)])+".fasta")


        return index_file

    @staticmethod
    def build_index(db_file:str,index_file:str):
        tmp_file = index_file+".tmp"
        index_writer = open(tmp_file,"w+")
        for record in SeqIO.parse(db_file,'fasta'):
            _id = str(record.id)
            _seq = str(record.seq)
            _peptides = CleavageManager.cleave(_seq)

            index_writer.write(">"+_id+"\n")
            index_writer.write("#".join(_peptides)+"\n")
        index_writer.close()

        os.rename(tmp_file,index_file)


    @staticmethod
    def set_miscleavage(ml:int):
        '''
        Override the default setting of miscleavages
        :param ml:
        :return:
        '''
        if ml <0 :
            warnings.warn("miss cleavage cannot be zero, nothing changes")
            return False
        MISS_CLEAVAGE_ALLOWED = ml

    @staticmethod
    def cleave(sequence:str):
        '''
        Cleave sequences into peptides
        :param sequence:
        :return peptides set:
        '''
        peptides = set()

        cleavage_pos = [0]
        for match in re.finditer(Param.CLEAVAGE_RULES[Param.ENZYME], sequence):
            cleavage_pos.append(match.span(0)[1])
        cleavage_pos.append(len(sequence))

        # get cleaved peptides
        for miscleavage in range(Param.MISS_CLEAVAGE_ALLOWED + 1):
            for i in range(len(cleavage_pos) - 1 - miscleavage):
                _start = cleavage_pos[i]
                _end = cleavage_pos[i + 1 + miscleavage]

                if _end - _start <4:
                    continue
                _pep = sequence[_start:_end]
                # _mass = mass.calculate_mass(sequence=_pep)
                # if  _mass < self.__min_mass or _mass > self.__max_mass:
                #     continue



                peptides.add(_pep)

        return peptides


