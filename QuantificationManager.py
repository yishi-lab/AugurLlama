from Base import *
from ProteinManager import Protein
from Peptide import Belonging
from functools import lru_cache
from pymsfilereader import MSFileReader
from scipy.cluster.hierarchy import ward, fcluster
from scipy.spatial.distance import pdist
from scipy.cluster import hierarchy
from scipy.stats import zscore
from Config import Param
from scipy.ndimage import gaussian_filter1d,maximum_filter1d
from scipy.signal import find_peaks, peak_widths
from tqdm import tqdm

from Experiment import Experiemnt

class RetentionTimeShiftArray(object):

    def __init__(self,median_rt_shift:np.array):
        self.median_rt_shift= median_rt_shift

    def __repr__(self):
        return str(self.median_rt_shift)

    def get_inferred_rt(self,experiment_target_index:int,experiemnt_reference_index:int,reference_rt:float):
        '''

        :param experiment_target_index:
        :param experiemnt_reference_index:
        :param reference_rt: peptide retention time identified in reference experiment
        :return:
        '''
        reversed = 1
        if experiment_target_index > experiemnt_reference_index:
            reversed = -1

        index = experiment_target_index+experiemnt_reference_index-1
        delta_rt = self.median_rt_shift[index]


        return reference_rt+delta_rt*reversed


class PeptideIon(object):
    def __init__(self, seq,mz, charge):
        self.seq = seq
        self.charge = charge
        self.mz = mz

    def add_rt_array(self,experiments_num:int):
        '''
        Used to determine NaN retention time among experiments in quantification
        initialize empty retention time matrix
        :param experiments_num: number of experiments
        :return:
        '''
        self.rt_inferred = np.zeros(experiments_num,int)
        self.rt_array          = np.zeros(experiments_num,float)


    def add_belonging(self,belonging):
        self.belonging = belonging


    def add_rt(self,experiment_index:int,rt:float):
        '''
        Add a retention time value to a specific condition-replicate cell
        :param experiment_index:
        :param rt:
        :return:
        '''
        self.rt_array[experiment_index] = rt


    def add_abundance_array(self,abundance_array:np.array):
        '''
        After finishing quantification, fill the data into object
        :param abundance_array:
        :return:
        '''
        self.abundance_array = abundance_array

    def add_zscores(self,zscores):
        '''
        Add zscores array generated from mean abundance array across 3 conditions
        :param zscores:
        :return:
        '''
        self.zscores = zscores

    def assign_type(self,type:int):
        self.type = type

    def assign_type_2(self,type:int):
        self.type_2 = type


    def add_conditionwise_abundance(self,conditionwise_abundance):
        self.conditionwise_abundance = conditionwise_abundance

    def __hash__(self):
        return hash((self.seq, self.charge))

    def __ne__(self, other):
        return not self.__eq__(other)

    def __eq__(self, other):
        return self.seq == other.seq and self.charge == other.charge

    def __str__(self):
        _str_peptide_ion = self.seq+"," + str(self.belonging) +","+str(self.mz)+","+str(self.charge)+","
        _str_peptide_ion += ",".join(self.rt_array.astype(str)) +","
        _str_peptide_ion += ",".join(self.rt_inferred.astype(bool).astype(str)) + ","
        _str_peptide_ion += ",".join(self.abundance_array.astype(str)) +","
        _str_peptide_ion += ",".join(self.conditionwise_abundance.astype(str))+","
        _str_peptide_ion += ",".join(self.zscores.astype(str))+","
        _str_peptide_ion += str(self.type) +"," +str(self.type_2)
        _str_peptide_ion +="\n"
        return _str_peptide_ion

    def __repr__(self):
        return "seq:{},mz:{},charge:{},rt:{},abudance:{}".format(self.seq,str(self.mz), str(self.charge),str(self.rt_array),str(self.abundance_array))



class ExperimentManager(object):
    def __init__(self):
        self.experiments_num  = 0
        self.experiments = list()
        self.read_experiments_conf()

    def initialize(self):
        '''
        Initialization process
        Read in quantification setting
        Recalibrate the retention time by calculating
        the general shift through constant peptides(can be identified in all conditions)
        :return:
        '''

        if Param.QUANTIFICATION:
            self.generate_rt_shift_array()


    def read_experiments_conf(self):

        try:
            rawfile_folder = Param.RAW_DIR
            peptide_folder = Param.PEP_DIR
            fragment_spectra_folder= Param.FRAG_SPEC_DIR

            for _obj in Param.EXPERIMENTS_SET:
                _index,_raw_file,_id_file, _condition_num = _obj
                raw_file_full_path = os.path.join(rawfile_folder, _raw_file)
                id_file_full_path = os.path.join(peptide_folder, _id_file)
                experiment = Experiemnt(identification_file=id_file_full_path, raw_ms_file=raw_file_full_path,
                                        fragment_spectra_folder=fragment_spectra_folder,
                                        index=int(_index), condition=int(_condition_num))
                self.experiments.append(experiment)
            self.experiments.sort()
            self.experiments_num = len(self.experiments)
        except:
            raise AttributeError("Error in reading experiment setting")



    def get_shared_peptides(self,experiments:list):
        _peptide_id_1 = pd.read_excel(experiments[0].peptide_identification_file)
        _peptide_id_2 = pd.read_excel(experiments[1].peptide_identification_file)

        _peptide_id_1['sequence'] = _peptide_id_1['Annotated Sequence'].apply(lambda x: x.upper())
        _peptide_id_2['sequence'] = _peptide_id_2['Annotated Sequence'].apply(lambda x: x.upper())

        _peptide_id_1['m/z'] = _peptide_id_1['m/z [Da]'].apply(lambda x: round(float(x),4))
        _peptide_id_1['charge'] = _peptide_id_1['Charge'].apply(lambda x: int(x))
        _peptide_id_2['m/z'] = _peptide_id_2['m/z [Da]'].apply(lambda x: round(float(x), 4))
        _peptide_id_2['charge'] = _peptide_id_2['Charge'].apply(lambda x: int(x))

        _peptide_id_1 = _peptide_id_1.sort_values(['sequence', 'charge', 'Isolation Interference [%]']).drop_duplicates(['sequence', 'charge'], 'first')
        _peptide_id_2 = _peptide_id_2.sort_values(['sequence', 'charge', 'Isolation Interference [%]']).drop_duplicates( ['sequence', 'charge'], 'first')

        joined = _peptide_id_1.merge(_peptide_id_2, on=["sequence", 'charge'])

        joined = joined[['RT [min]_x', 'RT [min]_y']]

        return joined
    '''
    def get_shared_peptides(self,experiments:list):
        sharedpeptides = None

        for experiment in experiments:
            data = pd.read_excel(experiment.peptide_identification_file)
            data['sequence'] = data['Annotated Sequence'].apply(lambda x: x.upper())
            pepset = set()
            for i, row in data.iterrows():
                p = PeptideIon(row['sequence'],round(float(row['m/z [Da]']),4), int(row['Charge']))
                pepset.add(p)
            if sharedpeptides == None:
                sharedpeptides = pepset
            else:
                sharedpeptides = sharedpeptides.intersection(pepset)

        shared = pd.DataFrame([[s.seq, s.charge] for s in sharedpeptides], columns=['peptide', 'charge'])
        for index in range(len(experiments)):
            experiment = experiments[index]
            data = pd.read_excel(experiment.peptide_identification_file)
            data['sequence'] = data['Annotated Sequence'].apply(lambda x: x.upper())
            data['Charge'] = data['Charge'].astype(int)
            data = data.sort_values(['sequence', 'Charge', 'RT [min]']).drop_duplicates(['sequence', 'Charge'],'first')
            data = data.set_index(['sequence','Charge'])
            shared[index] = shared.apply(lambda row: data.loc[(row['peptide'],row['charge'])]['RT [min]'],axis=1)
        return shared
    '''


    def generate_rt_shift_array(self):
        shift = list()
        for i in range(len(self.experiments) - 1):
            for j in range(i+1, len(self.experiments)):
                shared_peptide_data = self.get_shared_peptides([self.experiments[i],self.experiments[j]])
                shift_median = np.median(shared_peptide_data['RT [min]_x'] - shared_peptide_data['RT [min]_y'])
                shift.append(shift_median)

        self.rt_shift =  RetentionTimeShiftArray(shift)


    def get_inferred_rt(self, experiment_target_index: int, experiemnt_reference_index: int, reference_rt: float):

        return self.rt_shift.get_inferred_rt(experiment_target_index,experiemnt_reference_index,reference_rt)

    def fill_rt_array_for_peptide(self,peptide:PeptideIon):
        #fill missing rt
        new_rt_array = np.copy(peptide.rt_array)
        new_rt_inferred = np.copy(peptide.rt_inferred)
        for target_experiment_index in np.where(peptide.rt_array ==0)[0]:
            inferred_rts = list()
            for reference_experiment_index in np.nonzero(peptide.rt_array)[0]:
                rt = self.get_inferred_rt(target_experiment_index, reference_experiment_index,
                                          peptide.rt_array[reference_experiment_index])
                inferred_rts.append(rt)
            new_rt_array[target_experiment_index] = np.median(inferred_rts)
            new_rt_inferred[target_experiment_index] = 1

        #update info in peptide
        peptide.rt_array = new_rt_array
        peptide.rt_inferred = new_rt_inferred





class QuantificationManager(object):
    def __init__(self):
        '''
        Initilize experiment_manager and calculate general retention time shift pattern
        '''

        self.peptides = dict()
        self.em = None

    def register_experiment_manager(self,em:ExperimentManager):
        self.em = em
        self.em.generate_rt_shift_array()

    def load_raw_files(self):
        experiments = self.em.experiments

        experiments.sort()


        readers = list()

        for experiment in experiments:
            raw_file_reader = MSFileReader(experiment.raw_ms_file)

            readers.append(raw_file_reader)
        self.raw_file_readers = np.array(readers)

    def close_raw_files(self):
        for raw_file_reader in self.raw_file_readers.flatten():
            raw_file_reader.Close()






    def add_peptide(self,peptide:str,mz:float,charge:int,rt:float,index:int,belonging:Belonging):
        '''
        Add CDR3 peptides that needed  to be quantified
        :param peptide:
        :param mz:
        :param charge:
        :param rt:
        :param index:
        :return:
        '''
        if not peptide in self.peptides.keys():
            self.peptides[peptide] = list()
        is_set = False
        for peptide_ion in self.peptides[peptide]:
            if peptide_ion.charge == charge:
                peptide_ion.add_rt(index,rt)
                is_set = True
                break
        if not is_set:
            peptide_ion = PeptideIon(peptide,mz,charge)
            peptide_ion.add_belonging(belonging)
            peptide_ion.add_rt_array(self.em.experiments_num)
            peptide_ion.add_rt(index,rt)
            self.peptides[peptide].append(peptide_ion)

    def fill_rt_array_for_peptides(self):
        '''
        Fill the missing RT for peptide
        :return:
        '''
        for peptide_ion_list in self.peptides.values():
            for peptide_ion in peptide_ion_list:

                self.em.fill_rt_array_for_peptide(peptide_ion)

    def quantify_peptides(self):
        pbar = tqdm(desc="Peptides quantification",total=len(self.peptides))
        for peptide_ion_list in self.peptides.values():
            for peptide_ion in peptide_ion_list:
                abudances=list()
                for index in range(self.em.experiments_num):
                    mz = peptide_ion.mz
                    charge = peptide_ion.charge
                    rt = peptide_ion.rt_array[index]
                    rt_inferred = peptide_ion.rt_inferred[index]
                    abudances.append(self.calculate_abundance(mz, charge, rt, rt_inferred,index))

                abudances = np.array(list(abudances))
                peptide_ion.add_abundance_array(abudances)
                self.calculate_conditionwise_abundance(peptide_ion)
            pbar.update(1)
        pbar.close()

    def generate_xic(self,mz:float,charge:int,rt:float,rt_inferred:int,index:int):
        '''
        get extracted intensity chromotogram
        :param mz: mass-to-charge ratio of the peptide
        :param charge: charge state of the peptide
        :param rt: retention time
        :param rt_inferred: [1,0] wheter the retention time is inferred(determine the XIC window width)
        :param index:
        :return:
        '''
        mz_min = round(mz - Param.QUANTIFICATION_MS1_ION_MZ_TOLERANCE * mz / 1000000, 4)
        mz_max = round(mz + Param.QUANTIFICATION_MS1_ION_MZ_TOLERANCE * mz / 1000000, 4)



        raw_file_reader = self.raw_file_readers[index]
        '''
        ensure retention time is within the LC time
        '''
        if rt < raw_file_reader.GetStartTime():
            rt = 1.0
        if rt > raw_file_reader.GetEndTime():
            rt = raw_file_reader.GetEndTime()-1

        rt_start = rt - (0.25 if rt_inferred ==0 else 2)

        rt_end = rt + (0.25 if rt_inferred ==0 else 2)

        if rt_start <raw_file_reader.GetStartTime() :
            rt_start = raw_file_reader.GetStartTime()+0.1
            rt_end = rt_start +4
        if rt_end >raw_file_reader.GetEndTime():
            rt_end = raw_file_reader.GetEndTime() - 0.1
            rt_start = rt_end -4

        start_scan = raw_file_reader.ScanNumFromRT(rt_start)
        end_scan = raw_file_reader.ScanNumFromRT(rt_end)


        RTs = []
        Intensities = []
        for scan in range(start_scan,end_scan+1):
            has = False
            scan += 1
            RT = raw_file_reader.RTFromScanNum(scan)
            if raw_file_reader.GetMSOrderForScanNum(scan) == 1:
                label = raw_file_reader.GetLabelData(scan)
                '''
                labels:      mass (double),
                             intensity (double),
                             resolution (float),
                             baseline (float),
                             noise (float)
                             charge (int)
                '''

                mzs = label[0][0]
                intensities = label[0][1]
                charges = label[0][5]

                for i in range(len(mzs)):
                    _mz = mzs[i]
                    _intensity = intensities[i]
                    _charge = charges[i]
                    if (_mz <= mz_max and _mz >= mz_min and _charge == charge):
                        has = True
                        Intensities.append(_intensity)
                        RTs.append(RT)
                if not has:
                    Intensities.append(0)
                    RTs.append(RT)

        return [RTs, Intensities]

    def calculate_abundance(self,mz:float,charge:int,rt:float,rt_inferred:int,index:int):
        chromotogram = self.smooth(self.generate_xic(mz,charge,rt,rt_inferred,index))
        if rt_inferred == 0:
            return  self.area_integration(chromotogram)
        else:
            quantifiable, narrowed_chromotograms = self.detect_peaks(chromotogram)

            if not quantifiable:
                return np.nan
            if narrowed_chromotograms == None:
                return 0

            abundances = [self.area_integration(chrom) for chrom in narrowed_chromotograms]
            return np.mean(abundances)



    def detect_peaks(self,chromotogram:list):



        _peaks,_properties = find_peaks(chromotogram[1], prominence=1,width=2)

        if len(_peaks) == 0:
            return [True,None]

        _peak_width = peak_widths(chromotogram[1],_peaks,rel_height=1)

        peak_regions=list()

        indices = sorted(np.array(_peak_width[2:4]).flatten())
        for _peak in _peaks:
            peak_regions.append(self.find_peak_width(_peak,indices))

        upper_limit_abundance = np.max(_properties['prominences'])
        lower_limit_abundance = 0.1 * upper_limit_abundance

        interference_peaks_index = [i for i in range(len(_peaks)) if _peaks[i]>= lower_limit_abundance and _peaks[i]<= upper_limit_abundance]

        if len(interference_peaks_index)>=5:
            #Not quantifiable
            return [False,None]

        narrowed_chromotograms = [
            [
                chromotogram[0][peak_region[0]:peak_region[1]+1],
                chromotogram[1][peak_region[0]:peak_region[1]+1]
            ]
            for peak_region in peak_regions]

        return [True,narrowed_chromotograms]

    def find_peak_width(self,peak,peak_indices):

            left = 0
            right = len(peak_indices)-1

            while left <= right:
                mid = left + int(( right - left)/2)
                if peak_indices[mid] < peak:
                    left=mid+1
                elif peak_indices[mid] > peak:
                    right=mid-1
                else:
                    return [peak_indices[mid],peak_indices[mid+1]]
            return [int(round(peak_indices[left-1],0)),int(round(peak_indices[left],0))]


    def area_integration(self,chromotogram:list):
        area = 0

        for i in range(len(chromotogram[1]) - 1):
            diff_rt = chromotogram[0][i + 1] - chromotogram[0][i]
            start_int = chromotogram[1][i]
            end_int = chromotogram[1][i + 1]
            area += 60 * diff_rt * (end_int + start_int) / 2
        return area

    def smooth(self,chromotogram:list):
        maxi =maximum_filter1d(chromotogram[1],3)
        gauss = gaussian_filter1d(maxi,2)
        return [chromotogram[0],gauss]

    def calculate_conditionwise_abundance(self,peptide:PeptideIon):
        condition_bins = np.diff([experiment.condition for experiment in self.em.experiments]).nonzero()[0]+1
        abundance_groups = np.split(peptide.abundance_array,condition_bins)
        mean_abundance_array = np.array([abundances.mean() for abundances in abundance_groups])

        peptide.add_conditionwise_abundance(mean_abundance_array)



    @lru_cache(maxsize=500)
    def get_peptide_abudance(self,peptide_sequence:str):
        conditionwise_abundances = np.zeros(3,float)

        for peptide_ion in self.peptides[peptide_sequence]:
            conditionwise_abundances += peptide_ion.conditionwise_abundance

        return conditionwise_abundances

    def generate_abundance_for_protein(self,protein:Protein):
        cdr3_conditionwise_abundances = np.zeros(3,float)

        for _peptide_matching in protein.peptides_matchings:
            if _peptide_matching.belonging == Belonging.CDR3:
                cdr3_conditionwise_abundances += self.get_peptide_abudance(_peptide_matching.peptide.sequence)



        protein.add_conditionwise_abundance(cdr3_conditionwise_abundances)


    def classify_protein(self,protein:Protein):
        type_set= set()

        for _peptide_matching in protein.peptides_matchings:
            if _peptide_matching.belonging == Belonging.CDR3:
                _type = self.peptide_type[_peptide_matching.peptide.sequence]
                type_set.add(_type)

        if len(type_set)==0:
            print("None CDR3 peptides :{}".format(str(protein.peptides_matchings)))
            protein.assign_cdr3_type(None)
        elif len(type_set) >1:
            protein.assign_cdr3_type(-1)
        else:
            protein.assign_cdr3_type(type_set.pop())

    def classify_protein_2(self,protein:Protein):
        type_set= set()

        for _peptide_matching in protein.peptides_matchings:
            if _peptide_matching.belonging == Belonging.CDR3:
                _type = self.peptide_type_2[_peptide_matching.peptide.sequence]
                type_set.add(_type)

        if len(type_set)==0:
            print("Sequence:{}\n".format(protein.sequence))
            print("None CDR3 peptides :{}\n".format(str(protein.peptides_matchings)))
            print("Protein CDR3:{}-{}\n".format(protein.domain_cdr3.start,protein.domain_cdr3.end))
            protein.assign_cdr3_type_2(None)
        elif len(type_set) >1:
            protein.assign_cdr3_type_2(-1)
        else:
            protein.assign_cdr3_type_2(type_set.pop())

    def is_decreasing_zscores(self,zscores):
            y1, y2, y3 = zscores

            if y1 == np.nan or y2 == np.nan or y3 == np.nan:
                return -1
            if y1 == y2 and y2 == y3:
                return 0

            b = y2 - y1
            c = y3 - y1

            if b > 0 and c >= 0 or b >= 0 and c > 0:
                return 1
            elif b <= 0 and c <= 0:
                return 0
            elif b * c < 0:
                if b < 0:
                    return int(abs(b) < (c) / 2)
                else:
                    return int((b) > abs(c) / 2)
            return 0

    def annotate_cluster(self,row):
        c1 = row['z1']
        c2 = row['z2']
        c3 = row['z3']
        if c3 >= c2 and c2 > c1:
            return 2
        return 1

    def classify_peptides_2(self):
        FOLD = 1
        self.peptide_type_2 = dict()
        for _peptide_ion_list in self.peptides.values():
            _peptide_seq = _peptide_ion_list[0].seq
            ab_1,ab_2,ab_3 = self.get_peptide_abudance(_peptide_seq)

            flags = []
            type = 0

            if ab_1 > FOLD*(ab_2 + ab_3):
                flags.append(0)
            if ab_2 > FOLD*(ab_1 + ab_3):
                flags.append(1)
            if ab_3 > FOLD*(ab_1 + ab_2) :
                flags.append(2)
            if len(flags) == 2:
                if 1 in flags and 2 in flags:
                    type =  2
            if len(flags) == 1:
                type = flags[0]

            self.peptide_type_2[_peptide_seq] = type
            for _peptide_ion in _peptide_ion_list:

                _peptide_ion.assign_type_2(type)

    def classify_peptides(self):
        decreasing_peptides_sequences = list()
        nan_peptides_sequence = list()
        to_cluster_peptides= list()
        for _peptide_ion_list in self.peptides.values():
            _peptide_seq = _peptide_ion_list[0].seq
            zscore_abundance_array = zscore(self.get_peptide_abudance(_peptide_seq), ddof=2)
            for _peptide_ion in _peptide_ion_list:
                _peptide_ion.add_zscores(zscores=zscore_abundance_array)
                ret = self.is_decreasing_zscores(zscore_abundance_array)
            if  ret== 0:
                decreasing_peptides_sequences.append(_peptide_seq)
            elif ret == -1:
                nan_peptides_sequence.append(_peptide_seq)
            else:
                record = [_peptide_seq,zscore_abundance_array[0],zscore_abundance_array[1],zscore_abundance_array[2]]
                to_cluster_peptides.append(record)

        peptides_data = pd.DataFrame(to_cluster_peptides,columns=['sequence','z1','z2','z3'])
        print(peptides_data)
        peptides_data.index = peptides_data['sequence']
        peptides_data = peptides_data[['z1','z2','z3']]
        dist = pdist(peptides_data[['z1', 'z2', 'z3']], metric='correlation')
        print(dist)
        hcluster = hierarchy.linkage(dist, 'ward')
        peptides_data['cluster_id'] = fcluster(hcluster, 6, 'maxclust')
        peptides_data =peptides_data.groupby('cluster_id')[['z1','z2','z3']].transform(lambda x: x.mean())
        peptides_data['label'] = peptides_data.apply(lambda row: self.annotate_cluster(row), axis=1)
        #peptides_data.drop('cluster_id',axis=1, inplace=True)

        self.peptide_type = dict()

        for _peptide_seq in decreasing_peptides_sequences:
            self.peptide_type[_peptide_seq] = 0
            for _peptide_ion in self.peptides[_peptide_seq]:
                _peptide_ion.assign_type(0)
        for _peptide_seq in nan_peptides_sequence:
            self.peptide_type[_peptide_seq] = -1
            for _peptide_ion in self.peptides[_peptide_seq]:
                _peptide_ion.assign_type(-1)

        for _peptide_seq, row in peptides_data.iterrows():
            self.peptide_type[_peptide_seq] = row['label']
            for _peptide_ion in self.peptides[_peptide_seq]:
                _peptide_ion.assign_type(row['label'])



    def output_peptide_quantification_table(self):

        output_name = Param.TASK_NAME + "_peptide_quantification.csv"
        with open(output_name,"w+") as fout:
            fout.write("peptide,belong, m/z,charge,")
            raw_names = [os.path.basename(experiment.raw_ms_file) for experiment in self.em.experiments]
            fout.write(",".join(["RT_" + raw_name for raw_name in raw_names]))
            fout.write(",")
            fout.write(",".join(["RT_Inferred_" + raw_name for raw_name in raw_names]))
            fout.write(",")
            fout.write(",".join(["Abundance_" + raw_name for raw_name in raw_names]))
            fout.write(",")
            fout.write("Condition 1,Condition 2 ,Condition3,")
            fout.write("Z-score 1,Z-score 2,Z-score 3,")
            fout.write("type,type_2\n")
            for peptide_ion_list in self.peptides.values():
                for peptide_ion in peptide_ion_list:
                    fout.write(str(peptide_ion))

























if __name__ == "__main__":
    em = ExperimentManager()
    em.initialize()
    print(em.rt_shift)