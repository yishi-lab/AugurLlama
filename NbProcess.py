from PeptideManager import PeptideManager
from ProteinManager import Protein
from Peptide import Belonging
from Experiment import Experiemnt
from QuantificationManager import ExperimentManager,QuantificationManager
from PeptideManager import CleavageManager
import sys
import logging
from functools import partial
import tempfile
import pickle
from Base import *
from Config import Param,read_params

class Processor(object):
    def __init__(self,experiment:Experiemnt,param_file:str):
        self.peptide_m = PeptideManager()
        self.proteins = list()
        self.experiment = experiment
        read_params(param_file)


    def assembly(self,fpath:str):
        '''
        Protein assembly
        Determine regions that peptide maps
        Calculate coverage
        Collect CDR3 peptides need to be quantified
        :param fpath:
        :return cdr3_peptides_set:
        '''
        num=0
        skipped_sequence_length =0
        skipped_cdr3_length = 0
        skipped_non_cdr3= 0
        peptides_set_to_quantify=set()
        for record in SeqIO.parse(fpath,'fasta'):

            _id = str(record.id)
            _seq = str(record.seq)

            if len(_seq) < Param.MIN_SEQUENCE_LEN+38:
                skipped_sequence_length+=1
                continue

            _protein = Protein(_id,_seq,self.experiment.peptide_identification_file)


            if len(_protein.domain_cdr3)<Param.MIN_CDR3_LEN:
                skipped_cdr3_length+=1
                continue

            is_contain_cdr3_peptides = False

            for _peptide in self.peptide_m.peptides.values():
                _peptide_sequence = _peptide.sequence

                if _id in _peptide.proteins_ids:
                    _peptide_fragmentation = self.peptide_m.get_peptide_fragmentation(_peptide_sequence)
                    if _peptide_fragmentation == None:
                        logging.error(
                            msg="[{}] None peptide fragmentation object:{}".format(os.path.basename(self.experiment.peptide_identification_file),_peptide_sequence))
                        continue
                    _added,belonging = _protein.add_peptide(_peptide,_peptide_fragmentation)
                    if _added and not is_contain_cdr3_peptides and belonging == Belonging.CDR3:
                        is_contain_cdr3_peptides = True

            if not is_contain_cdr3_peptides:
                skipped_non_cdr3+=1
                continue
            _protein.calculate_coverage()

            cdr3_peptides = _protein.get_peptides_to_quantify()


            peptides_set_to_quantify.update(cdr3_peptides)
            self.proteins.append(_protein)
        logging.info(
                msg="[{}] {} proteins skipped with length < {} aas".format(os.path.basename(self.experiment.peptide_identification_file),
                                                         skipped_sequence_length,Param.MIN_SEQUENCE_LEN))
        logging.info(
                msg="[{}] {} proteins skipped with cdr3 length < {} aas".format(os.path.basename(self.experiment.peptide_identification_file),
                                                         skipped_cdr3_length,Param.MIN_CDR3_LEN))

        logging.info(
                msg="[{}] {} proteins skipped with no cdr3 coverage".format(os.path.basename(self.experiment.peptide_identification_file),
                                                         skipped_non_cdr3))

        logging.info(msg="[{}] {} proteins identified".format(os.path.basename(self.experiment.peptide_identification_file),len(self.proteins)))
        return peptides_set_to_quantify

    def process(self):


        start = time.time()
        self.peptide_m.load_peptides_from_pd(self.experiment)
        end = time.time()
        logging.info(msg = "[{}]enzyme:{}],miscleav:{}".format(os.path.basename(self.experiment.peptide_identification_file),Param.ENZYME,Param.MISS_CLEAVAGE_ALLOWED))
        logging.info(msg="[{}]load peptide: ".format(os.path.basename(self.experiment.peptide_identification_file),len(self.peptide_m.peptides))+time.strftime("%H:%M:%S", time.gmtime(end - start)))
        sys.stdout.flush()


        self.peptide_m.load_peptide_fragments_annotation(self.experiment)
        end2 = time.time()
        logging.info(msg="[{}]load fragments: ".format(os.path.basename(self.experiment.peptide_identification_file))+time.strftime("%H:%M:%S", time.gmtime(end2 - end)))
        sys.stdout.flush()

        self.peptide_m.find_parent_protein_id_for_peptides(Param.DB_PATH)
        end3 = time.time()
        logging.info(msg="[{}]find parent proteins: ".format(os.path.basename(self.experiment.peptide_identification_file))+time.strftime("%H:%M:%S", time.gmtime(end3 - end2)))
        sys.stdout.flush()


        cdr3_peptides = self.assembly(Param.DB_PATH)
        end4 = time.time()
        logging.info(msg="[{}]assemble: ".format(os.path.basename(self.experiment.peptide_identification_file))+time.strftime("%H:%M:%S", time.gmtime(end4 - end3)))
        sys.stdout.flush()

        self.peptide_m.output_peptides_information()


        return cdr3_peptides





def runProcess(experiment:Experiemnt,param_file:str):

    p = Processor(experiment,param_file)

    peptides_to_quantify =  p.process()

    logging.info(msg="[{}] get {} CDR3 peptides: ".format(os.path.basename(experiment.peptide_identification_file), len(peptides_to_quantify)))

    tmp_filename = serialize_proteins(p.proteins)


    return (tmp_filename,peptides_to_quantify)

def serialize_proteins(proteins:list):
    with tempfile.NamedTemporaryFile(delete=False) as tmp_file:
        pickle.dump(proteins, open(tmp_file.name, 'wb'))
        tmp_filename = tmp_file.name
    logging.info("Serialized process")
    return tmp_filename

def deserialize_proteins(tmp_filename):
    _proteins = pickle.load(open(tmp_filename, 'rb'))
    os.remove(tmp_filename)
    logging.info("Deserialized process")

    return _proteins

def merge_proteins(serialized_proteins):
    final_proteins = dict()

    for tmp_file in serialized_proteins:
        _proteins = deserialize_proteins(tmp_file)
        for _protein in _proteins:
            if _protein.id in final_proteins.keys():
                '''
                Only keep the best covered protein
                overrided __lt__
                '''
                if compare_protein_coverage(final_proteins[_protein.id], _protein):
                    final_proteins[_protein.id] = _protein
            else:
                final_proteins[_protein.id] = _protein

    return final_proteins

def output_proteins_without_quantification(final_proteins):
    with open(Param.TASK_NAME+"_out_proteins_noquantification.csv","w+") as out:
        out.write("id,sequence,coverage,cdr coverage,")
        out.write("cdr1,cdr1 coverage,cdr1 peptides,")
        out.write("cdr2,cdr2 coverage,cdr2 peptides,")
        out.write("cdr3,cdr3 coverage,cdr3 peptides,")
        out.write("source\n")

        for _protein in final_proteins.values():
            out.write(_protein.simple_str())



logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')

def pipeline(param_file):

    read_params(param_file)

    start = time.time()

    '''
    Read experiment setting from config file
    '''
    experiment_m = ExperimentManager()
    CleavageManager.verify_index(Param.DB_PATH)
    '''
    Multiprocessing Protein aseembly
    '''

    pool = mp.Pool(processes=Param.PROCESS_NUM)

    _partial = partial(runProcess,param_file =param_file)
    results = pool.map(_partial,experiment_m.experiments)

    pool.close()
    pool.join()

    #results = [runProcess(experiment) for experiment in experiment_m.experiments]

    '''
    Quantify CDR3 peptides
    '''
    serialized_proteins = list()
    cdr_peptides_result = list()
    for result in results:
        serialized_proteins.append(result[0])
        cdr_peptides_result.append(result[1])


    end = time.time()
    logging.info(msg="Protein Assembly all end: " + time.strftime("%H:%M:%S", time.gmtime(end - start)))

    if Param.QUANTIFICATION == False:
        final_proteins_without_quantification = merge_proteins(serialized_proteins)
        output_proteins_without_quantification(final_proteins_without_quantification)
        return



    logging.info(msg="start quantification")
    quantification_m = QuantificationManager()
    quantification_m.register_experiment_manager(experiment_m)

    end2 = time.time()
    logging.info(msg="{} RT shift matrix:".format(str(quantification_m.em.rt_shift))+time.strftime("%H:%M:%S", time.gmtime(end2 - end)))
    for cdr_peptides in cdr_peptides_result:
        _cdr3_peptides = cdr_peptides

        for _peptide in _cdr3_peptides:
            quantification_m.add_peptide(peptide = _peptide.sequence,mz = _peptide.mz,charge = _peptide.charge,
                                         rt=_peptide.retention_time,index=_peptide.index,belonging=Belonging.CDR3)
            for _additional_peptide in _peptide.additional:
                quantification_m.add_peptide(peptide=_additional_peptide.sequence,mz = _additional_peptide.mz, charge=_additional_peptide.charge,
                                             rt=_additional_peptide.retention_time, index=_additional_peptide.index,belonging=Belonging.CDR3)


    quantification_m.load_raw_files()
    quantification_m.fill_rt_array_for_peptides()
    quantification_m.quantify_peptides()
    quantification_m.classify_peptides()
    quantification_m.classify_peptides_2()

    quantification_m.close_raw_files()
    quantification_m.output_peptide_quantification_table()

    '''
    Merge protein identifications from this batch of data
    Assign protein abudance based on CDR3 peptides
    Assign protein affinity group 
    '''

    final_proteins = dict()

    for tmp_file in serialized_proteins:
        _proteins = deserialize_proteins(tmp_file)
        for _protein in _proteins:
            if _protein.id in final_proteins.keys():
                '''
                Only keep the best covered protein
                overrided __lt__
                '''
                if compare_protein_coverage(final_proteins[_protein.id],_protein):
                    final_proteins[_protein.id] = _protein
            else:
                final_proteins[_protein.id] = _protein

    for _protein in final_proteins.values():
        quantification_m.generate_abundance_for_protein(_protein)
        quantification_m.classify_protein(_protein)
        quantification_m.classify_protein_2(_protein)

    with open(Param.TASK_NAME+"_out_proteins.csv","w+") as out:
        out.write("id,sequence,coverage,cdr coverage,")
        out.write("cdr1,cdr1 coverage,cdr1 peptides,")
        out.write("cdr2,cdr2 coverage,cdr2 peptides,")
        out.write("cdr3,cdr3 coverage,cdr3 peptides,")
        out.write("cdr3 abundance 1,cdr3 abundance 2,cdr3 abundance 3,cdr3 type,cdr3 type2,")
        out.write("source\n")

        for _protein in final_proteins.values():
            out.write(str(_protein))
        #
        #Classify protein affinity
        #

    '''
    Group proteins based on CDR1,CDR2,CDR3 and classified affinity
    '''





    logging.info(msg="end")

if __name__ == "__main__":
    paramas = ["tmp//"+f for f in os.listdir("tmp//")]
    print(paramas)
    for param in paramas:
        try:
            pipeline(param)
        except Exception as e:
            print(param)
            raise IOError("end")


