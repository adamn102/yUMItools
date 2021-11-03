# The data simulator will be used to generate test data for yUMItools testing

import numpy as np
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import random
from Bio import SeqIO, bgzf


class TestTube(object):

    def __init__(self, refernce_sequence_filepath):

        # assign the file path for the reference sequence
        self.refence_sequence_filepath = refernce_sequence_filepath

        # read in reference sequence
        self._read_fasta_file()

        # find barcodes in reference sequence
        self.find_barcodes(self.reference_sequence_dict[self.record_id])

        # list of library sequences
        self.barcode_library_list = list()
        self.sequence_pool = list()

    def _read_fasta_file(self):
        from Bio import SeqIO

        self.reference_sequence_dict = dict()

        for record in SeqIO.parse(self.refence_sequence_filepath, "fasta"):
            self.record_id = record.id
            self.reference_sequence_dict[record.id] = record.seq
            # self.reference_sequence = record

    def find_barcodes(self, library_sequence):
        """
        takes in a library reference sequence (with N-mers) and returns a dict
        containing the positions of the barcode
        :param library_sequence:
        :return:
        """

        # locate the barcodes within the sequence
        cache = 0
        barcode_dict = dict()
        tmp_umi_sequence = library_sequence + "X"
        i = 0
        barcode_end = 0
        while tmp_umi_sequence.find('N') != -1:

            barcode_start = tmp_umi_sequence.find('N')

            # find barcode length
            length = 0
            while tmp_umi_sequence[barcode_start + length] == "N":
                length += 1
            barcode_end = barcode_start + length

            barcode_dict['barcode' + str(i)] = (barcode_start + cache, barcode_end + cache)

            i += 1
            cache += barcode_end
            tmp_umi_sequence = tmp_umi_sequence[barcode_end:]

        self.barcode_dict = barcode_dict

    def generate_barcode_library(self, clones=10):
        """
        Take the reference sequence and return a realized barcode library
        generated a sequence set of n-clones
        :return:
        """

        barcode_library_list = []

        for i in range(clones):
            # resolve barcode
            umi_set = resolve_barcode(self.reference_sequence_dict[self.record_id],
                                      self.barcode_dict,
                                      self.record_id + "_" + str(i))
            # mutation umi_seq
            # umi_set = reverse_transcription_mutate(umi_set, mutation_rate)

            barcode_library_list.append(umi_set)

        self.barcode_library_list = barcode_library_list


def reverse_transcribe_library(barcode_library_list, clones=2000, mut_type='random', mutation_rate=1. / 10000):
    """

    :param barcode_library_list:
    :param clones:
    :param mutation_rate:
    :return:
    """
    # this a reverse transcribed subset of the library
    rt_library = []

    # random subset of library
    idx = np.random.randint(len(barcode_library_list), size=clones)
    umi_sets = np.array(barcode_library_list)[idx]

    for clone in umi_sets:
        if mut_type == 'random':
            mut_clone = reverse_transcription_mutate(clone, mutation_rate)
        elif mut_type == 'count':
            mut_clone = reverse_transcription_mutate_count(clone, mutation_rate)
        rt_library.append(mut_clone)

    return rt_library


def resolve_barcode(reference_sequence, barcode_dict, id_name):
    """
    Take one instance of a reference sequence and population the barcodes
    :param id_name:
    :param barcode_dict:
    :param reference_sequence:
    :return:
    """

    final_sequence = str(reference_sequence)
    for k, v in barcode_dict.items():
        umi_len = int(v[1] - v[0])
        final_sequence = final_sequence.replace(umi_len * "N", generate_barcode(umi_len), 1)

    # convert to SeqRecord format
    # record = SeqRecord(Seq(final_sequence), id_name, "")
    # return as a string
    record = final_sequence

    return record


def generate_barcode(n=15):
    """
    Create a n-mer seq object barcode

    To Do: Model barcode params from actual barcode data - have the fraction of each base match
    observed distributions
    """

    bases = ["A", "C", "T", "G"]
    barcode = ''.join(np.random.choice(bases, size=n, replace=True))

    return barcode


def write_fasta_sequence_list(sequence_list, outfile):
    # need to convert data to SeqRecord
    sequence_record_list = []
    count = 0

    for read in sequence_list:
        sequence_record_list.append(
            SeqRecord(Seq(read),
                      str(count), str(count)))
        count += 1

    with open(outfile, "w") as output_handle:
        SeqIO.write(sequence_record_list, output_handle, "fasta")


def reverse_transcription_mutate(umi_set, mutation_rate=1. / 10000):
    template = str(umi_set)

    mu = mutation_rate / 3
    mutation_matrix = np.array(
        [[1 - 3 * mu, mu, mu, mu],
         [mu, 1 - 3 * mu, mu, mu],
         [mu, mu, 1 - 3 * mu, mu],
         [mu, mu, mu, 1 - 3 * mu]]
    )

    # iterate though each base pair
    # need to change to a matrix based method

    mutation_list = []

    for i in range(len(template)):
        ref = str(template[i])
        mut = _mutation(ref, mutation_rate)
        if ref != mut:
            mutation_list.append((i, mut))

    # update sequence
    for mutation in mutation_list:
        template = _update_sequence(template, mutation[0], mutation[1])

    # make a SeqRecord
    # updated_template = SeqRecord(Seq(template), umi_set.id, umi_set.description)
    updated_template = template

    return updated_template


def reverse_transcription_mutate_count(umi_set, mutation_rate=1. / 10000):
    template = str(umi_set)

    mu = 1 / 3
    mutation_matrix = np.array(
        [[1 - 3 * mu, mu, mu, mu],
         [mu, 1 - 3 * mu, mu, mu],
         [mu, mu, 1 - 3 * mu, mu],
         [mu, mu, mu, 1 - 3 * mu]]
    )

    # iterate though each base pair
    # need to change to a matrix based method

    mutation_list = []

    count = 1
    for i in range(len(template)):
        if count % int(1 / mutation_rate) == 0:
            ref = str(template[i])
            mut = _mutation_count(ref)
            mutation_list.append((i + 1, mut))
        count += 1
    # update sequence
    for mutation in mutation_list:
        template = _update_sequence(template, mutation[0], mutation[1])

    # make a SeqRecord
    # updated_template = SeqRecord(Seq(template), umi_set.id, umi_set.description)
    updated_template = template

    return updated_template


def _mutation(ref, mutation_rate=1. / 10000):
    import numpy as np

    """
    Input: a specific basepair
    """
    index = dict({'A': 0, 'C': 1, 'T': 2, 'G': 3})
    mu = mutation_rate / 3
    mutation_matrix = np.array(
        [[1 - 3 * mu, mu, mu, mu],
         [mu, 1 - 3 * mu, mu, mu],
         [mu, mu, 1 - 3 * mu, mu],
         [mu, mu, mu, 1 - 3 * mu]]
    )

    # sample for mutation from mutation matrix
    mutation = np.array(['A', 'C', 'T', 'G'])[np.random.multinomial(1, mutation_matrix[index[ref]]) == 1][0]
    return mutation


def _mutation_count(ref):
    import numpy as np

    """
    Input: a specific basepair
    """
    index = dict({'A': 0, 'C': 1, 'T': 2, 'G': 3})
    mu = 1
    mutation_matrix = np.array(
        [[0, mu / 3, mu / 3, mu / 3],
         [mu / 3, 0, mu / 3, mu / 3],
         [mu / 3, mu / 3, 0, mu / 3],
         [mu / 3, mu / 3, mu / 3, 0]]
    )
    # sample for mutation from mutation matrix
    mutation = np.array(['A', 'C', 'T', 'G'])[np.random.multinomial(1, mutation_matrix[index[ref]]) == 1][0]
    return mutation


def _update_sequence(sequence, site, mutation):
    """
    for a specific site and mutation, update sequence
    """
    return (sequence[:site - 1] + mutation + sequence[site:])


def pcr_amplificaiton(template, cycles=30, p=0.5):
    # convert to np array
    # temps = np.array(template)
    temps = np.array([list(seq) for seq in template])

    # cycle
    cycle_number = 0
    while cycle_number < cycles:
        # select random fragments to amplify
        idx = np.random.randint(len(temps), size=round(len(temps) * p))

        # index template pool from random fragments
        duplicate = temps[idx]

        # mutate with Q5
        duplicate = pcr_mutate(duplicate, mut_rate=1. / 100000)

        # update total pool
        temps = np.append(temps, duplicate, axis=0)

        # update cycle number
        cycle_number += 1
    # convert into list of strigs
    output_list = ["".join(seq) for seq in temps]

    return output_list


def library_amp(template, cycles=30, p=0.5):
    """
    Take a barcoded library and amplifies it
    :param template:
    :param cycles:
    :param p:
    :return:
    """

    temps = np.array(template)

    # cycle
    cycle_number = 0
    while cycle_number < cycles:
        # select random fragments to amplify
        idx = np.random.randint(len(temps), size=round(len(temps) * p))

        # index template pool from random fragments
        duplicate = temps[idx]

        # update total pool
        temps = np.append(temps, duplicate, axis=0)

        # update cycle number
        cycle_number += 1

    # # convert into list of strigs
    # output_list = ["".join(seq) for seq in temps]

    return temps


def pcr_mutate(duplicate, mut_rate=1. / 1000):
    """
    PCR mutation function
    :param duplicate: s
    :param mut_rate:
    :return:
    """

    duplicate = np.array([list(seq) for seq in duplicate])

    # shape data_len
    shape = duplicate.shape[1]
    # flatten data
    flat_data = duplicate.flatten()

    # select nts
    data_len = flat_data.shape[0]
    indx = np.random.randint(low=0, high=data_len, size=round(data_len * mut_rate))

    for i in range(len([indx][0])):
        flat_data[[indx][0][i]] = q5_mutate(flat_data[[indx][0][i]])
    return flat_data.reshape(-1, shape)


def q5_mutate(ref):
    """
    Mutation function for the simulated PCR reaction

    Input: A nt that has been selected to mutate

    Output: Returns a nt based on the mutation mutation_matrix

    Notes:
    currently, the prob of all mutations are equal (1/3)
    this should be updated in the future.
    """
    mutation_matrix = np.array(
        [[0, 1. / 3, 1. / 3, 1. / 3],
         [1. / 3, 0, 1. / 3, 1. / 3],
         [1. / 3, 1. / 3, 0, 1. / 3],
         [1. / 3, 1. / 3, 1. / 3, 0]]
    )
    index = dict({'A': 0, 'C': 1, 'T': 2, 'G': 3})
    return np.array(['A', 'C', 'T', 'G'])[np.random.multinomial(1, mutation_matrix[index[ref]]) == 1][0]


def random_str_index(string, length):
    # randomly subset a string
    index = random.randrange(0, len(string) - length + 1)
    return string[index: (index + length)]


def tagment(template_pool, ave_size=500, std=50):
    """
    the tagment function will take in a numpy array of sequence data
    and prune the fragments to random sizes - following a normal distribution
    in order to simulate the tagmentation process

    todo - figure out how to index each row with a different slice size
    maybe need to convert into a string first
    :param std:
    :param template_pool:
    :param ave_size:
    :return:
    """
    output_library = []

    for sequence in template_pool:
        # select a random length
        rand_length = int(np.round(np.random.normal(loc=ave_size, scale=std, size=1)))

        # random subset
        output_library.append(random_str_index(sequence, rand_length))

    return output_library


def deep_sequence(template_pool, output_file, read_length=300, coverage=10000):
    # select a subset of templates to sequence
    indx = np.random.randint(low=0, high=len(template_pool), size=coverage)
    selected_reads = np.array(template_pool)[indx]

    # output lists
    R1_reads = []
    R2_reads = []

    count = 0
    for read in selected_reads:
        if len(read) >= read_length:
            R1_name = _generate_description(count, '1')
            R1 = SeqRecord(Seq(''.join(read[:read_length])),
                           id=R1_name,
                           name=R1_name,
                           description='')

            # add qual info
            R1.letter_annotations['phred_quality'] = _qual_score(R1)
            R1_reads.append(R1)

            R2_name = _generate_description(count, '2')
            R2 = SeqRecord(Seq(''.join(read[-read_length:])),
                           id=R2_name,
                           name=R2_name,
                           description='')
            R2.letter_annotations['phred_quality'] = _qual_score(R2)
            R2.seq = R2.seq.reverse_complement()
            R2_reads.append(R2)
        count += 1

    # write R1 reads
    with bgzf.BgzfWriter(output_file + "_S00_L001_R1_001.fastq.gz", "wb") as outgz:
        SeqIO.write(sequences=R1_reads, handle=outgz, format="fastq")

    # write R2 reads
    with bgzf.BgzfWriter(output_file + "_S00_L001_R2_001.fastq.gz", "wb") as outgz:
        SeqIO.write(sequences=R2_reads, handle=outgz, format="fastq")


def hamming(s1, s2):
    if isinstance(s1, str) == 0:
        s1 = str(s1)
    if isinstance(s2, str) == 0:
        s2 = str(s2)

    s1 = list(s1)
    s2 = list(s2)

    dist = len([i for i, j in zip(s1, s2) if i != j and i != 'N' and j != 'N'])

    return dist


def _qual_score(seqRecord_file, quality=38):
    qual_lengeth = len(seqRecord_file)
    phred_quality = []

    for i in range(qual_lengeth):
        phred_quality.append(quality)

    return phred_quality


def _generate_name(count):
    xcord = str(count)
    ycord = "0000"
    name = ':'.join(['ERRSEQSIM', "14", "000000000-C2CG5", "1", "1101", xcord, ycord])
    return name


def _generate_description(count, read):
    filter = 'N'
    control_number = '0'
    sample_sheet_number = '11'

    name_1 = _generate_name(count)
    name_2 = ':'.join([str(read), filter, control_number, sample_sheet_number])

    description = ' '.join([name_1, name_2])
    return description
