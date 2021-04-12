import pysam
from collections import defaultdict, Counter
import numpy as np
import pandas as pd

chars = 'ACGTN'
char_to_int = dict((c, i) for i, c in enumerate(chars))
int_to_char = dict((i, c) for i, c in enumerate(chars))


class YUMISet:
    """
    yUMI set Data Class

    This is the primary data structure for the yUMItools software

    Main functionality will include:

    * Parsing multiple bam data sets for any one bed formatted reference data set
    * error correcting UMIs
    * phasing multiple UMIs
    * generating consensus sequences
    * comparing multiple yUMI sets to compute mutation data
    * generate quality metrics
    """

    def __init__(self, reference_sequence, bamfile):
        self.reference_sequence = reference_sequence
        self.reference_sequence_name = self.reference_sequence.split("/")[-1].split(".")[0]
        self.barcode_dict = dict()
        self.flankseq_dict = dict()

        # parse reference sequence - extract barcodes and barcode locations
        self.parse_reference_sequence()

        # read in the samfile
        self.samfile = pysam.AlignmentFile(bamfile, "rb")

        # read dict for sam file
        self.read_dict = defaultdict(lambda: [None, None])
        self.read_pair_dict()

    def parse_reference_sequence(self):
        """
        This function will parse reference sequences for UMI locations
        and reference sequences

        :param reference_sequence:
        :return:
        * dict of UMI locations
        * UMI flanking sequences
        * reference sequence
        """
        from Bio import SeqIO
        cache = 0
        barcode_dict = dict()
        flankseq_dict = dict()

        for record in SeqIO.parse(self.reference_sequence, "fasta"):
            # add the X to prevent index overflow
            # this is for barcodes at the end of the sequence
            sequence = record.seq + "X"
            i = 0
            barcode_end = 0
            while sequence.find('N') != -1:

                barcode_start = sequence.find('N')

                # find barcode length
                length = 0
                while sequence[barcode_start + length] == "N":
                    length += 1
                barcode_end = barcode_start + length

                barcode_dict['barcode' + str(i)] = (barcode_start + cache, barcode_end + cache)

                # find the 15bp flanking sequences of the barcodes
                # will not work properly if barcode distance is < 15 bp
                five_prime_flank = str(record.seq[barcode_start + cache - 15: barcode_start + cache])
                three_prime_flank = str(record.seq[barcode_end + cache: barcode_end + cache + 15])
                flankseq_dict['barcode' + str(i)] = (five_prime_flank, three_prime_flank)

                i += 1
                cache += barcode_end
                sequence = sequence[barcode_end:]

        # update object variables
        self.barcode_dict = barcode_dict
        self.flankseq_dict = flankseq_dict

    def parse_bam_file(self, bam_file):
        """
        create a barcode dict where the keys are the barcode
        and the values are the ID values for reads that map
        to the key barcode

        for paired end reads - if one read maps to the barcode, make sure to keep both reads
        :param bam_file:
        :return:
        """
        pass

    def read_pair_dict(self):
        # returns a dict with read names as keys and read files as values
        # may get kinda large - be sure to del when done

        for read in self.samfile.fetch(None):
            if read.is_paired:
                read_name = read.query_name

            if read_name not in self.read_dict:
                if read.is_read1:
                    self.read_dict[read_name][0] = read
                else:
                    self.read_dict[read_name][1] = read
            else:
                if read.is_read1:
                    self.read_dict[read_name][0] = read
                else:
                    self.read_dict[read_name][1] = read

    def umi_phase(self, barcode_sequence, position_A, position_B, distance=2):

        # take the reads mapped to the barcode_sequence at position_A and
        # identify the phased barcode at position_B

        # fetch reads to barcode in position A
        umi_data = fetch_barcode_reads(barcode=barcode_sequence,
                                       samfile=self.samfile,
                                       reference=self.reference_sequence.split("/")[-1].split(".")[0],
                                       barcode_location=self.barcode_dict['barcode' + str(position_A)],
                                       barcode_flank=self.flankseq_dict['barcode' + str(position_A)],
                                       distance=distance)

        # find reads from position A that overlap position B
        umi_list = barcode_extraction_read_list(umi_data,
                                                self.flankseq_dict['barcode' + str(position_B)],
                                                self.read_dict)

        # return the most common barcode
        # todo: return umi metrics
        return Counter(umi_list).most_common(1)[0][0]

    def stitching(self, seed):

        # iterate the phasing of multiple barcodes together from
        # barcode 0 to barcode 5

        umi_set = defaultdict(lambda: [None, None, None, None, None, None])
        x = 0
        umi_set[x] = seed
        while x < 5:
            seed = self.umi_phase(seed, x, x + 1)
            x = x + 1
            umi_set[x] = seed
        return umi_set


def barcode_extraction(samfile, reference, barcode_location, barcode_flank):
    # list for collecting barcodes
    barcodes = []

    bam_iter = samfile.fetch(reference, barcode_location[0], barcode_location[1])
    for x in bam_iter:
        if barcode_flank[0][10:] in x.seq:
            if barcode_flank[1][:5] in x.seq:
                start = x.seq.find(barcode_flank[0][10:]) + 5
                end = x.seq.find(barcode_flank[1][:5])
                if start < end and len(x.seq[start:end]) == 15:
                    barcodes.append(x.seq[start:end])
    return barcodes


def barcode_extraction_dict(samfile, reference, barcode_location, barcode_flank):
    # similar to barcode extraction, but returns a dict of barcodes with reads as the values
    # that can be used to parse data later on

    # list for collecting barcodes
    barcodes_dict = dict()

    bam_iter = samfile.fetch(reference, barcode_location[0], barcode_location[1])
    input_reads = 0
    total_reads = 0
    for x in bam_iter:
        input_reads += 1
        if (barcode_flank[0][12:] in x.seq) & (barcode_flank[1][:3] in x.seq):
            if (barcode_location[0] in x.positions) & (barcode_location[1] in x.positions):
                # total_reads += 1

                start = x.positions.index(barcode_location[0])
                end = x.positions.index(barcode_location[1])

                if start < end and len(x.seq[start:end]) == 15:
                    total_reads += 1
                    barcode = x.seq[start:end]

                    if barcode not in barcodes_dict.keys():
                        barcodes_dict[barcode] = [x.query_name]
                    else:
                        barcodes_dict[barcode].append(x.query_name)
    return barcodes_dict


def error_correct_barcode_dict(barcodes_dict, method='hard_limit'):
    if method == "hard_limit":
        error_corrected_barcode_dict = correct_barcodes_cutoff(barcodes_dict, cutoff=10)

    if method == "grouping":
        error_corrected_barcode_dict = correct_barcodes_grouping(barcodes_dict, cutoff=10)

    return (error_corrected_barcode_dict)


def correct_barcodes_cutoff(barcodes_dict, cutoff=10):
    """
    fillters keys,values from barcode dictionary where the number of reads 
    in the key is less than the cutoff
    
    :param barcodes_dict: 
    :param cutoff: 
    :return: 
    """
    # output dict
    error_corrected_barcode_dict = dict()

    for k, v in barcodes_dict.items():
        if len(v) >= cutoff:
            error_corrected_barcode_dict[k] = v

    return (error_corrected_barcode_dict)


def correct_barcodes_grouping(barcodes_dict, cutoff=10):
    pass


def fetch_reads_from_dict(barcodes_dict, barcode):
    reads = barcodes_dict[barcode]

    return reads


def barcode_extraction_read_list(read_list, barcode_flank, read_dict):
    # list for collecting barcodes
    barcodes = []

    read_data_list = []
    [read_data_list.append(read_dict[x][0]) for x in read_list]
    [read_data_list.append(read_dict[x][1]) for x in read_list]

    for x in read_data_list:
        if barcode_flank[0][10:] in x.seq:
            if barcode_flank[1][:5] in x.seq:
                start = x.seq.find(barcode_flank[0][10:]) + 5
                end = x.seq.find(barcode_flank[1][:5])
                if start < end and len(x.seq[start:end]) == 15:
                    barcodes.append(x.seq[start:end])
    return barcodes


def fetch_barcode_reads(barcode, samfile, reference, barcode_location, barcode_flank, distance=2):
    # goal of this function is to fetch paired reads to specific barcode regions
    # for any specific barcode will fetch the paired reads

    # list for collecting barcodes
    barcode_dict = dict()
    read_list = []

    bam_iter = samfile.fetch(reference, barcode_location[0], barcode_location[1])
    for x in bam_iter:
        if barcode_flank[0][10:] in x.seq:
            if barcode_flank[1][:5] in x.seq:
                start = x.seq.find(barcode_flank[0][10:]) + 5
                end = x.seq.find(barcode_flank[1][:5])
                if start < end and len(x.seq[start:end]) == 15:
                    umi = x.seq[start:end]
                    if hamming(umi, barcode) <= distance:
                        if x.query_name not in read_list:
                            read_list.append(x.query_name)
                        # if umi not in barcode_dict.keys():
                        #     x.query_name
                        #     barcode_dict[barcode] = [x.query_name]
                        #
                        # else:
                        #     barcode_dict[barcode].append(x.query_name)
    return read_list


def write_read_list(read_list, read_dict, output_file, samfile):
    pairedreads = pysam.AlignmentFile(output_file, "w", template=samfile)
    [pairedreads.write(read_dict[x][0]) for x in read_list]
    [pairedreads.write(read_dict[x][1]) for x in read_list]
    pairedreads.close()


def hamming(s1, s2):
    if isinstance(s1, str) == 0:
        s1 = str(s1)
    if isinstance(s2, str) == 0:
        s2 = str(s2)

    s1 = list(s1)
    s2 = list(s2)

    dist = len([i for i, j in zip(s1, s2) if i != j and i != 'N' and j != 'N'])

    return dist


def sequences_to_one_hot(sequences, chars='ACGTN'):
    """

    :param sequences:
    :param chars:
    :return:
    """
    seqlen = len(sequences[0])
    char_to_int = dict((c, i) for i, c in enumerate(chars))

    one_hot_encoded = []
    for seq in sequences:

        onehot_seq = []

        integer_encoded = [char_to_int[base] for base in seq]
        for value in integer_encoded:
            letter = [0 for _ in range(len(chars))]
            letter[value] = 1
            onehot_seq.append(letter)
        one_hot_encoded.append(onehot_seq)

    one_hot_encoded = np.array(one_hot_encoded)

    return one_hot_encoded


def consensus_caller(split_arr, min_coverage=5):
    sequence_list = []
    position_list = []
    coverage_list = []
    fraction_list = []

    for pos in range(len(split_arr)):
        coverage = len(split_arr[pos])
        if coverage > min_coverage:
            onehot_group = split_arr[pos][:, 0:5]
            base_averages = np.average(onehot_group, axis=0)
            top_base = np.argmax(base_averages, axis=0)
            top_base_acc = np.amax(base_averages, axis=0)

            consensus_bases = int_to_char[top_base]

            position_list.append(split_arr[pos][:, 5][0])
            sequence_list.append(consensus_bases)
            coverage_list.append(coverage)
            fraction_list.append(top_base_acc)

    return position_list, sequence_list, coverage_list, fraction_list


def UMI_consensus(barcode, corrected_barcode_dict, y_file):
    barcode = barcode
    sequences = corrected_barcode_dict[str(barcode)]

    # parse reads and bases
    positions = []
    bases = []
    quality = []
    for read_name in sequences:
        for read in y_file.read_dict[read_name]:
            # make sure the read have the same length for position and sequence
            if len(read.positions) == len(read.seq):
                positions.extend(read.positions)
                bases.extend(read.seq)
                quality.extend(read.qual)
    if len(positions) == 0:
        return 0

    # one hot encode
    integer_encoded = [char_to_int[base] for base in bases]

    onehot_seq = []
    for value in integer_encoded:
        letter = [0 for _ in range(len(chars))]
        letter[value] = 1
        onehot_seq.append(letter)

    one_hot_encoded = np.array(onehot_seq)
    # print(len(one_hot_encoded))
    positions_array = np.array(positions)
    # print(len(positions_array))
    # positions_list = set(positions)

    # for pos in positions_list:
    #     site_data = one_hot_encoded[int(pos) == positions_array]

    # sorted_positions_array = positions_array.argsort()
    # sorted_one_hot_encoded = one_hot_encoded[sorted_positions_array]
    # sorted_positions_array = sorted_positions_array.reshape((len(sorted_positions_array),1))

    # merge arrays
    positions_array = positions_array.reshape((len(positions_array), 1))
    arr = np.hstack((one_hot_encoded, positions_array))

    # sort array on position
    arr = arr[arr[:, 5].argsort()]

    # split array by position
    split_arr = np.array_split(arr, np.where(np.diff(arr[:, 5]))[0] + 1)

    position_list, sequence_list, coverage_list, fraction_list = consensus_caller(split_arr)

    df = pd.DataFrame({
        'position': position_list,
        'base': sequence_list,
        'coverage': coverage_list,
        'fraction': fraction_list
    })
    df['UMI'] = barcode
    return df


def consensus_df(barcode_dict, y_file, min_coverage=5):
    count = 0
    for umi in barcode_dict.keys():
        # print(umi)
        if len(barcode_dict[umi]) >= min_coverage:
            if count == 0:
                df = UMI_consensus(umi, barcode_dict, y_file=y_file)
                if type(df) != int:
                    count += 1
            else:
                df2 = UMI_consensus(umi, barcode_dict, y_file=y_file)
                if type(df2) != int:
                    df = df.append(df2)

    # sample_df = df
    df['position'] = df['position'].astype(int)
    df['UMI_pos'] = df['UMI'] + "_" + df['position'].astype(str)
    return df


def pipeline_function(y_lib, y_sample, cutoff=5):
    barcodes = list(y_lib.barcode_dict.keys())

    count = 0

    for barcode in barcodes:

        # select one barcode location
        lib_barcodes_dict = barcode_extraction_dict(samfile=y_lib.samfile,
                                                    reference='Rp0-reference',
                                                    barcode_location=y_lib.barcode_dict[barcode],
                                                    barcode_flank=y_lib.flankseq_dict[barcode])

        sample_barcodes_dict = barcode_extraction_dict(samfile=y_sample.samfile,
                                                       reference='Rp0-reference',
                                                       barcode_location=y_sample.barcode_dict[barcode],
                                                       barcode_flank=y_sample.flankseq_dict[barcode])

        lib_corrected_barcode_dict = correct_barcodes_cutoff(lib_barcodes_dict, cutoff=cutoff)
        print(barcode, "library barcodes:", len(lib_corrected_barcode_dict.keys()))

        sample_corrected_barcode_dict = correct_barcodes_cutoff(sample_barcodes_dict, cutoff=cutoff)
        print(barcode, 'sample barcodes:', len(sample_corrected_barcode_dict.keys()))

        # overlap between the two libraries
        overlap = len(list(set(sample_corrected_barcode_dict.keys()) & set(lib_corrected_barcode_dict.keys())))
        print(barcode, "sample and library overlap", overlap)

        lib_df = consensus_df(lib_corrected_barcode_dict, y_file=y_lib,min_coverage=cutoff)
        samp_df = consensus_df(sample_corrected_barcode_dict, y_file=y_sample,min_coverage=cutoff)

        if count == 0:
            # merge sample and lib df
            merged_df = pd.merge(samp_df, lib_df, how="inner", on='UMI_pos')
        else:
            tmp_df = pd.merge(samp_df, lib_df, how="inner", on='UMI_pos')
            merged_df = merged_df.append(tmp_df)
        count += 1

    total_sites = merged_df.shape[0]
    total_mutations = merged_df[merged_df['base_x'] != merged_df['base_y']].shape[0]

    print(total_sites, total_mutations)

    return merged_df