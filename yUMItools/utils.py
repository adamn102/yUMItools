import pysam
from collections import defaultdict, Counter
import numpy as np
import pandas as pd
from collections import Counter

chars = 'ACGTN'
char_to_int = dict((c, i) for i, c in enumerate(chars))
int_to_char = dict((i, c) for i, c in enumerate(chars))


def char_to_int_multiple(string):
    output = ''
    for char in range(len(string)):
        output += str(char_to_int[string[char]])
    return output


def int_to_char_multiple(int_string):
    output = ''
    for char in range(len(int_string)):
        output += str(int_to_char[int(int_string[char])])
    return output


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
        self.reference = self.reference_sequence.split("/")[-1].split(".")[0]
        self.barcode_dict = dict()
        self.flankseq_dict = dict()

        # parse reference sequence - extract barcodes and barcode locations
        self.parse_reference_sequence()

        # read in the samfile
        self.samfile = pysam.AlignmentFile(bamfile, "rb")

        # read dict for sam file
        self.read_dict = defaultdict(lambda: [None, None])
        self.read_pair_dict()

        # create variable name for output file
        self.umi_df = pd.DataFrame()

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

    def barcode_extraction_dict(self, barcode):
        # similar to barcode extraction, but returns a dict of barcodes with reads as the values
        # that can be used to parse data later on
        # assumes all barcode are 15 mers - this should be update to be more dynamic

        # dict for collecting barcodes
        barcodes_dict = dict()

        # load barcode location
        barcode_location = self.barcode_dict[barcode]
        barcode_flank = self.flankseq_dict[barcode]

        bam_iter = self.samfile.fetch(self.reference, barcode_location[0], barcode_location[1])
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

    def library_pipeline(self, cutoff=5):
        count = 0
        barcodes = list(self.barcode_dict.keys())

        for barcode in barcodes:

            lib_barcodes_dict = self.barcode_extraction_dict(barcode)

            corrected_barcode_dict = correct_barcodes_cutoff(lib_barcodes_dict, cutoff=cutoff)
            print(barcode, "library barcodes:", len(corrected_barcode_dict.keys()))

            umi_df = self.consensus_df(corrected_barcode_dict, min_coverage=cutoff)

            umi_df['barcode'] = barcode

            if self.umi_df.empty == True:
                # merge df across multiple
                self.umi_df = umi_df
            else:
                self.umi_df = self.umi_df.append(umi_df)

        return self.umi_df

    def consensus_df(self, barcode_dict, min_coverage=5):
        count = 0
        for umi in barcode_dict.keys():
            # print(umi)
            if len(barcode_dict[umi]) >= min_coverage:
                if count == 0:
                    df = self.UMI_consensus(umi, barcode_dict)
                    if type(df) != int:
                        count += 1
                else:
                    df2 = self.UMI_consensus(umi, barcode_dict)
                    if type(df2) != int:
                        df = df.append(df2)

        # sample_df = df
        df['position'] = df['position'].astype(int)
        df['UMI_pos'] = df['UMI'] + "_" + df['position'].astype(str)
        return df

    def UMI_consensus(self, barcode, corrected_barcode_dict):

        sequences = corrected_barcode_dict[str(barcode)]

        # parse reads and bases
        positions = []
        bases = []
        quality = []
        for read_name in sequences:
            for read in self.read_dict[read_name]:
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
        positions_array = np.array(positions)

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

    def UMI_consensus_indel(self, barcode, corrected_barcode_dict, min_coverage=5):

        # find the reads (sequences) that correspond to a specific barcode
        read_list = corrected_barcode_dict[str(barcode)]
        #print(read_list)
        outer_count = 0
        for read_name in read_list:
            print(read_name)
            for read in self.read_dict[read_name]:
                # print(len(read.positions), len(read.seq))
                align_array = np.array(read.aligned_pairs)
                # print(align_array)
                # find deletions
                deletions = list(np.where(align_array[:, 0] == None)[0])
                insertions = list(np.where(align_array[:, 1] == None)[0])

                # convert bases and quals into an array
                quals = np.array(list(read.qual), dtype=object)
                bases = np.array(list(read.seq), dtype=object)
                # quals = read.qual
                # bases = read.seq

                #account for deletions
                for deletion in deletions:
                    print(deletion)
                    quals = np.append(np.append(quals[0:deletion], np.array(["N"], dtype=object),0),quals[deletion:],0)
                    bases = np.append(np.append(bases[0:deletion], np.array(["N"], dtype=object), 0),
                                      bases[deletion:], 0)


                # account for insertions

                insertion_lengths = list()
                insertion_positions = list()
                #print(insertions)
                for insertion in insertions:
                    count = 0
                    insertion_length = 0
                    align_array[insertion, ][1] = int(align_array[insertion + 1, ][1])
                    # check if insertions are sequential
                    if count < len(insertions)-1:
                        if insertion + 1 == insertions[count + 1]:
                            insertion_length += 1
                            count += 1
                    else:
                        insertion_lengths.append(insertion_length)
                        insertion_positions.append(insertion)

                    # insert new bases ans q-scores
                for i in range(len(insertion_positions)):
                    print(insertion_positions[i])
                    insert_bases = "".join(bases[insertion_positions[i]:insertion_positions[i] + insertion_lengths[i] + 1])
                    bases[insertion_positions[i]] = insert_bases

                    insert_quals = "".join(
                        quals[insertion_positions[i]:insertion_positions[i] + insertion_lengths[i] + 1])
                    quals[insertion_positions[i]] = insert_quals

                # print(bases)
                # remove old bases
                deletion_sites = [i + 1 for i in insertion_positions]
                bases = np.delete(bases, deletion_sites)
                quals = np.delete(quals, deletion_sites)
                align_array = np.delete(align_array, deletion_sites, axis=0)

                # print(bases)

                # qual_array = np.array([x for x in quals]).reshape(len(quals), 1)
                quals = quals.reshape(len(quals), 1)
                # print(qual_array.shape)
                # integer_encoded = np.array([int(char_to_int[base]) for base in bases]).reshape(len(bases), 1)
                bases = bases.reshape(len(bases), 1)
                # print(integer_encoded.shape)
                if outer_count == 0:
                    arr = np.concatenate((align_array, bases, quals), axis=1)
                    outer_count += 1
                else:
                    tmp_arr = np.concatenate((align_array, bases, quals), axis=1)
                    arr = np.concatenate((arr, tmp_arr), axis=0)

        # sort array on position
        arr[:, 1] = arr[:, 1].astype(int)
        arr = arr[arr[:, 1].argsort()]
        print(arr)
       #arr[:, 2] = arr[:, 2].astype(int)

        # np.unique(onehot_group, return_counts=True)

        # split array by position
        # find unique positions
        sequence_list = []
        position_list = []
        coverage_list = []
        fraction_list = []

        unique_positions, position_index = np.unique(arr[:, 1], return_index=True)
        split_arr = np.split(arr, position_index[1:])

        for position in range(len(split_arr)):

            coverage = len(split_arr[position])

            if coverage >= min_coverage:
                onehot_group = split_arr[position][:, 0:4]
                unique, counts = np.unique(onehot_group[:, 2], return_counts=True)

                # base_count = np.bincount(onehot_group[:, 2].astype(int))
                max_base = unique[np.argmax(counts)]
                top_base_acc = counts[np.argmax(counts)]
                #consensus_bases = int_to_char[max_base]

                position_list.append(split_arr[position][:, 1][0])
                sequence_list.append(max_base)
                coverage_list.append(coverage)
                fraction_list.append(top_base_acc)

        return position_list, sequence_list, coverage_list, fraction_list

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

    def library_pipeline_v2(self, cutoff=5):

        positions = []
        bases = []
        quality = []
        umi_list = []
        barcode_list = []

        for barcode in list(self.barcode_dict.keys()):
            print(barcode)
            lib_barcodes_dict = self.barcode_extraction_dict(barcode)
            umi_input_list = list(lib_barcodes_dict.keys())

            for umi in umi_input_list:
                sequences = lib_barcodes_dict[umi]
                if len(sequences) >= cutoff:
                    for read_name in sequences:
                        for read in self.read_dict[read_name]:
                            # make sure the read have the same length for position and sequence
                            if len(read.positions) == len(read.seq):
                                positions.extend(read.positions)
                                bases.extend(read.seq)
                                quality.extend(read.qual)
                                umi_list.extend([umi] * len(read.positions))
                                barcode_list.extend([barcode] * len(read.positions))

        df = pd.DataFrame({
            'position': positions,
            'base': bases,
            'coverage': quality,
            'UMI': umi_list,
            'barcode': barcode_list
        })

        group_df = df.groupby(['UMI', 'position', 'barcode'])['base'].apply(max_count)
        group_df = group_df.reset_index()
        group_df = group_df.fillna(0)

        group_df[['base', 'fraction', 'coverage']] = pd.DataFrame(group_df['base'].tolist(), index=group_df.index)
        group_df['UMI_pos'] = group_df['UMI'] + "_" + group_df['position'].astype(str)

        self.umi_df = group_df

        return self.umi_df

    def write_umi_data_to_bam(barcode, umi, output_file):
        """
        Write out the mapped bam data for any specific UMI -
        This is useful for visually inspecting mutations in a bam file
        :param umi:
        :param output_file:
        :return:
        """
        lib_barcodes_dict = self.barcode_extraction_dict(barcode)
        write_read_list(lib_barcodes_dict[umi], self.read_dict, output_file + umi + '.sam', self.samfile)


def max_count(string):
    out = Counter(string)
    value = out.most_common(1)[0][0]
    freq = out[value] / sum(out.values())
    out.values()
    return value, freq, out[value]
    # def barcode_extraction(samfile, reference, barcode_location, barcode_flank):


#     # list for collecting barcodes
#     barcodes = []
#
#     bam_iter = samfile.fetch(reference, barcode_location[0], barcode_location[1])
#     for x in bam_iter:
#         if barcode_flank[0][10:] in x.seq:
#             if barcode_flank[1][:5] in x.seq:
#                 start = x.seq.find(barcode_flank[0][10:]) + 5
#                 end = x.seq.find(barcode_flank[1][:5])
#                 if start < end and len(x.seq[start:end]) == 15:
#                     barcodes.append(x.seq[start:end])
#     return barcodes


# def barcode_extraction_dict(samfile, reference, barcode_location, barcode_flank):
#     # similar to barcode extraction, but returns a dict of barcodes with reads as the values
#     # that can be used to parse data later on
#
#     # list for collecting barcodes
#     barcodes_dict = dict()
#
#     bam_iter = samfile.fetch(reference, barcode_location[0], barcode_location[1])
#     input_reads = 0
#     total_reads = 0
#     for x in bam_iter:
#         input_reads += 1
#         if (barcode_flank[0][12:] in x.seq) & (barcode_flank[1][:3] in x.seq):
#             if (barcode_location[0] in x.positions) & (barcode_location[1] in x.positions):
#                 # total_reads += 1
#
#                 start = x.positions.index(barcode_location[0])
#                 end = x.positions.index(barcode_location[1])
#
#                 if start < end and len(x.seq[start:end]) == 15:
#                     total_reads += 1
#                     barcode = x.seq[start:end]
#
#                     if barcode not in barcodes_dict.keys():
#                         barcodes_dict[barcode] = [x.query_name]
#                     else:
#                         barcodes_dict[barcode].append(x.query_name)
#     return barcodes_dict

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

# def UMI_consensus(barcode, corrected_barcode_dict, y_file):
#     barcode = barcode
#     sequences = corrected_barcode_dict[str(barcode)]
#
#     # parse reads and bases
#     positions = []
#     bases = []
#     quality = []
#     for read_name in sequences:
#         for read in y_file.read_dict[read_name]:
#             # make sure the read have the same length for position and sequence
#             if len(read.positions) == len(read.seq):
#                 positions.extend(read.positions)
#                 bases.extend(read.seq)
#                 quality.extend(read.qual)
#     if len(positions) == 0:
#         return 0
#
#     # one hot encode
#     integer_encoded = [char_to_int[base] for base in bases]
#
#     onehot_seq = []
#     for value in integer_encoded:
#         letter = [0 for _ in range(len(chars))]
#         letter[value] = 1
#         onehot_seq.append(letter)
#
#     one_hot_encoded = np.array(onehot_seq)
#     # print(len(one_hot_encoded))
#     positions_array = np.array(positions)
#     # print(len(positions_array))
#     # positions_list = set(positions)
#
#     # for pos in positions_list:
#     #     site_data = one_hot_encoded[int(pos) == positions_array]
#
#     # sorted_positions_array = positions_array.argsort()
#     # sorted_one_hot_encoded = one_hot_encoded[sorted_positions_array]
#     # sorted_positions_array = sorted_positions_array.reshape((len(sorted_positions_array),1))
#
#     # merge arrays
#     positions_array = positions_array.reshape((len(positions_array), 1))
#     arr = np.hstack((one_hot_encoded, positions_array))
#
#     # sort array on position
#     arr = arr[arr[:, 5].argsort()]
#
#     # split array by position
#     split_arr = np.array_split(arr, np.where(np.diff(arr[:, 5]))[0] + 1)
#
#     position_list, sequence_list, coverage_list, fraction_list = consensus_caller(split_arr)
#
#     df = pd.DataFrame({
#         'position': position_list,
#         'base': sequence_list,
#         'coverage': coverage_list,
#         'fraction': fraction_list
#     })
#     df['UMI'] = barcode
#     return df

#
# def consensus_df(barcode_dict, y_file, min_coverage=5):
#     count = 0
#     for umi in barcode_dict.keys():
#         # print(umi)
#         if len(barcode_dict[umi]) >= min_coverage:
#             if count == 0:
#                 df = UMI_consensus(umi, barcode_dict, y_file=y_file)
#                 if type(df) != int:
#                     count += 1
#             else:
#                 df2 = UMI_consensus(umi, barcode_dict, y_file=y_file)
#                 if type(df2) != int:
#                     df = df.append(df2)
#
#     # sample_df = df
#     df['position'] = df['position'].astype(int)
#     df['UMI_pos'] = df['UMI'] + "_" + df['position'].astype(str)
#     return df


# def library_pipeline(y_lib, reference="Rp0-reference", cutoff=5):
#     count = 0
#     barcodes = list(y_lib.barcode_dict.keys())
#     for barcode in barcodes:
#
#         lib_barcodes_dict = barcode_extraction_dict(samfile=y_lib.samfile,
#                                                     reference=reference,
#                                                     barcode_location=y_lib.barcode_dict[barcode],
#                                                     barcode_flank=y_lib.flankseq_dict[barcode])
#
#         lib_corrected_barcode_dict = correct_barcodes_cutoff(lib_barcodes_dict, cutoff=cutoff)
#         print(barcode, "library barcodes:", len(lib_corrected_barcode_dict.keys()))
#
#         lib_df = consensus_df(lib_corrected_barcode_dict, y_file=y_lib, min_coverage=cutoff)
#
#         lib_df['barcode'] = barcode
#
#         if count == 0:
#             # merge df across multiple
#             merged_df = lib_df
#             count += 1
#         else:
#             merged_df = merged_df.append(lib_df)
#
#     return merged_df


# def pipeline_function(lib_df,lib_corrected_barcode_dict,y_sample, reference="Rp0-reference", cutoff=5):
#     """
#     TODO: write bam file for specific UMIset data
#     :param y_lib:
#     :param y_sample:
#     :param cutoff:
#     :return:
#     """
#     barcodes = list(y_lib.barcode_dict.keys())
#
#     count = 0
#
#     for barcode in barcodes:
#
#         sample_barcodes_dict = barcode_extraction_dict(samfile=y_sample.samfile,
#                                                        reference=reference,
#                                                        barcode_location=y_sample.barcode_dict[barcode],
#                                                        barcode_flank=y_sample.flankseq_dict[barcode])
#
#         sample_corrected_barcode_dict = correct_barcodes_cutoff(sample_barcodes_dict, cutoff=cutoff)
#         print(barcode, 'sample barcodes:', len(sample_corrected_barcode_dict.keys()))
#
#         # overlap between the two libraries
#         overlap_barcodes = list(set(sample_corrected_barcode_dict.keys()) & set(lib_corrected_barcode_dict.keys()))
#         print(barcode, "sample and library overlap", len(overlap_barcodes))
#
#         # print(len(lib_corrected_barcode_dict))
#         # print(len(sample_corrected_barcode_dict))
#
#         # select barcodes that are in both samples
#         lib_corrected_barcode_dict = {k: lib_corrected_barcode_dict[k] for k in overlap_barcodes}
#         sample_corrected_barcode_dict = {k: sample_corrected_barcode_dict[k] for k in overlap_barcodes}
#
#         # print(len(lib_corrected_barcode_dict))
#         # print(len(sample_corrected_barcode_dict))
#
#         lib_df = consensus_df(lib_corrected_barcode_dict, y_file=y_lib, min_coverage=cutoff)
#         samp_df = consensus_df(sample_corrected_barcode_dict, y_file=y_sample, min_coverage=cutoff)
#
#         if count == 0:
#             # merge sample and lib df
#             merged_df = pd.merge(samp_df, lib_df, how="inner", on='UMI_pos')
#         else:
#             tmp_df = pd.merge(samp_df, lib_df, how="inner", on='UMI_pos')
#             merged_df = merged_df.append(tmp_df)
#         count += 1
#
#     # total_sites = merged_df.shape[0]
#     # total_mutations = merged_df[merged_df['base_x'] != merged_df['base_y']].shape[0]
#     #
#     # print(total_sites, total_mutations)
#
#     return merged_df

# TODO - add plotting functions
