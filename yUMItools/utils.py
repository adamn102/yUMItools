import pysam
from collections import defaultdict, Counter
import numpy as np
import pandas as pd
from collections import Counter

chars = 'ACGTN'


class YUMISet:
    """
    yUMI set Data Class

    This is the primary data structure for the yUMItools software

    Main functionality will include:

    * Parsing multiple bam data sets for any one bed formatted reference data set
    * error correcting UMIs
    * phasing multiple UMIs (maybe)
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
        This function parses reference sequences for UMI locations
        and reference mapping sequences

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

    def read_pair_dict(self):
        """
        Takes in the bam file passed into the yUMISet object
        and parses the reads into a dict object
        with the read names as the key and the read data as the values
        this can get rather large - maybe worth deleting once done.
        :return:
        """

        for read in self.samfile.fetch(None):
            # filter out unpaired reads
            if read.is_paired:
                # filer out non mapped reads
                if not read.is_unmapped:
                    # make sure the read object has the aligned_pairs attribute
                    if hasattr(read, "aligned_pairs"):
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


    def library_pipeline(self, cutoff=5, subset=False, sub_value=1000):
        count = 0
        barcodes = list(self.barcode_dict.keys())

        for barcode in barcodes:

            lib_barcodes_dict = self.barcode_extraction_dict(barcode)

            corrected_barcode_dict = correct_barcodes_cutoff(lib_barcodes_dict, cutoff=cutoff)
            print(barcode, "library barcodes:", len(corrected_barcode_dict.keys()))

            umi_df = self.consensus_df(corrected_barcode_dict, min_coverage=cutoff, subset=subest, sub_value=sub_value)

            umi_df['barcode'] = barcode

            if self.umi_df.empty == True:
                # merge df across multiple
                self.umi_df = umi_df
            else:
                self.umi_df = self.umi_df.append(umi_df)

        return self.umi_df

    def consensus_df(self, barcode_dict, min_coverage=5,subset=False, sub_value_=1000):
        count = 0
        for umi in barcode_dict.keys():
            # print(umi)
            if len(barcode_dict[umi]) >= min_coverage:
                if count == 0:
                    df = self.consensus_caller(umi, barcode_dict, min_coverage, subest, sub_value)
                    if type(df) != int:
                        count += 1
                else:
                    df2 = self.consensus_caller(umi, barcode_dict, min_coverage, subest, sub_value)
                    if type(df2) != int:
                        df = df.append(df2)

        # sample_df = df
        df['position'] = df['position'].astype(int)
        # df['UMI_pos'] = df['UMI'] + "_" + df['position'].astype(str)
        return df

    def consensus_caller(self, barcode, corrected_barcode_dict, min_coverage=5, subset=False, sub_value=1000):

        # find the reads (sequences) that correspond to a specific barcode
        read_list = corrected_barcode_dict[str(barcode)]
        # print(len(read_list))

        # if there are too many reads, take a subset of them
        if subset:
            if len(read_list) >= sub_value:
                read_list = np.random.choice(read_list, sub_value, replace=False)

        # print("modified_read_list: ", len(read_list))
        outer_count = 0
        # print(len(read_list))
        for read_name in read_list:
            # print(read_name)
            for read in self.read_dict[read_name]:

                # print(len(read.positions), len(read.seq))

                # print(hasattr(read, "aligned_pairs"))
                if hasattr(read, "aligned_pairs"):

                    align_array = np.array(read.aligned_pairs)

                    if align_array.ndim == 2:

                        # print(barcode, align_array.ndim ,align_array.shape, len(read.query_sequence), read.is_unmapped)

                        # find deletions
                        deletions = list(np.where(align_array[:, 0] == None)[0])

                        # find insertions
                        insertions = list(np.where(align_array[:, 1] == None)[0])

                        # print(insertions)

                        # convert bases and quals into arrays
                        quals = np.array(list(read.qual), dtype=object)
                        bases = np.array(list(read.seq), dtype=object)

                        # handle deletions
                        for deletion in deletions:
                            # print(deletion)
                            quals = np.append(
                                np.append(quals[0:deletion], np.array(
                                    ["N"], dtype=object), 0), quals[deletion:], 0)
                            bases = np.append(
                                np.append(bases[0:deletion], np.array(
                                    ["N"], dtype=object), 0), bases[deletion:], 0)

                        # handle insertions

                        #
                        # for each insertion, count the insertion length, store as a dict
                        index_array = np.array([(i - x) for i, x in enumerate(insertions)]).reshape(-1)
                        data_array = np.array(insertions).reshape(-1)
                        a = np.stack((index_array, data_array), axis=1)
                        insertion_lengths = []
                        insertion_positions = []
                        n = np.unique(a[:, 0])
                        for i in n:
                            insertion_positions.append(a[a[:, 0] == i, 1][0])
                            insertion_lengths.append(len(a[a[:, 0] == i, 1]))
                        insertion_dict = dict(zip(insertion_positions, insertion_lengths))
                        # print(insertion_dict)

                        for k, v in insertion_dict.items():
                            # print(k, v)
                            # print(align_array)
                            align_array[k,][1] = int(align_array[k + v,][1])
                            # print(align_array)

                            insert_bases = "".join(bases[k:k + v + 1])
                            # print(insert_bases)
                            bases[k] = insert_bases

                            insert_quals = "".join(
                                quals[k:k + v])
                            quals[k] = insert_quals

                        deletion_sites = []
                        for k, v in insertion_dict.items():
                            count = 1
                            for i in range(v):
                                deletion_sites.append(k + count)
                                count += 1

                        # print(deletion_sites)
                        if len(deletion_sites) > 0:
                            # print(deletion_sites)
                            bases = np.delete(bases, deletion_sites)
                            quals = np.delete(quals, deletion_sites)
                            align_array = np.delete(align_array, deletion_sites, axis=0)

                        quals = quals.reshape(len(quals), 1)
                        bases = bases.reshape(len(bases), 1)

                        if outer_count == 0:
                            arr = np.concatenate((align_array, bases, quals), axis=1)
                            outer_count += 1
                        else:
                            tmp_arr = np.concatenate((align_array, bases, quals), axis=1)
                            arr = np.concatenate((arr, tmp_arr), axis=0)

        # sort array on position
        arr[:, 1] = arr[:, 1].astype(int)
        arr = arr[arr[:, 1].argsort()]

        # split array by position
        # find unique positions
        sequence_list = []
        position_list = []
        coverage_list = []
        fraction_list = []

        unique_positions, position_index = np.unique(arr[:, 1], return_index=True)
        split_arr = np.split(arr, position_index[1:])
        # print(split_arr)

        for position in range(len(split_arr)):

            coverage = len(split_arr[position])

            if coverage >= min_coverage:
                onehot_group = split_arr[position][:, 0:4]
                unique, counts = np.unique(onehot_group[:, 2], return_counts=True)

                # print(position, unique, counts)
                # base_count = np.bincount(onehot_group[:, 2].astype(int))
                max_base = unique[np.argmax(counts)]
                top_base_acc = counts[np.argmax(counts)] / sum(counts)
                # print(max_base, top_base_acc)
                # consensus_bases = int_to_char[max_base]

                position_list.append(split_arr[position][:, 1][0])
                sequence_list.append(max_base)
                coverage_list.append(coverage)
                fraction_list.append(top_base_acc)

        # return position_list, sequence_list, coverage_list, fraction_list

        df = pd.DataFrame({
            'position': position_list,
            'base': sequence_list,
            'coverage': coverage_list,
            'fraction': fraction_list
        })
        df['UMI'] = barcode
        return df

    def write_umi_data_to_bam(seread_namelf, barcode, umi, output_file):
        """
        Write out the mapped bam data for any specific UMI -
        This is useful for visually inspecting mutations in a specific bam file
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
