import pysam
from collections import defaultdict, Counter


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
    """

    def __init__(self, reference_sequence, bamfile):
        self.reference_sequence = reference_sequence
        self.barcode_dict = dict()
        self.flankseq_dict = dict()
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


def fetch_barcode_reads(barcode, samfile, reference, barcode_location, barcode_flank, distance=2):
    # goal of this function is to fetch paired reads to specific barcode regions
    # for any specific barcode will fetch the paired reads

    # list for collecting barcodes
    barcode_dict = dict()

    bam_iter = samfile.fetch(reference, barcode_location[0], barcode_location[1])
    for x in bam_iter:
        if barcode_flank[0] in x.seq:
            if barcode_flank[1] in x.seq:
                start = x.seq.find(barcode_flank[0][10:]) + 5
                end = x.seq.find(barcode_flank[1][:5])
                if start < end and len(x.seq[start:end]) == 15:
                    umi = x.seq[start:end]
                    if hamming(umi, barcode) <= distance:
                        if umi not in barcode_dict.keys():
                            barcode_dict[barcode] = [x.query_name]

                        else:
                            barcode_dict[barcode].append(x.query_name)
    return barcode_dict[barcode]


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
