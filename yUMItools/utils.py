import pysam


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

    def __init__(self, reference_sequence):
        self.reference_sequence = reference_sequence
        self.barcode_dict = dict()
        self.flankseq_dict = dict()
        self.value = 0

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


def read_bamfile(file_name):

    samfile = pysam.AlignmentFile(file_name, "rb")

    return samfile


def barcode_extraction(samfile, reference, barcode_location, barcode_flank):
    # list for collecting barcodes
    barcodes = []

    bam_iter = samfile.fetch(reference, barcode_location[0], barcode_location[1])
    for x in bam_iter:
        if barcode_flank[0] in x.seq:
            if barcode_flank[1] in x.seq:
                start = x.seq.find(barcode_flank[0][10:]) + 5
                end = x.seq.find(barcode_flank[1][:5])
                if start < end and len(x.seq[start:end]) == 15:
                    barcodes.append(x.seq[start:end])
    return barcodes

def fetch_barcodes(barcode, samfile, reference, barcode_location, barcode_flank):

    # goal of this function is to fetch paired reads to specific barcode regions
    # for any specific barcode will fetch the paired reads


    pass