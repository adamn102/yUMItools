from Bio import SeqIO


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

    def __init__(self, bed_file, reference_sequence):
        self.bed_file = bed_file

        # parse reference sequence
        for record in SeqIO.parse(reference_sequence, "fasta"):
            self.reference_sequence = record
