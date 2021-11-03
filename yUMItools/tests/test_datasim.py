from yUMItools.datasim import *


def test__read_fasta_file():
    import Bio.Seq

    # filepath
    test_data_fasta_filepath = 'test_data/reference_sequence/Rp0-reference.fa'

    # create sequence object and load fasta file
    s = testtube(test_data_fasta_filepath)

    # test loading one file
    assert type(s.reference_sequence_dict['Rp0-reference']) is Bio.Seq.Seq


def test__read_fasta_file_multiple():
    # filepath
    test_data_fasta_filepath = 'test_data/reference_sequence/mixed-reference.fa'

    # create sequence object and load fasta file
    s = testtube(test_data_fasta_filepath)

    # test loading multiple files
    assert len(s.reference_sequence_dict.keys()) == 2


def test_find_barcodes():
    # filepath
    test_data_fasta_filepath = 'test_data/reference_sequence/mixed-reference.fa'

    # create sequence object and load fasta file
    s = testtube(test_data_fasta_filepath)

    library_sequence = s.reference_sequence_dict['Rp0-reference']

    s.find_barcodes(library_sequence=library_sequence)

    assert s.barcode_dict['barcode0'] == (265, 280)
    assert s.barcode_dict['barcode5'] == (3149, 3164)
