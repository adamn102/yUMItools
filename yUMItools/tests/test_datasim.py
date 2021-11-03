from yUMItools.datasim import *


def test__read_fasta_file():
    import Bio.Seq

    # filepath
    test_data_fasta_filepath = 'test_data/reference_sequence/Rp0-reference.fa'

    # create sequence object and load fasta file
    s = TestTube(test_data_fasta_filepath)

    # test loading one file
    assert type(s.reference_sequence_dict['Rp0-reference']) is Bio.Seq.Seq

    del s


def test__read_fasta_file_multiple():
    # filepath
    test_data_fasta_filepath = 'test_data/reference_sequence/mixed-reference.fa'

    # create sequence object and load fasta file
    s = TestTube(test_data_fasta_filepath)

    # test loading multiple files
    assert len(s.reference_sequence_dict.keys()) == 2

    del s


def test_find_barcodes():
    # filepath
    test_data_fasta_filepath = 'test_data/reference_sequence/mixed-reference.fa'

    # create sequence object and load fasta file
    s = TestTube(test_data_fasta_filepath)

    library_sequence = s.reference_sequence_dict['Rp0-reference']

    s.find_barcodes(library_sequence=library_sequence)

    assert s.barcode_dict['barcode0'] == (265, 280)
    assert s.barcode_dict['barcode5'] == (3149, 3164)

    del s


def test_generate_barcode_library():
    # filepath
    test_data_fasta_filepath = 'test_data/reference_sequence/Rp0-reference.fa'

    # create sequence object and load fasta file
    s = TestTube(test_data_fasta_filepath)

    library_sequence = s.reference_sequence_dict['Rp0-reference']
    s.find_barcodes(library_sequence=library_sequence)

    s.generate_barcode_library(clones=10)

    assert len(s.barcode_library_list) == 10
    assert hamming(s.barcode_library_list[0], s.barcode_library_list[1]) > 0

    del s


def test_reverse_transcribe_library():
    # filepath
    test_data_fasta_filepath = 'test_data/reference_sequence/Rp0-reference.fa'

    # create sequence object and load fasta file
    s = TestTube(test_data_fasta_filepath)

    library_sequence = s.reference_sequence_dict['Rp0-reference']
    s.find_barcodes(library_sequence=library_sequence)

    # only generate one clone - this allows us to easily count mutations
    s.generate_barcode_library(clones=1)

    # test with no mutations
    rt_library = reverse_transcribe_library(s.barcode_library_list,
                                            clones=2,
                                            mut_type='random',
                                            mutation_rate=0)
    assert hamming(rt_library[0], rt_library[1]) == 0

    # test with low mutation rate
    rt_library = reverse_transcribe_library(s.barcode_library_list,
                                            clones=2,
                                            mut_type='random',
                                            mutation_rate=1. / 10000)
    assert hamming(rt_library[0], rt_library[1]) <= 10

    # test with high mutation rate
    rt_library = reverse_transcribe_library(s.barcode_library_list,
                                            clones=2,
                                            mut_type='random',
                                            mutation_rate=1. / 100)
    assert hamming(rt_library[0], rt_library[1]) >= 10
