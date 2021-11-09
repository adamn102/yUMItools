import sys
import os

#sys.path.append(os.path.realpath('.'))
from yUMItools.utils import *


def test_parse_reference_sequence():
    y = YUMISet('test_data/reference_sequence/Rp0-reference.fa',
                bamfile='test_data/sorted_reads/LRT-TEST-Rp0-INDELMUT-A_S00.bam')
    y.parse_reference_sequence()
    assert len(y.barcode_dict) == 6
    assert y.barcode_dict['barcode0'] == (265, 280)


def test_umi_consensus_indel():
    y = YUMISet('test_data/reference_sequence/Rp0-reference.fa',
                bamfile='test_data/sorted_reads/LRT-TEST-Rp0-INDELMUT-A_S00.bam')
    barcodes = list(y.barcode_dict.keys())
    barcode = barcodes[0]

    lib_barcodes_dict = y.barcode_extraction_dict(barcode)

    corrected_barcode_dict = correct_barcodes_cutoff(lib_barcodes_dict, cutoff=4)

    #barcode = 'ATGTTTCTGGGAGCT'
    barcode = 'ATTACTCGTATGACA'
    #barcode = 'CCACCACGTGCGAAT'
    df = y.consensus_caller(barcode, corrected_barcode_dict,min_coverage=5)

    df = pd.DataFrame({
        'position': position_list,
        'base': sequence_list,
        'coverage': coverage_list,
        'fraction': fraction_list
    })
    df['UMI'] = barcode

    df.to_csv('test_data/test_UMI_consensus.csv')
    print(df[df['position'] == 692]['base'] == 'TT')
    #assert df[df['position'] == 692]['base'] == 'TT'


def test_consensus_df():
    y = YUMISet('test_data/reference_sequence/Rp0-reference.fa',
                bamfile='test_data/sorted_reads/LRT-TEST-Rp0-INDELMUT-A_S00.bam')
    barcodes = list(y.barcode_dict.keys())
    barcode = barcodes[0]

    lib_barcodes_dict = y.barcode_extraction_dict(barcode)

    corrected_barcode_dict = correct_barcodes_cutoff(lib_barcodes_dict, cutoff=5)

    df = y.consensus_df(corrected_barcode_dict, min_coverage=5)
    df.to_csv('test_data/test_UMI_consensus.csv')
    print(df.shape)

def test_library_pipeline():
    y = YUMISet('test_data/reference_sequence/Rp0-reference.fa',
                bamfile='test_data/sorted_reads/LRT-TEST-Rp0-INDELMUT-A_S00.bam')

    df = y.library_pipeline(cutoff=5)
    df.to_csv('test_data/test_UMI_consensus.csv')
    print(df.shape)