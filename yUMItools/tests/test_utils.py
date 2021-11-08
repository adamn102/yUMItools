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


def test_UMI_consensus_indel():
    y = YUMISet('test_data/reference_sequence/Rp0-reference.fa',
                bamfile='test_data/sorted_reads/LRT-TEST-Rp0-INDELMUT-A_S00.bam')
    barcodes = list(y.barcode_dict.keys())
    barcode = barcodes[0]

    lib_barcodes_dict = y.barcode_extraction_dict(barcode)

    corrected_barcode_dict = correct_barcodes_cutoff(lib_barcodes_dict, cutoff=4)

    barcode = 'ATGTTTCTGGGAGCT'
    position_list, sequence_list, coverage_list, fraction_list = y.UMI_consensus_indel(barcode, corrected_barcode_dict,min_coverage=1)

    df = pd.DataFrame({
        'position': position_list,
        'base': sequence_list,
        'coverage': coverage_list,
        'fraction': fraction_list
    })
    df['UMI'] = barcode

    df.to_csv('test_data/test_UMI_consensus.csv')
    print(df)

