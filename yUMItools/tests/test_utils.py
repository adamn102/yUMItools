import sys
import os

sys.path.append(os.path.realpath('.'))
from yUMItools.utils import *

def test_parse_reference_sequence():
    y = YUMISet('test_data/reference_sequence/LRT-REJOc-reference.fa')
    y.parse_reference_sequence()
    assert len(y.barcode_dict) == 6
    assert y.barcode_dict['barcode0'] == (265, 280)