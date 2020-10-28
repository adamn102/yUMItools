#!/bin/bash

#import packages
from yUMItools.utils import *
import os
from umi_tools import UMIClusterer
import io
from contextlib import redirect_stderr


#create YUMIset object - with refernece sequnce and bam file
y = YUMISet(reference_sequence = 'data/reference_sequence/2880-reference.fa',
           bamfile ='data/sorted_reads/p2880-insert-v1_S70.bam')