#!/usr/bin/env python
# coding: utf-8

import Bio
from Bio import SeqIO
import pandas as pd
from pandas import DataFrame
import re
import argparse


parser = argparse.ArgumentParser()
parser.add_argument(
    '-f',help='input fasta file providing protein sequence. e.g.UP000005640_9606.fasta')
parser.add_argument(
    '-a',help='requried aa type')
parser.add_argument(
    '-o',help='output file with specific aa position in sequence')
args = vars(parser.parse_args())


getpos = []
for record in SeqIO.parse(args['f'],
                             "fasta"):
    for i in re.finditer(args['a'],str(record.seq)):
        pos = i.span()[1]
        seq_id = record.id.split('|')
        seq_gn = [x[3:] for x in record.description.split(' ') if x.find ('GN=') >= 0]
        getpos.append([seq_id[1],seq_id[2],seq_gn[0],pos])
aa_pos = DataFrame(getpos,columns = ['id','pro_id','gene_name','position'])
aa_pos.to_csv(args['o'],sep='\t',index=False)

