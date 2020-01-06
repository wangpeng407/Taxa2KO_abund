#!/usr/bin/python
# -*- coding: UTF-8 -*-

import re,os,sys
from itertools import islice
args=sys.argv
def select_key_taxa(kegg_anno, kegg_taxa, keyword):
    with open(kegg_anno, 'r') as f:
        sel_KOs = [line.rstrip().split('\t')[2] for line in islice(f, 1, None) if re.search((keyword), line, re.I)]
        KO_sign = dict(zip(sel_KOs, range(len(sel_KOs))))

    with open(kegg_taxa, 'r') as f:
        lnum = 0
        for line in f:
            if lnum == 0:
                print(line, end='')
            else:
                KO_ID = line.strip().split('\t')[1]
                print(line, end='') if KO_ID in KO_sign else ''
            lnum += 1

def main():
    if len(args) <= 3:
        sys.exit('''Usage: python {} kegg_anno kegg_taxa 'Methane' > sel.Unigenes.KEGG.tax.xls
        '''.format(args[0]))
    #select_key_taxa(kegg_anno=kegg_anno, kegg_taxa=kegg_taxa, keyword='Methane')
    select_key_taxa(kegg_anno=args[1], kegg_taxa=args[2], keyword=args[3])

if __name__ == '__main__':
    main()
