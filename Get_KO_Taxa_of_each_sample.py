#!/usr/bin/python
# -*- coding: UTF-8 -*-
# import argparse
import os
import sys
import time
from itertools import islice
from multiprocessing import Process

args=sys.argv

def one_sam_func(sample, id_value, allKOs, allTaxon, taxa_gene, KO_gene, outpath):
    t0 = time.time()
    abs_path = os.path.abspath(outpath)
    if not os.path.isdir(abs_path):
        os.mkdir(abs_path)
    out_file_name1 = abs_path + '/' + sample + '-KO-TAXA.list'
    otu_file_name2 = abs_path + '/' + sample + '-KO-TAXA.xls'
    #out1 = open(out_file_name1, 'w')
    #out1.write("Taxa\tKO_ID\tabundance\n")
    out2 = open(otu_file_name2, 'w')
    out2.write("Taxa" + '\t' + '\t'.join(allKOs) + "\n")
    for Taxa in allTaxon:
        out2.write(Taxa)
        for KO in allKOs:
            uniq_gene_id = set(taxa_gene[Taxa]).intersection(set(KO_gene[KO]))
            gene_abun = [float(id_value[g]) for g in uniq_gene_id if g in id_value]
            sum_abun = sum(gene_abun)
     #       out1.write(Taxa + "\t" + KO + "\t" + str(sum_abun) + "\n")
            out2.write("\t" + str(sum_abun))
        out2.write("\n")
   # out1.close()
    out2.close()
    print("...{} DONE, using {} s.".format(sample, time.time() - t0))


def Get_KOs_for_Taxa(kegg_taxa, even_reads, outpath):
    with open(kegg_taxa, 'r') as f:
        taxa_gene = {}
        taxa_KOs = {}
        KO_gene = {}
        for line in islice(f, 1, None):
            l = line.rstrip().split('\t')
            Taxa, Gene_id, KO_id = l[-1], l[-2], l[1]

            if Taxa not in taxa_gene:
                taxa_gene[Taxa] = []
                taxa_gene[Taxa].append(Gene_id)
            else:
                taxa_gene[Taxa].append(Gene_id)

            if KO_id not in KO_gene:
                KO_gene[KO_id] = []
                KO_gene[KO_id].append(Gene_id)
            else:
                KO_gene[KO_id].append(Gene_id)

            if Taxa not in taxa_KOs:
                taxa_KOs[Taxa] = []
                taxa_KOs[Taxa].append(KO_id)
            else:
                taxa_KOs[Taxa].append(KO_id)

    allKOs = [KO for KO in KO_gene.keys()]
    allTaxon = [Taxa for Taxa in taxa_gene]

    with open(even_reads, 'r') as f:
        sam_abun = {}
        lnum = 0
        for line in f:
            # print(lnum)
            if lnum == 0:
                head = line.rstrip().split('\t')
                samples = head[1:]
            else:
                l = line.rstrip().split('\t')
                for index in range(len(samples)):
                    if samples[index] not in sam_abun:
                        sam_abun[samples[index]] = {}
                        sam_abun[samples[index]][l[0]] = l[index + 1]
                    else:
                        sam_abun[samples[index]][l[0]] = l[index + 1]
            lnum += 1

    for sam, id_value in sam_abun.items():
        p = Process(target=one_sam_func, args=(sam,id_value, allKOs, allTaxon, taxa_gene, KO_gene, outpath,))
        p.start()
    p.join()

def main():
    if len(args) < 4:
        sys.exit('''Usage: python {} test_Unigenes.KEGG.tax.xls test_Unigenes.readsNum.even.xls outpath'''.format(args[0]))
    #abs_path = 'C:\\Users\\wangpeng\\Desktop'
    #f1 = abs_path + '/' + 'test_Unigenes.KEGG.tax.xls'
    #f2 = abs_path + '/' + 'test_Unigenes.readsNum.even.xls'
    #Get_KOs_for_Taxa(f1, f2, abs_path)
    Get_KOs_for_Taxa(args[1], args[2], args[3])

if __name__ == '__main__':
    t0 = time.time()
    main()
    print('All completed using {} s'.format(time.time()-t0))
