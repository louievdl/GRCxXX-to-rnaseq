# mainly cribbed from https://www.biostars.org/p/9562542/

import sys

gene_names = {}
tx_parents = {}

# problems to solve:
# 1. 1645 'gene' records have no 'Name' flag in 9th col
# 2. many non-coding genes of interest exist, should find these too.
# 3. also handle 'transcript' records

# idea, do 2 passes, once to get gene- & tx-like records.
# second pass to 1) filter out exons not part of a named gene, 2) put gene symbols in exon records,
# 3) put gene names on tx records

with open(sys.argv[1], 'r') as gff_input:
    for line in gff_input:
        line = line.rstrip('\n')
        if line.startswith('#'):
            continue
        col = line.split('\t')
        if 'ID=gene:' in str(line):
            try:
                gene_id = col[8].split(';')[0][3:]
                gene_name = col[8].split('Name=')[1].strip()
                gene_name = gene_name.split(';')[0]
                gene_names[gene_id] = gene_name
            except IndexError:
                # if no gene name supplied, use ens gid
                gene_name = gene_id.split(':')[1]
                gene_names[gene_id] = gene_name
        elif 'ID=transcript:' in str(line):
            tx_id = col[8].split(';')[0][3:]
            tx_parent = col[8].split('Parent=')[1].strip()
            tx_parent = tx_parent.split(';')[0]
            tx_parents[tx_id] = tx_parent

with open(sys.argv[1], 'r') as gff_input_2:
    for line in gff_input_2:
        line = line.rstrip('\n')
        col = line.split('\t')
        # if a gene record, print it out
        if line.startswith('#'):
            print(line)
        elif 'ID=transcript:' in str(line):
            tx_id = col[8].split(';')[0][3:]
            parent_id = tx_parents[tx_id]
            gene_name = gene_names[parent_id]
            print(line + ';gene_name=' + gene_name)
        elif col[2] == 'exon':
            parent_id = col[8].split(';')[0][7:] # first element is 'Parent=transcript:ENS...'
            tx_parent = tx_parents[parent_id]
            gene_name = gene_names[tx_parent]
            print(line + ';gene_id=' + tx_parent + ';gene_name=' + gene_name)
        else:
            print(line)

