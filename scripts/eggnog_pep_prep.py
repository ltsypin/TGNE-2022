from Bio import SeqIO

with open('./pep.fasta', 'r') as f:
    seq_list = list(SeqIO.parse(f, 'fasta'))

SeqIO.write(seq_list, './pep.fasta', 'fasta')

# Reading and writing the file ensures that the gene names are correctly annotated within the file

exit()