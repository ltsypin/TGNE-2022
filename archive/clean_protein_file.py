from Bio import SeqIO

with open('./Tetrahymena_Genome_annotation_V2024_Protein_addAnno.fasta', 'r') as f:
    pep_records = list(SeqIO.parse(f, 'fasta'))

    num_total_records = len(pep_records)

    odd_pep_record_idxs = []

for idx, r in enumerate(pep_records):
    r.id = r.id.split('|')[0]

    if '*' in r.seq:

        astr_idxs = [idx for idx, c in enumerate(r.seq) if c == '*']

        if len(astr_idxs) > 1 or len(r.seq) - 1 != astr_idxs[len(astr_idxs) - 1]:
            print(r.id, f'has a \'*\' at amino acid {", ".join([str(i) for i in astr_idxs])}.')
            odd_pep_record_idxs.append(idx)

        else:
            r.seq = r.seq[:len(r.seq) - 1]

odd_pep_records = []

for i in sorted(odd_pep_record_idxs, reverse=True):
    odd_pep_records.append(pep_records.pop(i))

print(' TOTAL RECORDS:', num_total_records)
print('NORMAL RECORDS:', len(pep_records))
print('   ODD RECORDS:', len(odd_pep_records))

SeqIO.write(pep_records, './pep.fasta', 'fasta')
SeqIO.write(odd_pep_records, './odd_pep.fasta', 'fasta')
