from Bio import SeqIO

with open('pep.fasta', 'r') as f:
    pep_records = list(SeqIO.parse(f, 'fasta'))

    num_total_records = len(pep_records)

    odd_pep_record_idxs = []

    for idx, r in enumerate(pep_records):
        r.id = r.id.split('|')[0]  # Clean up the record ID

        if '*' in r.seq:
            astr_idxs = [i for i, c in enumerate(r.seq) if c == '*']

            if len(astr_idxs) > 1 or len(r.seq) - 1 != astr_idxs[-1]:
                print(r.id, f'has a \'*\' at amino acid {", ".join([str(i) for i in astr_idxs])}.')
                odd_pep_record_idxs.append(idx)
            else:
                # Remove trailing '*' from the sequence
                r.seq = r.seq[:-1]

    for i in sorted(odd_pep_record_idxs, reverse=True):
        odd_pep_record_seq_split = pep_records[i].seq.split('*')
        # Find the longest subsequence
        pep_records[i].seq = max(odd_pep_record_seq_split, key=len)

    print(' TOTAL RECORDS:', num_total_records)
    print('NORMAL RECORDS:', len(pep_records) - len(odd_pep_record_idxs))
    print('   ODD RECORDS:', len(odd_pep_record_idxs))

    SeqIO.write(pep_records, 'pep_cleaned.fasta', 'fasta')

exit()