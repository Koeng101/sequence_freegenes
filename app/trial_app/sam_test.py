import time
import jobs
import pandas as pd

csv_file = '/home/koeng/gits/sequence_freegenes/example/merge-all.csv'
df = pd.read_csv(csv_file)
big_seq = []
for index,row in df.iterrows():
    big_seq.append(row['forward_insert'])
    big_seq.append(row['seq'])
    big_seq.append(row['reverse_insert'])
big_seq = ''.join(big_seq)

start = time.time()
jobs.fastq_to_sam('sqlite:///sql.db','TAGGCATG','CTCCTTAC',1,big_seq)
end = time.time()
print(end - start)
