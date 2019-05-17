import requests
from redis import Redis
from rq import Queue
import pandas as pd
import os
import shutil
import subprocess
import time

from seq import sequence, test_print
# Init rq
q = Queue(connection=Redis(), default_timeout=7200)




csv_file='/home/koeng/gits/sequence_freegenes/merge-all.csv'
df = pd.read_csv(csv_file)
big_seq = []
for index,row in df.iterrows():
    big_seq.append(row['forward_insert'])
    big_seq.append(row['seq'])
    big_seq.append(row['reverse_insert'])
big_seq = ''.join(big_seq)

reads = ['/home/koeng/gits/sequence_freegenes/merge-all_S1_L001_R2_001.fastq.gz','/home/koeng/gits/sequence_freegenes/merge-all_S1_L001_R1_001.fastq.gz']

job = q.enqueue(sequence, args=('test',big_seq,reads))
#job = q.enqueue(test_print, args=('hello world',))
#print(sequence('test',big_seq,reads))

count = 0
print(job.get_status())
while job.get_status() != 'finished':
    time.sleep(1)
    count+=1
    if count%10 == 0:
        print(job.get_status())
        print(count)

print(job.result)


