import time
import jobs

with open('/dev/shm/seq/tmp.fa','r') as seq_file:
    seq = [row for row in seq_file][1]
start = time.time()
with open('/dev/shm/seq/alignment.sam','r') as input_file:
    jobs.fastq_to_sam('sqlite:///sql.db','TAGGCATG','CTCCTTAC',1,seq)
end = time.time()
print(end - start)
