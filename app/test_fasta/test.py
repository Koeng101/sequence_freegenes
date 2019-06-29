import gzip
from sqlalchemy import Column, ForeignKey, Integer, String,DateTime
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import relationship
from sqlalchemy import create_engine
from sqlalchemy.sql import func
from sqlalchemy.orm import sessionmaker

Base = declarative_base()

from redis import Redis
from rq import Queue
q = Queue(connection=Redis(), default_timeout=7200)

from fastq_process import fasta_to_db

json_file = {'file_name': 'merge-all_S1_L001_R2_001.fastq.gz'}
input_file = '/home/koeng/gits/sequence_freegenes/example/merge-all_S1_L001_R1_001.fastq.gz'
with gzip.open(input_file,'rt') as fin:
    full_file = fin.read()
job = q.enqueue(fasta_to_db, args=('sqlite://sql.db',full_file,json_file['file_name']))

count = 0
print(job.get_status())
while job.get_status() != 'finished':
    time.sleep(1)
    count+=1
    if count%10 == 0:
        print(job.get_status())
        print(count)
print(job.result)


