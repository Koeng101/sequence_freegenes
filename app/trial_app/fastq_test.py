import gzip
from sqlalchemy import Column, ForeignKey, Integer, String,DateTime
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import relationship
from sqlalchemy import create_engine
from sqlalchemy.sql import func
from sqlalchemy.orm import sessionmaker

import sys
from redis import Redis
from rq import Queue
q = Queue(connection=Redis(), default_timeout=7200)

from jobs import fasta_to_db
from models import *

# Add seqrun
engine = create_engine('sqlite:///sql.db')
Base.metadata.create_all(engine)
DBSession = sessionmaker(bind=engine)
session = DBSession()
new_seqrun = SeqRun(name='Test1')
session.add(new_seqrun)
session.commit()



# Fasta to db
json_file = {'file_name': 'merge-all_S1_L001_R2_001.fastq.gz'}
input_file = '/home/koeng/gits/sequence_freegenes/example/merge-all_S1_L001_R2_001.fastq.gz'
with gzip.open(input_file,'rt') as fin:
    full_file = fin.readlines()
    print(type(full_file))
    fasta_to_db('sqlite:///sql.db',full_file,json_file['file_name'],new_seqrun.id)
#job = q.enqueue(fasta_to_db, args=('sqlite:///sql.db',full_file,json_file['file_name']))

sys.exit()
count = 0
print(job.get_status())
while job.get_status() != 'finished':
    time.sleep(1)
    count+=1
    if count%10 == 0:
        print(job.get_status())
        print(count)
print(job.result)


