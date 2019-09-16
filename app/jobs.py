import gzip
from sqlalchemy import Column, ForeignKey, Integer, String,DateTime
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import relationship
from sqlalchemy import create_engine
from sqlalchemy.sql import func
from sqlalchemy.orm import sessionmaker
import pandas as pd
import os
import shutil
import subprocess
import os
import shutil
import subprocess

from .models import *

Base = declarative_base()

def flatten_list(l):
    return [item for sublist in l for item in sublist]

def fasta_to_db(url,input_file,file_name,seqrun_id,index_for,index_rev,type_sequencing='illumina',type_instrument='iseq'):
    engine = create_engine(url)
    Base.metadata.create_all(engine)
    DBSession = sessionmaker(bind=engine)
    session = DBSession()

    count = 0

    for i,line in enumerate(input_file):

        if count == 0: # Name, etc
            count+=1
            if type_sequencing == 'illumina':
                new_fastq = Fastq(file_name=file_name,type_instrument=type_instrument,type_sequencing=type_sequencing,docs=line.strip('\n'))
                doc = flatten_list([lst.split(':') for lst in line.split(' ')]) # Split between ' ', then split between ':', then flatten that resulting list
                # Add to class
                new_fastq.seqrun_uuid=seqrun_id
                new_fastq.instrument_id=doc[0]
                new_fastq.run_id=doc[1]
                new_fastq.flow_cell_id=doc[2]
                new_fastq.lane=doc[3]
                new_fastq.tile_number=doc[4]
                new_fastq.x_coord=doc[5]
                new_fastq.y_coord=doc[6]
                new_fastq.member_pair=doc[7]
                new_fastq.read_filter=doc[8]
                new_fastq.control_bits=doc[9]
                new_fastq.defined_index_for = index_for
                new_fastq.defined_index_rev = index_rev

                if len(doc) > 10:
                    indexs = doc[10].strip('\n').split("+")
                    if len(indexs) > 1:
                        new_fastq.index_for=indexs[0]
                        new_fastq.index_rev=indexs[1]
            continue
        elif count == 1: # Sequence
            count+=1
            new_fastq.sequence=line.strip('\n')
            continue
        elif count == 2: # +
            count+=1
            new_fastq.comments=line.strip('\n')
            continue
        elif count == 3: # Quality
            count=0
            new_fastq.read_quality=line.strip('\n')
            session.add(new_fastq)
            if ((i+1)/4)%10000 ==0:
                print('commit')
                session.commit()
    session.commit()
    return 'Completed without error'

def seqrun_bam_generation(url,seqrun_uuid):

    # Get minimap version
    align_vers = str(subprocess.check_output("minimap2 --version",shell=True).decode("utf-8")).rstrip()
    print(align_vers)

    # Create session
    engine = create_engine(url)
    Base.metadata.create_all(engine)
    DBSession = sessionmaker(bind=engine)
    session = DBSession()

    # Grab seqrun
    seqrun = session.query(SeqRun).filter_by(uuid=seqrun_uuid).first() 

    # Get samples
    for obj in seqrun.samples:
        sample_uuid = str(obj.uuid)
        index_for = obj.index_for
        index_rev = obj.index_rev
        bigseq = obj.full_seq

        # TODO: ONLY IMPORT UNIQUE INDEXS, THEN ALL SAMPLES WITH THOSE INDEXS
        print('sql_query')
        sql_query = "SELECT fastqs.uuid,fastqs.sequence,fastqs.comments,fastqs.read_quality FROM fastqs WHERE fastqs.defined_index_for='{}' AND fastqs.defined_index_rev='{}' AND fastqs.seqrun_uuid='{}'"
        with engine.connect() as con:
            rs = con.execute(sql_query.format(index_for,index_rev,seqrun_uuid))
    
            with open('/dev/shm/seq/tmp.fastq', 'w') as f:
                for row in rs:
                    for i,r in enumerate(row):
                        if i == 0:
                            f.write("@{}".format(str(r)) + '\n')
                        else:
                            f.write(r + '\n')
    

        with open('/dev/shm/seq/tmp.fa', 'w') as f:
            f.write('>{}\n'.format(str(obj.sample_uuid)))
            f.write(bigseq)

        print('Generate bam')
        bam_file = subprocess.check_output("minimap2 -a --cs /dev/shm/seq/tmp.fa /dev/shm/seq/tmp.fastq | samtools view -bS -F 4 - | samtools sort - -o /dev/shm/seq/example.bam",shell=True)

        with open("/dev/shm/seq/example.bam","rb") as f:
            obj.bam=f.read()
            session.commit()
        print('Upload complete for {}'.format(sample_uuid))
    seqrun.aligned = True
    session.commit()
    return 'Completed without error'



