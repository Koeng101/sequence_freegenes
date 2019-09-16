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

def seqrun_sam_generation(url,seqrun_uuid):

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
                    for i,obj in enumerate(row):
                        if i == 0:
                            f.write("@{}".format(str(obj)) + '\n')
                        else:
                            f.write(obj + '\n')
    

        with open('/dev/shm/seq/tmp.fa', 'w') as f:
            f.write('>test\n')
            f.write(bigseq)

        print('subprocess check')
        sam_file = subprocess.check_output("minimap2 -a --cs /dev/shm/seq/tmp.fa /dev/shm/seq/tmp.fastq | samtools view -b -F 4 - | samtools sort - -O sam",shell=True).decode("utf-8").rstrip()
        sams = []

        with open('/dev/shm/seq/example.sam', 'w') as f:
            f.write(sam_file)

        #print('start')
        #for i,line in enumerate(sam_file.split('\n')):
        #    lst = line.split('\t')
        #    if i > 4: # Skip headers, figure that out later
        #        new_sam = {
        #                "sample_uuid": sample_uuid,
        #                "fastq_uuid": lst[0],
        #                "alignment_tool": "minimap2",
        #                "alignment_tool_version": "",
        #                "flag": lst[1],
        #                "pos": lst[3],
        #                "mapq": lst[4],
        #                "cigar": lst[5],
        #                "rnext": lst[6],
        #                "pnext": lst[7],
        #                "tlen": lst[8],
        #                }
        #        for e in lst[11:]:
        #            new_sam[e[0:2].lower()] = e
        #        sams.append(new_sam)
        #    if i %5000 == 0:
        #        print('Working')
        #session.bulk_insert_mappings(Sam, sams)
        #session.commit()
        print('Upload complete for {}'.format(sample_uuid))
    seqrun.aligned = True
    session.commit()
    return 'Completed without error'



