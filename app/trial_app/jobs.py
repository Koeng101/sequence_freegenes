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

from models import *

Base = declarative_base()

def flatten_list(l):
    return [item for sublist in l for item in sublist]

def fasta_to_db(url,input_file,file_name,seqrun_id,type_sequencing='illumina',type_instrument='iseq'):
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
                new_fastq.seqrun_id=seqrun_id
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

                indexs = doc[10].strip('\n').split("+")
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

def fastq_to_sam(url,index_for,index_rev,seqrun_id,bigseq):
    engine = create_engine(url)
    Base.metadata.create_all(engine)
    DBSession = sessionmaker(bind=engine)
    session = DBSession()


    sql_query = "SELECT fastqs.docs,fastqs.sequence,fastqs.comments,fastqs.read_quality FROM fastqs WHERE fastqs.index_for='{}' AND fastqs.index_rev='{}' AND fastqs.seqrun_id={}"

    with engine.connect() as con:
        rs = con.execute(sql_query.format(index_for,index_rev,seqrun_id))

        with open('/dev/shm/seq/tmp.fastq', 'w') as f:
            for row in rs:
                for obj in row:
                    f.write(obj + '\n')

    with open('/dev/shm/seq/tmp.fa', 'w') as f:
        f.write('>test\n')
        f.write(bigseq)

    pileup_file = subprocess.check_output("minimap2 -a --cs /dev/shm/seq/tmp.fa /dev/shm/seq/tmp.fastq > /dev/shm/seq/alignment.sam",shell=True)
    print(pileup_file)

    new_samfile = SamFile(seqrun_id=seqrun_id,bigseq=bigseq,alignment_tool='minimap2',alignment_tool_version='2.17-r943-dirty',index_for=index_for,index_rev=index_rev)
    session.add(new_samfile)
    session.commit()

    with open('/dev/shm/seq/alignment.sam','r') as sam_file:
        for i,line in enumerate(sam_file):
            lst = line.split('\t')
            if i > 4: # Skip headers, figure that out later
                n = Sam(samfile_id=new_samfile.id)
                n.qname = lst[0]
                n.flag = lst[1]
                n.rname = lst[2]
                n.pos = lst[3]
                n.mapq = lst[4]
                n.cigar = lst[5]
                n.rnext = lst[6]
                n.pnext = lst[7]
                n.tlen = lst[8]
                n.seq = lst[9]
                n.qual = lst[10]

                for e in lst[11:]:
                    setattr(n,e[0:2],e)
                session.add(n)
                if i%10000==0:
                    print('commit')
                    session.commit()
                    
        session.commit()
    return 'Completed without error'


def sam_to_cls(url,sam_file):
    engine = create_engine(url)
    Base.metadata.create_all(engine)
    DBSession = sessionmaker(bind=engine)
    session = DBSession()

    #new_samfile = SamFile(seqrun_id=seqrun_id,bigseq=bigseq,alignment_tool='minimap2',alignment_tool_version='2.17-r943-dirty',index_for=index_for,index_rev=index_rev)
    for i,line in enumerate(sam_file):
        lst = line.split('\t')
        

        print(len(lst[11:]))




def sequence(run_name:str,big_seq:str,reads:list,tmp_location='./.tmp'):
    # Init by removing
    try:
        shutil.rmtree(tmp_location)
    except OSError:
        print('Failed to delete tmp')
    else:
        print('Deleted tmp')

    # Add directory
    try:
        os.mkdir(tmp_location)
    except OSError:
        print('Failed to create tmp')
    else:
        print('Created tmp')

    print(os.listdir('.'))

    # Run alignment and write combined pileup file
    pileup_loc = './.tmp/{}.pileup'.format(run_name)
    with open('./.tmp/tmp.fasta', 'w') as the_file:
        the_file.write('>tmp_fasta\n{}'.format(big_seq))
    command = 'bwa index ./.tmp/tmp.fasta && bwa mem ./.tmp/tmp.fasta {} | samtools view -bS - | samtools sort - | samtools mpileup -f ./.tmp/tmp.fasta - > {}'.format(' '.join(reads), pileup_loc)
    pileup_file = subprocess.check_output(command,shell=True)

    # Read pileup file into pandas
    combined_pileup = pd.read_table('{}'.format(pileup_loc), names = ["Sequence", "Position", "Reference Base", "Read Count", "Read Results", "Quality"])
    combined_pileup = combined_pileup.set_index('Position')
    return combined_pileup

