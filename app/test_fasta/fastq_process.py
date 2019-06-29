import gzip
from sqlalchemy import Column, ForeignKey, Integer, String,DateTime
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import relationship
from sqlalchemy import create_engine
from sqlalchemy.sql import func
from sqlalchemy.orm import sessionmaker
Base = declarative_base()

class Fastq(Base):
    __tablename__ = 'fastqs'
    id = Column(Integer, primary_key=True)
    time_created = Column(DateTime(timezone=True), server_default=func.now())

    file_name = Column(String)

    # Shared
    instrument_id = Column(String())
    flow_cell_id = Column(String())
    run_id = Column(String())

    # Types
    type_instrument = Column(String())
    type_sequencing = Column(String())

    # Fastq
    docs = Column(String)
    sequence = Column(String())
    comments = Column(String())
    read_quality = Column(String())

    ## Illumina   
    lane = Column(Integer())
    tile_number = Column(Integer())
    x_coord = Column(Integer())
    y_coord = Column(Integer())
    member_pair = Column(Integer())
    read_filter = Column(String()) # Y or N
    control_bits = Column(Integer())
    index_for = Column(String())
    index_rev = Column(String())


def flatten_list(l):
    return [item for sublist in l for item in sublist]

def fasta_to_db(url,input_file,file_name,type_sequencing='illumina',type_instrument='iseq'):
    engine = create_engine(url)
    Base.metadata.create_all(engine)
    DBSession = sessionmaker(bind=engine)
    session = DBSession()

    count = 0
    for i,line in enumerate(input_file):
        if count == 0: # Name, etc
            count+=1
            if type_sequencing == 'illumina':
                new_fastq = Fastq(file_name=file_name,type_instrument=type_instrument,type_sequencing=type_sequencing,docs=line)
                doc = flatten_list([lst.split(':') for lst in line.split(' ')]) # Split between ' ', then split between ':', then flatten that resulting list
                # Add to class
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
                session.commit()
    session.commit()
    return 'Completed without error'


