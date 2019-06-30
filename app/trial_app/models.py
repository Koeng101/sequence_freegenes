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

    seqrun_id = Column(Integer, ForeignKey('seqruns.id'), nullable=False)
    file_name = Column(String)

    # Shared
    instrument_id = Column(String)
    flow_cell_id = Column(String)
    run_id = Column(String)

    # Types
    type_instrument = Column(String)
    type_sequencing = Column(String)

    # Fastq
    docs = Column(String)
    sequence = Column(String)
    comments = Column(String)
    read_quality = Column(String)

    ## Illumina   
    lane = Column(Integer())
    tile_number = Column(Integer())
    x_coord = Column(Integer())
    y_coord = Column(Integer())
    member_pair = Column(Integer())
    read_filter = Column(String) # Y or N
    control_bits = Column(Integer())
    index_for = Column(String)
    index_rev = Column(String)

class SamFile(Base):
    __tablename__ = 'samfiles'
    id = Column(Integer, primary_key=True)
    time_created = Column(DateTime(timezone=True), server_default=func.now())

    bigseq = Column(String)
    alignment_tool = Column(String)
    alignment_tool_version = Column(String)

class Sam(Base):
    __tablename__ = 'sams'
    id = Column(Integer, primary_key=True)
    time_created = Column(DateTime(timezone=True), server_default=func.now())

    # https://en.wikipedia.org/wiki/SAM_(file_format)
    # https://samtools.github.io/hts-specs/SAMv1.pdf
    qname = Column(String)
    flag = Column(Integer())
    rname = Column(String)
    pos = Column(Integer())
    mapq = Column(Integer())
    cigar = Column(String)
    rnext = Column(String)
    pnext = Column(Integer())
    tlen = Column(Integer())
    seq = Column(String)
    qual = Column(String)
    index_for = Column(String)
    index_rev = Column(String)

    alignment_tool = Column(String)
    alignment_tool_version = Column(String)


class SeqRun(Base):
    __tablename__ = 'seqruns'
    id = Column(Integer, primary_key=True)
    time_created = Column(DateTime(timezone=True), server_default=func.now())
    time_updated = Column(DateTime(timezone=True), onupdate=func.now())

    name = Column(String) # run name
    notes = Column(String)
    run_id = Column(String) # Sequencing provider id
    machine_id = Column(String)
    sequencing_type = Column(String) # illumina, nanopore, etc
    machine = Column(String) # minion, iseq, etc
    provider = Column(String) # in-house

    fastqs = relationship('Fastq',backref='seqrun')

class PileupFile(Base):
    __tablename__ = 'pileupfiles'
    id = Column(Integer, primary_key=True)
    time_created = Column(DateTime(timezone=True), server_default=func.now())

    status = Column(String) # mutation,confirmed,etc
    full_search_sequence = Column(String)
    target_sequence = Column(String)
    index_for = Column(String)
    index_rev = Column(String)

    sample_uuid = Column(String)
    seqrun_id = Column(Integer, ForeignKey('seqruns.id'), nullable=False)

class PileupLine(Base):
    __tablename__ = 'pileups'
    id = Column(Integer, primary_key=True)
    time_created = Column(DateTime(timezone=True), server_default=func.now())

    # https://en.wikipedia.org/wiki/Pileup_format
    pileup_id = Column(Integer, ForeignKey('pileups.id'), nullable=False)
    sequence = Column(String)
    position = Column(Integer)
    reference_base = Column(String)
    read_count = Column(Integer)
    read_results = Column(String)
    quality = Column(String)

class Job(Base):
    __tablename__ = 'jobs'
    id = Column(Integer, primary_key=True)
    time_created = Column(DateTime(timezone=True), server_default=func.now())

    job_id = Column(String)
    status = Column(String)
    description = Column(String)

