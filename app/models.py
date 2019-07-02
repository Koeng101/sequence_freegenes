import gzip
#from sqlalchemy import Column, db.ForeignKey, db.Integer, db.String,db.DateTime
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import relationship
from sqlalchemy import create_engine
from sqlalchemy.sql import func
from sqlalchemy.orm import sessionmaker

from flask_sqlalchemy import SQLAlchemy

db = SQLAlchemy()

from jsonschema import validate
import json
import string

def schema_generator(properties,required,additionalProperties=False):
    return {"$schema": "http://json-schema.org/schema#",
            "type": "object",
            "properties": properties,
            "required": required,
            "additionalProperties": additionalProperties}

uuid_regex = '^[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}$'
uuid_schema = {'type': 'string','pattern': uuid_regex}
generic_string = {'type': 'string'}
generic_num = { "type": "number" }

def base_dict(cls):
    return {"id": cls.id, "time_created": cls.time_created.isoformat()}



class Fastq(db.Model):
    __tablename__ = 'fastqs'
    id = db.Column(db.Integer, primary_key=True)
    time_created = db.Column(db.DateTime(timezone=True), server_default=func.now())

    seqrun_id = db.Column(db.Integer, db.ForeignKey('seqruns.id'), nullable=False)
    file_name = db.Column(db.String)

    # Shared
    instrument_id = db.Column(db.String)
    flow_cell_id = db.Column(db.String)
    run_id = db.Column(db.String)

    # Types
    type_instrument = db.Column(db.String)
    type_sequencing = db.Column(db.String)

    # Fastq
    docs = db.Column(db.String)
    sequence = db.Column(db.String)
    comments = db.Column(db.String)
    read_quality = db.Column(db.String)

    ## Illumina   
    lane = db.Column(db.Integer())
    tile_number = db.Column(db.Integer())
    x_coord = db.Column(db.Integer())
    y_coord = db.Column(db.Integer())
    member_pair = db.Column(db.Integer())
    read_filter = db.Column(db.String) # Y or N
    control_bits = db.Column(db.Integer())
    index_for = db.Column(db.String)
    index_rev = db.Column(db.String)

class SamFile(db.Model):
    __tablename__ = 'samfiles'
    id = db.Column(db.Integer, primary_key=True)
    time_created = db.Column(db.DateTime(timezone=True), server_default=func.now())

    seqrun_id = db.Column(db.Integer, db.ForeignKey('seqruns.id'), nullable=False)

    bigseq = db.Column(db.String)
    alignment_tool = db.Column(db.String)
    alignment_tool_version = db.Column(db.String)
    index_for = db.Column(db.String)
    index_rev = db.Column(db.String)


class Sam(db.Model):
    __tablename__ = 'sams'
    id = db.Column(db.Integer, primary_key=True)
    time_created = db.Column(db.DateTime(timezone=True), server_default=func.now())

    samfile_id = db.Column(db.Integer, db.ForeignKey('samfiles.id'), nullable=False)
    # https://en.wikipedia.org/wiki/SAM_(file_format)
    # https://samtools.github.io/hts-specs/SAMv1.pdf
    qname = db.Column(db.String)
    flag = db.Column(db.Integer())
    rname = db.Column(db.String)
    pos = db.Column(db.Integer())
    mapq = db.Column(db.Integer())
    cigar = db.Column(db.String)
    rnext = db.Column(db.String)
    pnext = db.Column(db.Integer())
    tlen = db.Column(db.Integer())
    seq = db.Column(db.String)
    qual = db.Column(db.String)

    # In spec
    tp = db.Column(db.String)
    cm = db.Column(db.String)
    s1 = db.Column(db.String)
    s2 = db.Column(db.String)
    NM = db.Column(db.String)
    MD = db.Column(db.String)
    AS = db.Column(db.String)
    ms = db.Column(db.String)
    nn = db.Column(db.String)
    ts = db.Column(db.String)
    cg = db.Column(db.String)
    cs = db.Column(db.String)
    dv = db.Column(db.String)
    de = db.Column(db.String)
    rl = db.Column(db.String)

    # Not in spec
    SA = db.Column(db.String)

seqrun_schema = {
        "name": generic_string,
        "notes": generic_string,
        "run_id": generic_string,
        "machine_id": generic_string,
        "sequencing_type": {"type": "string", "enum": ['nanopore','illumina']},
        "machine": {"type": "string", "enum": ['iSeq100']},
        "provider": generic_string,
        "well": uuid_schema}
seqrun_required = ['name','run_id','machine_id','sequencing_type','machine','provider']
class SeqRun(db.Model):
    validator = schema_generator(seqrun_schema,seqrun_required)
    put_validator = schema_generator(seqrun_schema,[])

    __tablename__ = 'seqruns'
    id = db.Column(db.Integer, primary_key=True)
    time_created = db.Column(db.DateTime(timezone=True), server_default=func.now())
    time_updated = db.Column(db.DateTime(timezone=True), onupdate=func.now())

    name = db.Column(db.String) # run name
    notes = db.Column(db.String)
    run_id = db.Column(db.String) # Sequencing provider id
    machine_id = db.Column(db.String)
    sequencing_type = db.Column(db.String) # illumina, nanopore, etc
    machine = db.Column(db.String) # minion, iseq, etc
    provider = db.Column(db.String) # in-house
    well = db.Column(db.String)

    fastqs = relationship('Fastq',backref='seqrun')

    def toJSON(self,full=None):
        dictionary = {**base_dict(self), **{'name':self.name,'notes':self.notes,'run_id':self.run_id,'machine_id':self.machine_id,'sequencing_type':self.sequencing_type,'machine':self.machine,'provider':self.provider}}
        if full=='full':
            pass
        return dictionary

class PileupFile(db.Model):
    __tablename__ = 'pileupfiles'
    id = db.Column(db.Integer, primary_key=True)
    time_created = db.Column(db.DateTime(timezone=True), server_default=func.now())

    status = db.Column(db.String) # mutation,confirmed,etc
    full_search_sequence = db.Column(db.String)
    target_sequence = db.Column(db.String)
    index_for = db.Column(db.String)
    index_rev = db.Column(db.String)

    sample_uuid = db.Column(db.String)
    seqrun_id = db.Column(db.Integer, db.ForeignKey('seqruns.id'), nullable=False)

class PileupLine(db.Model):
    __tablename__ = 'pileups'
    id = db.Column(db.Integer, primary_key=True)
    time_created = db.Column(db.DateTime(timezone=True), server_default=func.now())

    # https://en.wikipedia.org/wiki/Pileup_format
    pileup_id = db.Column(db.Integer, db.ForeignKey('pileups.id'), nullable=False)
    sequence = db.Column(db.String)
    position = db.Column(db.Integer)
    reference_base = db.Column(db.String)
    read_count = db.Column(db.Integer)
    read_results = db.Column(db.String)
    quality = db.Column(db.String)

class Job(db.Model):
    __tablename__ = 'jobs'
    id = db.Column(db.Integer, primary_key=True)
    time_created = db.Column(db.DateTime(timezone=True), server_default=func.now())

    job_id = db.Column(db.String)
    status = db.Column(db.String)
    description = db.Column(db.String)

