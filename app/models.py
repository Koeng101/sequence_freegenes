ifrom sqlalchemy.dialects.postgresql import UUID
import sqlalchemy
from sqlalchemy.sql import func
from flask_sqlalchemy import SQLAlchemy
from flask_httpauth import HTTPBasicAuth
from flask import Flask, abort, request, jsonify, g, url_for, Response
import uuid


from .config import SPACES
from .config import BUCKET

from itsdangerous import (TimedJSONWebSignatureSerializer
                          as Serializer, BadSignature, SignatureExpired)
from passlib.apps import custom_app_context as pwd_context

db = SQLAlchemy()
auth = HTTPBasicAuth()

##################
### Validators ###
##################

from jsonschema import validate
import json
import string

# Shared
uuid_regex = '^[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}$'
null = {'type': 'null'}

uuid_schema = {'type': 'string','pattern': uuid_regex}
optional_uuid = {'oneOf': [uuid_schema,null]}

generic_string = {'type': 'string'}
optional_string ={'oneOf': [generic_string,null]}

dna_string = {"type": "string", "pattern": "^[ATGC]*$"}

generic_num = { "type": "number" }
optional_num = {'oneOf': [generic_num,null]}

generic_date = {'type': 'string','format':'date-time'}
optional_date = {'oneOf': [generic_date,null]}

name = {'type': 'string','minLength': 3,'maxLength': 30}
tags = {'type': 'array', 'items': optional_string}
force_to_many = {'type': 'array', 'items': uuid_schema}
to_many = {'type': 'array', 'items': {'oneOf': [uuid_schema,null]}}
#many_to_many = {'anyOf': [{'type': 'array','items': uuid},{'type': 'array','items': null}]}

def schema_generator(properties,required,additionalProperties=False):
    return {"$schema": "http://json-schema.org/schema#",
            "type": "object",
            "properties": properties,
            "required": required,
            "additionalProperties": additionalProperties}


#################
### FG import ###
#################

seqrun_schema = {
        "uuid": uuid_schema,
        "name": generic_string,
        "notes": generic_string,
        "machine_id": generic_string,
        "sequencing_type": {"type": "string", "enum": ["Nanopore","Illumina","PacBio"]},
        "machine": {"type": "string", "enum": ["MinION","iSeq"]},
        "provider": {"type": "string", "enum": ["in-house"]},
        "job": generic_string,
        }
seqrun_required = ['name','machine_id','sequencing_type','machine','provider']
class Seqrun(db.Model):
    validator = schema_generator(seqrun_schema,seqrun_required)
    put_validator = schema_generator(seqrun_schema,[])

    __tablename__ = 'seqruns'
    uuid = db.Column(UUID(as_uuid=True), unique=True, nullable=False,default=sqlalchemy.text("uuid_generate_v4()"), primary_key=True)
    time_created = db.Column(db.DateTime(timezone=True), server_default=func.now())
    time_updated = db.Column(db.DateTime(timezone=True), onupdate=func.now())

    name = db.Column(db.String) # run name
    notes = db.Column(db.String)
    run_id = db.Column(db.String) # Sequencing provider id
    machine_id = db.Column(db.String)
    sequencing_type = db.Column(db.String) # illumina, nanopore, etc
    machine = db.Column(db.String) # minion, iseq, etc
    provider = db.Column(db.String) # in-house

    job = db.Column(db.String) # the job id of the redis job

    fastqs = db.relationship('Fastq',backref='seqrun')

    def toJSON(self,full=None):
        dictionary= {'uuid':self.uuid,'time_created':self.time_created,'time_updated':self.time_updated,'name':self.name,'run_id':self.run_id,'machine_id':self.machine_id,'notes':self.notes,'sequencing_type':self.sequencing_type,'machine':self.machine,'provider':self.provider, 'job':self.job}
        if full=='full':
            dictionary['fastqs'] = [fastq.uuid for fastq in self.fastqs]
        return dictionary

pileup_schema = {
        "uuid": uuid_schema,
        "status": {"type": "string", "enum": ["Mutated","Confirmed","Failed"]},
        "full_search_sequence": dna_string,
        "target_sequence": dna_string,
        "sample_uuid": uuid_schema,
        "seqrun_uuid": uuid_schema,
        "index_for": generic_string,
        "index_rev": generic_String,
        }
pileup_required = ['full_search_sequence','target_sequence','sample_uuid','index_for','index_rev']
class Pileup(db.Model):
    validator = schema_generator(pileup_schema,pileup_required)
    put_validator = schema_generator(pileup_schema,[])

    __tablename__ = 'pileups'
    uuid = db.Column(UUID(as_uuid=True), unique=True, nullable=False,default=sqlalchemy.text("uuid_generate_v4()"), primary_key=True)
    time_created = db.Column(db.DateTime(timezone=True), server_default=func.now())
    time_updated = db.Column(db.DateTime(timezone=True), onupdate=func.now())

    status = db.Column(db.String) # mutation,confirmed,etc
    full_search_sequence = db.Column(db.String)
    target_sequence = db.Column(db.String) 
    index_for = db.Column(db.String)
    index_rev = db.Column(db.String)

    sample_uuid = db.Column(db.String())
    seqrun_uuid = db.Column(UUID, db.ForeignKey('seqruns.uuid'), nullable=False)

    
    def toJSON(self,full=None):
        dictionary= {'uuid':self.uuid,'time_created':self.time_created,'time_updated':self.time_updated,'status':self.status,'full_search_sequence':self.full_search_sequence,'target_sequence':self.target_sequence,'file_uuid':self.file_uuid, 'sample_uuid': self.sample_uuid}
        if full=='full':
            dictionary['fastqs'] = [fastq.uuid for fastq in self.fastqs]
        return dictionary

class Fastq(db.Model):
    __tablename__ = 'fastqs'
    uuid = db.Column(UUID(as_uuid=True), unique=True, nullable=False,default=sqlalchemy.text("uuid_generate_v4()"), primary_key=True)
    time_created = db.Column(db.DateTime(timezone=True), server_default=func.now())
    time_updated = db.Column(db.DateTime(timezone=True), onupdate=func.now())

    file_name = db.Column(db.String)
    seqrun_uuid = db.Column(UUID, db.ForeignKey('seqruns.uuid'), nullable=False)

    # Shared
    instrument_id = db.Column(db.String())
    flow_cell_id = db.Column(db.String())
    run_id = db.Column(db.String())

    # Types
    type_instrument = db.Column(db.String())
    type_sequencing = db.Column(db.String())

    # Fastq
    sequence = db.Column(db.String())
    comments = db.Column(db.String())
    read_quality = db.Column(db.String())

    # Nanopore
    read = db.Column(db.String())
    channel = db.Column(db.String())
    start_time = db.Column(db.String())
    basecaller = db.Column(db.String())
    basecaller_version = db.Column(db.String())

    # Illumina
    lane = db.Column(db.Integer())
    tile_number = db.Column(db.Integer())
    x_coord = db.Column(db.Integer())
    y_coord = db.Column(db.Integer())
    member_pair = db.Column(db.Integer())
    read_filter = db.Column(db.String()) # Y or N
    control_bits = db.Column(db.Integer())
    index_for = db.Column(db.String())
    index_rev = db.Column(db.String()) 

class Sam(db.Model):
    __tablename__ = 'sams'
    uuid = db.Column(UUID(as_uuid=True), unique=True, nullable=False,default=sqlalchemy.text("uuid_generate_v4()"), primary_key=True)
    time_created = db.Column(db.DateTime(timezone=True), server_default=func.now())
    time_updated = db.Column(db.DateTime(timezone=True), onupdate=func.now())

    # https://en.wikipedia.org/wiki/SAM_(file_format)
    # https://samtools.github.io/hts-specs/SAMv1.pdf
    qname = db.Column(db.String())
    flag = db.Column(db.Integer())
    rname = db.Column(db.String())
    pos = db.Column(db.Integer())
    mapq = db.Column(db.Integer())
    cigar = db.Column(db.String())
    rnext = db.Column(db.String())
    pnext = db.Column(db.Integer())
    tlen = db.Column(db.Integer())
    seq = db.Column(db.String())
    qual = db.Column(db.String())

class Pileup_lines(db.Model):
    __tablename__ = 'pileup_lines'
    uuid = db.Column(UUID(as_uuid=True), unique=True, nullable=False,default=sqlalchemy.text("uuid_generate_v4()"), primary_key=True)
    time_created = db.Column(db.DateTime(timezone=True), server_default=func.now())
    time_updated = db.Column(db.DateTime(timezone=True), onupdate=func.now())

    # https://en.wikipedia.org/wiki/Pileup_format
    pileup_uuid = db.Column(UUID, db.ForeignKey('pileups.uuid'), nullable=False)
    sequence = db.Column(db.String())
    position = db.Column(db.Integer())
    reference_base = db.Column(db.String())
    read_count = db.Column(db.Integer())
    read_results = db.Column(db.String())
    quality = db.Column(db.String())
