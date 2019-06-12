from sqlalchemy.dialects.postgresql import UUID
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
def get_total_bytes(s3, key):
    result = s3.list_objects(Bucket=BUCKET)
    for item in result['Contents']:
        if item['Key'] == key:
            return item['Size']

def get_object(s3, total_bytes,key):
    if total_bytes > 1000000:
        return get_object_range(s3, total_bytes, key)
    return s3.get_object(Bucket=BUCKET, Key=key)['Body'].read()

def get_object_range(s3, total_bytes, key):
    offset = 0
    while total_bytes > 0:
        end = offset + 999999 if total_bytes > 1000000 else ""
        total_bytes -= 1000000
        byte_range = 'bytes={offset}-{end}'.format(offset=offset, end=end)
        offset = end + 1 if not isinstance(end, str) else None
        yield s3.get_object(Bucket=BUCKET, Key=key, Range=byte_range)['Body'].read()

class Files(db.Model):
    def __init__(self,name,file,plate_type,order_uuid,status,plate_name,breadcrumb,plate_vendor_id):
        print(name)
        file_name = str(uuid.uuid4())
        def upload_file_to_spaces(file,file_name=file_name,bucket_name=BUCKET,spaces=SPACES):
            """
            Docs: http://boto3.readthedocs.io/en/latest/guide/s3.html
            http://zabana.me/notes/upload-files-amazon-s3-flask.html"""
            try:
                print('Attempting')
                spaces.upload_fileobj(file,bucket_name,file_name)
                print('Uploaded')
            except Exception as e:
                print("Failed: {}".format(e))
                return False
            return True
        if upload_file_to_spaces(file,file_name=file_name) == True:
            self.name = name
            self.file_name = file_name
            self.plate_type = plate_type
            self.order_uuid = order_uuid
            self.status = status
            self.breadcrumb = breadcrumb
            self.plate_name = plate_name
            self.plate_vendor_id = plate_vendor_id
            # Also include plate_name in json_file if uploading a glycerol stock or dna stock
    __tablename__ = 'files'
    uuid = db.Column(UUID(as_uuid=True), unique=True, nullable=False,default=sqlalchemy.text("uuid_generate_v4()"), primary_key=True)
    time_created = db.Column(db.DateTime(timezone=True), server_default=func.now())

    name = db.Column(db.String, nullable=False) # Name to be displayed to user
    file_name = db.Column(db.String, nullable=False) # Link to spaces

    
    def toJSON(self,full=None):
        return {'uuid':self.uuid,'name':self.name,'file_name':self.file_name,'plate_type':self.plate_type,'order_uuid':self.order_uuid,'breadcrumb':self.breadcrumb,'plate_name':self.plate_name,'status':self.status,'plate_vendor_id':self.plate_vendor_id}
    def download(self):
        s3 = SPACES
        key = self.file_name
        total_bytes = get_total_bytes(s3,key)
        return Response(
            get_object(s3, total_bytes, key),
            mimetype='text/plain',
            headers={"Content-Disposition": "attachment;filename={}".format(self.name)})

seqrun_schema = {
        "uuid": uuid_schema,
        "name": generic_string,
        "notes": generic_string,
        "basecaller": generic_string,
        "basecaller_version_number": generic_string,
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
    basecaller = db.Column(db.String)
    basecaller_version_number = db.Column(db.String)
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

pileup_fastq = db.Table('pileup_fastq',
    db.Column('pileup_uuid', UUID(as_uuid=True), db.ForeignKey('pileups.uuid'), primary_key=True),
    db.Column('fastq_uuid', UUID(as_uuid=True), db.ForeignKey('fastqs.uuid'),primary_key=True,nullable=True),
)

pileup_schema = {
        "uuid": uuid_schema,
        "status": {"type": "string", "enum": ["Mutated","Confirmed","Failed"]},
        "full_search_sequence": dna_string,
        "target_sequence": dna_string,
        "sample_uuid": uuid_schema,
        "fastqs": {"type": "array", "items": [uuid_schema]},
        "file_uuid": uuid_schema
        }
pileup_required = ['full_search_sequence','target_sequence','sample_uuid']
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

    sample_uuid = db.Column(db.String())
    fastqs = db.relationship('Fastq', secondary=pileup_fastq, lazy='subquery',
        backref=db.backref('pileups', lazy=True))

    file_uuid = db.Column(UUID, db.ForeignKey('files.uuid'),nullable=True)
    
    def toJSON(self,full=None):
        dictionary= {'uuid':self.uuid,'time_created':self.time_created,'time_updated':self.time_updated,'status':self.status,'full_search_sequence':self.full_search_sequence,'target_sequence':self.target_sequence,'file_uuid':self.file_uuid, 'sample_uuid': self.sample_uuid}
        if full=='full':
            dictionary['fastqs'] = [fastq.uuid for fastq in self.fastqs]
        return dictionary

fastq_schema = {
        "uuid": uuid_schema,
        "seqrun_uuid": uuid_schema,
        "file_uuid": uuid_schema,
        "name": generic_string,
        "index_for": dna_string,
        "index_rev": dna_string,
        }
fastq_required = ['seqrun_uuid','file_uuid','index_for','index_rev']
class Fastq(db.Model):
    validator = schema_generator(fastq_schema,fastq_required)
    put_validator = schema_generator(fastq_schema,[])

    __tablename__ = 'fastqs'
    uuid = db.Column(UUID(as_uuid=True), unique=True, nullable=False,default=sqlalchemy.text("uuid_generate_v4()"), primary_key=True)
    time_created = db.Column(db.DateTime(timezone=True), server_default=func.now())
    time_updated = db.Column(db.DateTime(timezone=True), onupdate=func.now())
    seqrun_uuid = db.Column(UUID, db.ForeignKey('seqruns.uuid'), nullable=False)
    name = db.Column(db.String)

    file_uuid = db.Column(UUID, db.ForeignKey('files.uuid'),nullable=False)

    index_for = db.Column(db.String)
    index_rev = db.Column(db.String)
    
    def toJSON(self,full=None):
        dictionary= {'uuid':self.uuid,'time_created':self.time_created,'time_updated':self.time_updated,'seqrun_uuid':self.seqrun_uuid,'file_uuid':self.file_uuid,'index_for':self.index_for,'index_rev':self.index_rev}
        if full=='full':
            dictionary['pileups'] = [pileup.uuid for pileup in self.pileups]
        return dictionary



