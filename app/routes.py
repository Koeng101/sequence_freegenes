from kg_flask_crud import create_crud,requires_auth
from .models import *
from redis import Redis
from rq.job import Job
from .jobs import fasta_to_db
from .config import *
from flask_restplus import Api, Resource, fields, Namespace
from flask import Flask, abort, request, jsonify, g, url_for, redirect

import threading

ns_seqrun = create_crud('seqrun','SeqRuns',SeqRun,db)
#ns_pileup = create_crud('pileup','Pileups',PileupFile)


ns_file = Namespace('fastq_upload', description='Fastq.gz file upload')

@ns_file.route('/')
class NewFile(Resource):
    @ns_file.doc('new_file',security='token')
    @requires_auth(['moderator','admin'])
    def post(self):
        json_file = json.loads(request.files['json'].read())
        file = request.files['file']
        
        file_name = json_file['file_name']
        seqrun_id = json_file['seqrun_id']
        with gzip.open(file,'rt') as fin:
            full_file = fin.readlines()
            new_job = threading.Thread(target=fasta_to_db, args=(URL, full_file,file_name,seqrun_id))
            new_job.start()
        return jsonify({'message': 'Successful, job queued'})



ns = [ns_seqrun,ns_file]

