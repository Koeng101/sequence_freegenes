from kg_flask_crud import create_crud,requires_auth
import os
from .models import *
import pandas
import subprocess

from .jobs import fasta_to_db, seqrun_bam_generation 
from .config import *
from flask_restplus import Api, Resource, fields, Namespace
from flask import Flask, abort, request, jsonify, g, url_for, redirect, make_response

import threading

ns_seqrun = create_crud('seqrun','SeqRuns',SeqRun,db)
ns_sample = create_crud('samples','Samples',Sample,db)
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
        seqrun_id = json_file['seqrun_uuid']
        with gzip.open(file,'rt') as fin:
            full_file = fin.readlines()
            new_job = threading.Thread(target=fasta_to_db, args=(URL, full_file,file_name,seqrun_id,json_file['index_for'],json_file['index_rev']))
            new_job.start()
        return jsonify({'message': 'Successful, job queued'})


@ns_seqrun.route('/generate_bam/<uuid>')
class GenerateSam(Resource):
    @ns_seqrun.doc('generate_bam',security='token')
    @requires_auth(['moderator','admin'])
    def get(self,uuid):
        new_sam = threading.Thread(target=seqrun_bam_generation, args=(URL,uuid))
        new_sam.start()
        return jsonify({'message':'Processing started'})

@ns_sample.route('/get_pileup/<uuid>')
class GetPileup(Resource):
    @ns_sample.doc('get_pileup',security='token')
    def get(self,uuid):
        obj = Sample.query.filter_by(uuid=uuid).first()

        with open('{}/bam_tmp.fa'.format(PATH), 'w') as f:
            f.write('>{}\n'.format(obj.sample_uuid))
            f.write(obj.full_seq)
        with open('{}/bam_tmp.bam'.format(PATH),'wb') as f:
            f.write(obj.bam)
        sam_file = subprocess.check_output("samtools mpileup -f {}/bam_tmp.fa {}/bam_tmp.bam -o {}/pileup_out.pileup".format(PATH,PATH,PATH),shell=True).decode("utf-8").rstrip()#.split('\n') 
        pileup = pandas.read_csv('{}/pileup_out.pileup'.format(PATH),sep='\t', names = ["Sequence", "Position", "Reference Base", "Read Count", "Read Results", "Quality"])
        pileup = pileup.set_index('Position')
        os.remove("{}/bam_tmp.fa.fai".format(PATH))


        length = len(obj.search_seq)
        start = obj.full_seq.find(obj.search_seq)
        indexs = start+1, start+length
        gene_df = pileup.loc[indexs[0]:indexs[1]]
        gene_df.loc[:,'Sequence'] = obj.sample_uuid
        gene_df = gene_df.reset_index(drop=True)
        gene_df.index += 1 
        gene_df['Position'] = gene_df.index

        # Format pileup file properly
        resp = make_response(gene_df.to_csv(sep='\t',index=False,columns = ["Sequence", "Position", "Reference Base", "Read Count", "Read Results", "Quality"]))
        resp.headers["Content-Disposition"] = "attachment; filename={}.pileup".format(obj.sample_uuid)
        resp.headers["Content-Type"] = "text/csv"
        return resp

ns = [ns_seqrun,ns_sample,ns_file]


