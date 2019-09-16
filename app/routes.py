from kg_flask_crud import create_crud,requires_auth
import pysam
from .models import *
import pandas
import subprocess

from .jobs import fasta_to_db, seqrun_sam_generation 
from .config import *
from flask_restplus import Api, Resource, fields, Namespace
from flask import Flask, abort, request, jsonify, g, url_for, redirect

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


@ns_seqrun.route('/generate_sam/<uuid>')
class GenerateSam(Resource):
    @ns_seqrun.doc('generate_sam',security='token')
    @requires_auth(['moderator','admin'])
    def get(self,uuid):
        new_sam = threading.Thread(target=seqrun_sam_generation, args=(URL,uuid))
        new_sam.start()
        return jsonify({'message':'Processing started'})

@ns_sample.route('/get_sam/<uuid>')
class GetSam(Resource):
    @ns_sample.doc('get_sam',security='token')
    @requires_auth(['moderator','admin'])
    def get(self,uuid):
        obj = Sample.query.filter_by(uuid=uuid).first()
        sql_query = """
        SELECT DISTINCT f.uuid as qname, s.flag as flag, sp.sample_uuid as rname, s.pos as pos, s.mapq as mapq, s.cigar as cigar, s.rnext as rnext, s.pnext as pnext, s.tlen as tlen, f.sequence as seq, f.read_quality as qual, s.nm as nm
        FROM "samples" as sp
        JOIN "sams" as s ON s.sample_uuid=sp.uuid
        JOIN "fastqs" as f on f.uuid=s.fastq_uuid
        WHERE sp.uuid='{}'
        """.format(uuid)

        result = db.engine.execute(sql_query)

        header = { 'HD': {'VN': '1.6'},
            'SQ': [{'LN': 2687, 'SN': '{}'}] }


        with pysam.AlignmentFile('/dev/shm/seq/pysam.bam', "wb", header=header) as outf:
            for s in result:
                a = pysam.AlignedSegment()
                a.query_name = str(s[0])
                a.query_sequence= s[9]
                a.flag = s[1]
                a.reference_id = 0
                a.reference_start = s[3]
                a.mapping_quality = s[4]
                a.cigarstring = s[5]
                a.next_reference_start= s[7]
                a.template_length= s[8]
                a.query_qualities = pysam.qualitystring_to_array(s[10])
                outf.write(a)


        
        #with open("/dev/shm/seq/new_sam_out.sam", "w") as f:
        #    f.write("@HD     VN:1.6  SO:coordinate\n@SQ     SN:test LN:2687\n" + open("/dev/shm/seq/sam_out.sam").read()) # mpileup requires header line (@HD) with a version number (1.5 because biostars example)


        with open('/dev/shm/seq/sam_tmp.fa', 'w') as f:
            f.write('>{}\n'.format(obj.sample_uuid))
            f.write(obj.full_seq)
        print(len(obj.full_seq))
        sam_file = subprocess.check_output(" samtools mpileup -f /dev/shm/seq/sam_tmp.fa /dev/shm/seq/sam_out.sam -o /dev/shm/seq/pileup_out.pileup",shell=True).decode("utf-8").rstrip()#.split('\n')
        
        pileup = pandas.read_csv('/dev/shm/seq/pileup_out.pileup',sep='\t')
        print(pileup)




        return jsonify({'message':'complete'})

ns = [ns_seqrun,ns_sample,ns_file]


