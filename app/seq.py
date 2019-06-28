import pandas as pd
import os
import shutil
import subprocess

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

def test_print(s):
    return s
