from kg_flask_crud import create_crud,requires_auth
from .models import *

ns_seqrun = create_crud('seqrun','SeqRuns',SeqRun,db)
#ns_pileup = create_crud('pileup','Pileups',PileupFile)

ns = [ns_seqrun]

