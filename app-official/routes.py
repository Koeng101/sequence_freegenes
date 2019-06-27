from kg_flask_crud import create_crud,requires_auth
from .models import *

ns_seqrun = create_crud('seqrun','Seqruns',Seqrun)
ns_pileup = create_crud('pileup','Pileups',Pileup)

ns = [ns_seqrun,ns_pileup]

