#!/bin/env python

from Bio import Entrez
from Bio import SeqIO
Entrez.email = "nedelman@g.harvard.edu"
handle = Entrez.efetch(db="nucleotide", id="101738423")
handle.close()
