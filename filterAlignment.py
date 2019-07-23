from Bio import AlignIO
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC

#Read alignment from gremlin
align = AlignIO.read('1561723510.fas', 'fasta')
lenInput = len(align[0])
lenAlign = len(align)

keepPos = []
dataPos = []

#The first alignment is dF106_dm, and positions in that sequence are interesting
for i in range(len(align[0])):
    if (not align[0][i] == '-'):
        keepPos.append(i)

#Generate the filtered alignment with only positions in dF106_dm
tempAlign = []
for record in align:
    newSeq = [record.seq[i] for i in keepPos]
    newSeq = "".join(newSeq)
    
    newRecord = SeqRecord(Seq(newSeq)
                          ,id=record.id
                          ,description='')
    tempAlign.append(newRecord)

#Write the alignment
SeqIO.write(tempAlign, "dF106TrimmedAlignment.fasta", "fasta")



