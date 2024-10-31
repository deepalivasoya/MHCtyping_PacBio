import sys
from Bio import AlignIO
from Bio.Align import AlignInfo
from Bio.Seq import Seq

def consensus(aln=sys.argv[1], gapped=True, threshold=0.5):
    with open(aln) as iH:
        aln_handle = AlignIO.read(iH, "clustal")

    summary_info = AlignInfo.SummaryInfo(aln_handle)
    return str(summary_info.gap_consensus(threshold=0.1, require_multiple=1))

if __name__ == "__main__":
    consensus()

f = open(sys.argv[2], "w")
check=consensus(sys.argv[1])
f.write (">")
f.write (sys.argv[3])
f.write ("\n")
f.write (check)
f.write ("\n")
