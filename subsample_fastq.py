import sys
import random
import HTSeq

# Usage: python subsample.py <fraction> <input file 1> <input file 2> <output file 1> <output file 2>
# Where <fraction> is a number between 0 and 1, giving the sampling faction

fraction = float(sys.argv[1])
in1 = iter(HTSeq.FastqReader(sys.argv[2]))
in2 = iter(HTSeq.FastqReader(sys.argv[3]))
out1 = open(sys.argv[4], "w")
out2 = open(sys.argv[5], "w")

while True:
    read1 = next(in1)
    read2 = next(in2)
    if random.random() < fraction:
        read1.write_to_fastq_file(out1)
        read2.write_to_fastq_file(out2)

out1.close()
out2.close()
