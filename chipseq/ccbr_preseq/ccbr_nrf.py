from __future__ import print_function
import sys
preseq_log=sys.argv[1]
with open(preseq_log, 'r') as fp:
        for line in fp:
            if line.startswith('TOTAL READS'):
                tot_reads = float(line.strip().split("= ")[1])
            elif line.startswith('DISTINCT READS'):
                distinct_reads = float(line.strip().split('= ')[1])
            elif line.startswith('1\t'):
                one_pair = float(line.strip().split()[1])
            elif line.startswith('2\t'):
                two_pair = float(line.strip().split()[1])
NRF = distinct_reads/tot_reads
PBC1 = one_pair/distinct_reads
PBC2 = one_pair/two_pair
print("NRF:%.3f\nPBC1:%.3f\nPBC2:%.3f"%(NRF,PBC1,PBC2))
