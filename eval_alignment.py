import sys
import argparse
import random
import math

if __name__ == '__main__':
   parser = argparse.ArgumentParser(description='Evaludate aligned reads.')
   parser.add_argument('reads', help='file containing all reads.')
   parser.add_argument('alignment', help='file containing aligned reads.')
   parser.add_argument('-g', help='Maximally allowable gap distance between positions of a read and aligned read.  Default=20', default=20)

   args = vars(parser.parse_args())
   read_file, alignment_file, G = args['reads'],args['alignment'],int(args['g'])
   read_count = 0
   pos_count = 0
   reads = {}
   for line in open(read_file):
      line = line.strip()
      if '\t' not in line:
         continue
      read, pos = line.split('\t')
      positions = set(int(p) for p in pos.split(' '))
      if positions:
         reads[read_count] = (read, set(positions))
         read_count += 1
         pos_count += len(positions)
   # print reads
   # print read_count, pos_count

   aligned_read_count = 0
   aligned_position_count = 0
   TP, FP, FN = 0.0, 0.0, float(pos_count)
   for line in open(alignment_file):
      line = line.strip()
      if '\t' not in line:
         continue
      idx, pos = line.split('\t')
      idx = int(idx)
      if idx not in reads:
         print "Error: non-existent index (%d)." % idx
         sys.exit(0)

      positions = set(int(p) for p in pos.split(' '))
      if positions:
         aligned_read_count += 1
         aligned_position_count += len(positions)
         read, correct_pos = reads[idx]
         tp, fp = 0, 0
         for p in positions:
            missing = True
            for q in correct_pos:
               if abs(p-q) <= G:
                  tp += 1
                  missing = False
                  break
            if missing:
               fp += 1
         # TP += len(correct_pos & positions)
         # FP += len(positions - correct_pos)
         # FN -= len(correct_pos & positions)
         print idx,read
         for p in correct_pos:
            print p,
         print
         for p in positions:
            print p,
         print
         print tp, fp, len(correct_pos)-tp
         TP += tp
         FP += fp
         FN -= tp
   print "G=%d" % G
   print "Reads\tPos\tAreads\tApos\tTP\tFP\tFN\tPrecision\tRecall"
   print "%d\t%d\t%d\t%d\t%d\t%d\t%d\t%.8f\t%.8f\n" % (
      read_count, pos_count, aligned_read_count, aligned_position_count,
      TP, FP, FN, TP/(TP+FP), TP/(TP+FN))

