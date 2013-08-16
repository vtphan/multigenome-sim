'''
Generate N reads with specified error rate from a given genome.

Usage: python generate_reads.py genome_file.fasta length N [-e ERROR_RATE]
'''
import sys
import argparse
import random
import math

ERR_RATE = 0.02
DEBUG = False

BASES = ['A','C','G','T']

#------------------------------------------------------------------------------
def generate_errors(dna):
   genome = list(dna)
   errors = int(ERR_RATE * len(genome))
   for i in range(errors):
      pos = random.randint(0, len(genome)-1)
      choices = list(BASES)
      choices.remove(genome[pos].upper())
      c = random.choice(choices)
      genome[pos] = c if not DEBUG else c.lower()
   return genome

#------------------------------------------------------------------------------
def generate_reads(genome, header, length, N):
   dna = ''.join(genome)
   with open('reads.fasta', 'w') as f:
      for idx in range(N):
         pos = random.randint(0, len(dna)-length)
         new_header = '>r%s %s:%s %s' % (idx+1, pos, pos+length+1, header[1:])
         f.write(new_header)
         f.write('\n')
         f.write(dna[pos: pos+length+1])
         f.write('\n')

   print 'Save output to reads.fasta'

#------------------------------------------------------------------------------
if __name__ == '__main__':
   parser = argparse.ArgumentParser(description='Generate reads based on a reference genome.')
   parser.add_argument('file_name', help='genome_file.fasta')
   parser.add_argument('length', type=int, help='read length')
   parser.add_argument('N', type=int, help='number of reads')
   parser.add_argument('-e','--error_rate', help='Base error rate (default %s)' % ERR_RATE, type=float, default=ERR_RATE)
   parser.add_argument('--debug', help='Turn on debug mode (default %s)' % DEBUG, action="store_true")

   args = vars(parser.parse_args())
   ERR_RATE = args['error_rate']
   file_content = open(args['file_name']).read().strip()
   N = args['N']
   length = args['length']
   DEBUG = args['debug']

   lines = file_content.splitlines()
   if lines[0][0] != '>':
      print 'Invalid fasta format.  First line must be a description'
      sys.exit(0)
   header = lines.pop(0)
   dna = ''.join(lines).upper()

   genome = generate_errors(dna)
   generate_reads(genome, header, length, N)