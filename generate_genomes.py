import sys
import argparse
import random
import math

SNP = [
   ('A','C'), ('A','G'),('A','T'),('C','G'),('C','T'),('G','T'),
   ('A','C','G'), ('A','C','T'), ('A','G','T'),('C','G','T'),
   ('A','C','G','T'),
]
BASES = ['A','C','G','T']

DEBUG = False
N = 10  # Number of genomes

MUT_RATE = 0.001

# compared with SNPs, indels occur at approximately eightfold lower rates (Lunter 2007; Cartwright 2009)
INDEL_FRAC = 1.0/9.0
INDEL_EXT = 0.3

FIRST_PROB = 0.75

#------------------------------------------------------------------------------
def sample_frequencies( snp_profile ):
   f = [ int(FIRST_PROB * N) ]
   limit = N - int(FIRST_PROB * N)
   sum = int(FIRST_PROB * N)
   for i in range(1, len(snp_profile)-1):
      v = random.randint(0, limit)
      sum += v
      limit = N - sum
      f.append(sum)
   f.append( N )
   return f

#------------------------------------------------------------------------------
def sample_positions(dna, rate):
   num = int( rate * len(dna) )
   pos = set()
   while len(pos) < num:
      x = random.randint(0, len(dna)-1)
      if x not in pos:
         pos.add(x)
   return pos

#------------------------------------------------------------------------------

def build_SNP_profile(dna):
   pos = sample_positions(dna, MUT_RATE)
   SNP_profile = []
   for p in pos:
      snp = random.choice(SNP)
      SNP_profile.append( (p, snp, sample_frequencies(snp)) )

   return SNP_profile

#------------------------------------------------------------------------------

def random_insertion():
   seq = random.choice(BASES)
   while random.random() <= INDEL_EXT:
      seq += random.choice(BASES)
   return seq

def random_deletion_length():
   cur_len = 1
   while random.random() <= INDEL_EXT:
      cur_len += 1
   return cur_len

#------------------------------------------------------------------------------

# def generate_errors( genome ):
#    errors = int(ERR_RATE * len(genome))
#    for i in range(errors):
#       pos = random.randint(0, len(genome)-1)
#       choices = list(BASES)
#       choices.remove(genome[pos].upper())
#       c = random.choice(choices)
#       genome[pos] = c if not DEBUG else c.lower()
#    return genome

#------------------------------------------------------------------------------

def generate_genome(dna, SNP_profile, header, idx):
   genome = list(dna)
   inserts = []
   deletes = []
   for pos, snp, freq in SNP_profile:
      r = random.random()
      if r < INDEL_FRAC:
         if random.random() <= 0.5:
            inserts.append(pos)
         else:
            deletes.append(pos)
      else:
         r2 = random.randint(0, N)
         for i,f in enumerate(freq):
            if r2 <= f:
               genome[pos] = snp[i] if not DEBUG else snp[i].lower()
               break

   for pos in deletes:
      p = pos
      for i in range(random_deletion_length()):
         if p < len(genome):
            genome[p] = '-'

   for pos in inserts:
      genome[pos] += (random_insertion() if not DEBUG else random_insertion().lower())
   genome = [ x for x in genome if x != '-']

   genome = ''.join(genome)

   new_header = '>%s.%s' % (idx+1, header[1:])

   if idx == 0:
      save_to_file('x.fasta', 'w', new_header, genome, idx)
   else:
      if idx == 1:
         save_to_file('reference.fasta', 'w', new_header, genome, idx)
      save_to_file('multigenome.fasta', 'a', new_header, genome, idx)

#------------------------------------------------------------------------------
def save_to_file(filename, mode, header, genome, idx):
   with open(filename, mode) as f:
      f.write(header)
      f.write('\n')
      f.write(genome)
      f.write('\n')
      # print 'Save genome #%s to %s' % (idx+1, filename)

#------------------------------------------------------------------------------
def log2(x):
   return math.log(x)/math.log(2)

#------------------------------------------------------------------------------
def entropy( freq ):
   N = float(freq[-1])
   p = [float(freq[0])/N]
   for i in range(1, len(freq)):
      p.append( float(freq[i]-freq[i-1]) / N )
   e = 0
   count = 0
   for x in p:
      if x > 0:
         e += - (x * log2(x))
         count += 1
   return e, log2(count)

#------------------------------------------------------------------------------

def output_SNP_profile( profile ):
   with open('snp_profile.txt', 'w') as f:
      if profile:
         average_e = 0.0

         f.write('Note: at each position, an indel occurs with default probability of 11%.\n')
         f.write('Pos\tSNP\tE\tMax E\tCumulative Expected Frequency\n')
         for pos, snp, freq in profile:
            entr, max_entr = entropy(freq)
            average_e += entr
            f.write('%s\t%s\t%.4f\t%.4f\t%s\n' %
               (pos, ''.join(snp), entr, max_entr, ','.join(str(i) for i in freq)))

         f.write('Average entropy: %.4f\n' %  (average_e / float(len(profile))))
      print 'Save SNP profile to snp_profile.txt'

#------------------------------------------------------------------------------
if __name__ == '__main__':
   parser = argparse.ArgumentParser(description='Generate multiple genomes based on a reference genome.')
   parser.add_argument('file_name', help='genome_file.fasta')
   parser.add_argument('-m','--mutation_rate',
      help='Base mutation rate (default %s)' % MUT_RATE, type=float, default=MUT_RATE)
   parser.add_argument('-p', '--first_prob', help='Probability of first base in a SNP profile (default %s)'%FIRST_PROB, type=float, default=FIRST_PROB)
   parser.add_argument('-i','--indel_frac',
      help='Indel fraction of SNP (default %.4f)' % INDEL_FRAC, type=float, default=INDEL_FRAC)
   parser.add_argument('-ie','--indel_ext',
      help='Indel extension probability (default %s)' % INDEL_EXT, type=float, default=INDEL_EXT)
   parser.add_argument('-n', help='Number of genomes (default %s)'%N, type=int, default=N)
   parser.add_argument('--debug', help='Turn on debug mode (default %s)' % DEBUG, action="store_true")

   args = vars(parser.parse_args())
   MUT_RATE, INDEL_FRAC, INDEL_EXT = args['mutation_rate'],args['indel_frac'],args['indel_ext']
   N = args['n']
   FIRST_PROB = args['first_prob']
   DEBUG = args['debug']

   file_content = open(args['file_name']).read().strip()
   lines = file_content.splitlines()
   if lines[0][0] != '>':
      print 'Invalid fasta format.  First line must be a description'
      sys.exit(0)
   header = lines.pop(0)
   dna = ''.join(lines).upper()
   if any( x for x in dna if x not in ['A','C','G','T','N'] ):
      print 'Content is not a valid DNA sequence, which must consist of only a,c,g,t.'
      sys.exit(0)

   #-------------------------------------
   SNP_profile = build_SNP_profile(dna)
   output_SNP_profile(SNP_profile)

   # clear old result
   with open('multigenome.fasta','w') as f:
      f.write('')

   for i in range(N):
      generate_genome(dna, SNP_profile, header, i)

   print 'Genomes are saved in x.fasta, reference.fasta and multigenome.fasta\n'

#------------------------------------------------------------------------------
# def w0(x): # Lambert W function using Newton's method
#    eps = 0.0000001 # max error allowed
#    w = x
#    while True:
#       print '>', w
#       ew = math.exp(w)
#       wNew = w - (w * ew - x) / (w * ew + ew)
#       if abs(w - wNew) <= eps: break
#       w = wNew
#    return w

# #------------------------------------------------------------------------------
# def find_prob( entropy ):
#    return math.e ** w0(entropy)

# #------------------------------------------------------------------------------
# # find p such that -k*p*log(p) = E, or p*log(p) = - E/k
# #------------------------------------------------------------------------------
# def find_p( entropy, k ):
#    return math.e ** w0( - float(entropy) / float(k) )
