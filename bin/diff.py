# Compare consensus genome fasta to reference

import argparse
import gzip
import sys

def getfp(filename):
	if   filename.endswith('.gz'):
		return gzip.open(filename, 'rt', encoding='ISO-8859-1')
	elif filename == '-':
		return sys.stdin
	else:
		return open(filename)

def readfasta(filename):
	name = None
	seqs = []
	fp = getfp(filename)
	while True:
		line = fp.readline()
		if line == '': break
		line = line.rstrip()
		if line.startswith('>'):
			if len(seqs) > 0:
				yield name, ''.join(seqs)
				name = line[1:]
				seqs = []
			else:
				name = line[1:]
		else:
			seqs.append(line)
	yield name, ''.join(seqs)
	fp.close()

# Argparse
parser = argparse.ArgumentParser(description="Compare reference and query consensus genomes base-by-base.")
parser.add_argument("-r", "--reference", required=True,
	help="Reference FASTA file")
parser.add_argument("-q", "--query", required=True,
	help="Query Consensus FASTA file")
parser.add_argument("-o", "--output", default="diff.fa",
	help="Output file name [%(default)s]")
parser.add_argument("-w", "--wrap", type=int, default=60,
	help="Line wrapping width [%(default)i]")
args = parser.parse_args()

# Load seqs
ref_seqs = {name: seq.upper() for name, seq in readfasta(args.reference)}
con_seqs = {name: seq.upper() for name, seq in readfasta(args.query)}

# Match consensus sequences to reference using accession prefix
def get_accession(name):
	return name.strip().split()[0]

# Build mapping from accession to full ref name
ref_map = {get_accession(name): name for name in ref_seqs}
con_map = {get_accession(name): name for name in con_seqs}

# Compare sequences and write results
with open(args.output, "w") as out:
	for acc in sorted(con_map.keys()):
		if acc not in ref_map:
			out.write(f"# Skipped: consensus sequence {acc} has no match in reference\n")
			continue
		refname = ref_map[acc]
		conname = con_map[acc]
		refseq = ref_seqs[refname]
		conseq = con_seqs[conname]

		if len(refseq) != len(conseq):
			reflen = len(refseq)
			conlen = len(conseq)
			out.write(f"# Skipped: sequence length mismatch for {acc} (ref={reflen}, query={conlen}, diff={reflen-conlen})\n")
			continue
		
		out.write(f">{refname}\n")
		for i in range(0, len(refseq), args.wrap):
			refseg = refseq[i:i+args.wrap]
			conseg = conseq[i:i+args.wrap]
			matchseg = ''.join(['|' if r == c else ' ' for r, c in zip(refseg, conseg)])
			out.write(f"Sbjct:\t{refseg}\n       \t{matchseg}\nQuery:\t{conseg}\n\n")
