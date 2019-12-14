from __future__ import print_function
import sys

k = int(sys.argv[3])

# from https://github.com/lh3/readfq/blob/master/readfq.py
def readfq(fp): # this is a generator function
    last = None # this is a buffer keeping the last unprocessed line
    while True: # mimic closure; is it a bad idea?
        if not last: # the first record or a record following a fastq
            for l in fp: # search for the start of the next record
                if l[0] in '>@': # fasta/q header line
                    last = l[:-1] # save this line
                    break
        if not last: break
        name, seqs, last = last[1:].partition(" ")[0], [], None
        for l in fp: # read the sequence
            if l[0] in '@+>':
                last = l[:-1]
                break
            seqs.append(l[:-1])
        if not last or last[0] != '+': # this is a fasta record
            yield name, ''.join(seqs), None # yield a fasta record
            if not last: break
        else: # this is a fastq record
            seq, leng, seqs = ''.join(seqs), 0, []
            for l in fp: # read the quality
                seqs.append(l[:-1])
                leng += len(l) - 1
                if leng >= len(seq): # have read enough quality
                    last = None
                    yield name, seq, ''.join(seqs); # yield a fastq record
                    break
            if last: # reach EOF before reading enough quality
                yield name, seq, None # yield a fasta record instead
                break

def get_cannocalization(kmer):
    rev_comp = {'A': 'T',
                'T': 'A',
                'C': 'G',
                'G': 'C',
                'a': 't',
                't': 'a',
                'c': 'g',
                'g': 'c',
                'n': 'n',
                'N': 'N'}
    rev_kmer = ''.join([ rev_comp[nt] for nt in kmer  ][::-1])
    if kmer > rev_kmer:
        return rev_kmer
    else:
        return kmer


#k = 31
mode = "full"
#mode = "3prime"
cano = False
level = "gene"
#level = "txp"

def main():
    from collections import defaultdict
    import bisect
    import pandas as pd

    print("kmer size: ", k, file=sys.stderr)
    print("Mode: ", mode, file=sys.stderr)
    print("Canonicalize: ", cano, file=sys.stderr)
    print("Level: ", level, file=sys.stderr)

    tgmap = {}
    for l in open(sys.argv[2]):
        toks = l.rstrip().split('\t')
        tgmap[toks[0].split("|")[0]] = toks[1]

    kdict = defaultdict(list)
    gene_ids = {}
    ufrac = {}
    gene_names = {}
    tfile = open(sys.argv[1])
    n, slen, qlen = 0, 0, 0
    for name, seq, qual in readfq(tfile):
        if n % 1000 == 0:
            print("\rprocessed {} transcripts".format(n), file=sys.stderr, end="")
        toks = name.split(' ')
        if level == "gene":
            gname = tgmap[toks[0].split("|")[0]]
        else:
            gname = toks[0].split("|")[0]

        if gname not in gene_ids:
            gid = len(gene_names)
            gene_ids[gname] = gid
            gene_names[gid] = gname
            ufrac[gname] = [0, 0]
        gnum = gene_ids[gname]
        n += 1
        slen += len(seq)

        start = 0
        seq_len = len(seq)
        if mode == "3prime":
            if seq_len > 1000:
                start = seq_len - 1000
        elif mode != "full":
            print ("ERROR")
            return 0

        for s in range(start, len(seq)-k+1):
            if cano:
                km = get_cannocalization( seq[s:s+k] )
            km = seq[s:s+k]
            l = kdict[km]
            if len(l) == 0:
                kdict[km] = [gnum]
            else:
                position = bisect.bisect(l, gnum)
                if position == 0 or (position > 0 and l[position - 1] != gnum):
                    bisect.insort(l, gnum)
                    kdict[km] = l

    nkmer = len(kdict)
    n = 0
    for km, glist in kdict.items():
        if n % 10000 == 0:
            print("\rprocessed {} of {} kmers".format(n, nkmer), file=sys.stderr, end="")
        n += 1
        is_unique = len(glist) == 1
        gnames = [gene_names[g] for g in glist]
        for gname in gnames:
            ufrac[gname][0] += 1 if is_unique else 0
            ufrac[gname][1] += 1

    print("gene\tunique\ttotal")
    for g, v in ufrac.items():
        print("{}\t{}\t{}".format(g, v[0], v[1]))

if __name__ == "__main__":
    main()
