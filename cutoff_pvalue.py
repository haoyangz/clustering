from __future__ import division
from scipy.special import comb, perm
from operator import mul
import numpy as np, argparse, string, Levenshtein as lev

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-L', type=int, help='the length of the sequences')
    parser.add_argument('-M', type=int, help='the number of alternatives per base')
    parser.add_argument('-c', type=int, help='the distance cutoff')
    parser.add_argument('-n', type=int, help='the num of sequences')
    parser.add_argument('-t', type=int, default=10000, help='the num of simulation trials')
    parser.add_argument('--mode', default='bound', help='analytical bound or simulation')
    return parser.parse_args()


def bound(n, L, M, N, cutoff):
    # The probability that any pair of sequence has a distance larger than cutoff
    P = 1 - sum([  comb(L, k) * (1/M)**(L-k) * (1-1/M)**k  for k in range(0, cutoff+1)])

    # Recursively calculate a upper bound on the probablity that any pair is far enough
    t_J = 1
    J = 1
    for i in range(n-1):
        t_J *= P
        J *= t_J

    # What we need is a conditional probability on all sequences are different
    p_alldifferent = np.exp(sum([np.log(float(x)) for x in range(N-n+1, N+1)]) - n*np.log(float(N)))
    return 1 - J / p_alldifferent


def simulate(n, L, M, N, n_trial, cutoff):
    dict = list(string.ascii_lowercase)[:M]
    result = []
    for trial in range(n_trial):
        # Sample n different sequences
	samples = []
        for idx in range(n):
            flag = False
            while not flag:
                sample = []
                for l in range(L):
                    sample.append(dict[np.random.randint(M)])
                sample = ''.join(sample)
                if sample not in samples:
                    samples.append(sample)
                    flag = True

		# Check if any pair has a distance within the cutoff
        t_result = 0
        for x in range(n-1):
            for y in range(x+1, n):
                if lev.distance(samples[x], samples[y]) <= cutoff:
                    t_result = 1
                    break
            if t_result == 1:
                break
        result.append(t_result)

    return sum(result) / float(n_trial)

args = parse_args()
if args.mode not in ['bound', 'simulate']:
	raise ValueError('Unrecognized args.mode {}'.format(args.mode))

pvalue = bound(args.n, args.L, args.M, args.M**args.L, args.c) if args.mode == 'bound' else simulate(args.n, args.L, args.M, args.M**args.L, args.t, args.c)
print pvalue
