import argparse, numpy as np, cPickle, h5py, tempfile
from os.path import join, dirname, basename, exists
from helper import *
from os import makedirs, system

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', dest='inputfile', type=str, help='a tsv file of sequences and counts')
    parser.add_argument('-o', dest='outdir', type=str, help='a tsv file of sequences and counts')
    parser.add_argument('--ns', dest='n_seed', default=100000, type=int, help='')
    parser.add_argument('-c', dest='cutoff', default=2, type=int, help='')
    parser.add_argument('-j', dest='n_jobs', default=16, type=int, help='')
    parser.add_argument('-b', dest='bs', default=1000, type=int, help='')
    parser.add_argument('--command', type=str, help='')
    parser.add_argument('--target_seq', type=str, help='')
    parser.add_argument('--target_outfile', type=str, help='')
    return parser.parse_args()


args = parse_args()
comp_kwargs = {'compression': 'gzip', 'compression_opts': 1}
seqs_pkl = join(args.outdir, 'seqs.pkl')
raw_seqs_pkl = join(args.outdir, 'raw_seqs.pkl')
assign_prep = join(args.outdir, 'assign_prep.pkl')
assign_after_seedclust = join(args.outdir, 'assign_after_seedclust_cutoff{}.pkl'.format(args.cutoff))
assign_after_nonseedclust = join(args.outdir, 'assign_after_nonseedclust_cutoff{}.pkl'.format(args.cutoff))
dist_seed = join(args.outdir, 'dist_seed.h5')
dist_nonseed = join(args.outdir, 'dist_nonseed.h5')
seed_seqs_file = join(args.outdir, 'seed_seqs.txt')
non_seed_seqs_file= join(args.outdir, 'non_seed_seqs.txt')
dist_seed_batch = join(args.outdir, 'dist_seed_batch')
dist_nonseed_batch = join(args.outdir, 'dist_nonseed_batch')
args.command = args.command.split(',')

if not exists(args.outdir):
    makedirs(args.outdir)

for x in args.command:
    if x not in ['preprocess', 'dist_seed', 'clust_seed', 'dist_nonseed', 'clust_nonseed', 'output', 'extract_dist_seed', 'single_seq_dist', 'find_nb_seed', 'find_nb_nonseed']:
        raise ValueError('args.command \'{}\' not recognized!'.format(x))

if 'preprocess' in args.command:
    # Load the sequences and counts
    seqs = []
    ids = []
    counts = []
    with open(args.inputfile) as f:
        for idx, x in enumerate(f):
            seqs.append(x.split()[0])
	    ids.append(idx)
            counts.append(float(x.split()[1]))
    seqs, ids, counts = np.asarray(seqs), np.asarray(ids), np.asarray(counts)

    # Get the seed (seqs with top counts)
    if len(seqs) > args.n_seed:
        reorder = np.argsort(-counts)
        seed_seqs = seqs[reorder[:args.n_seed]]
        seed_ids = ids[reorder[:args.n_seed]]
        non_seed_seqs = seqs[reorder[args.n_seed:]]
        non_seed_ids = ids[reorder[args.n_seed:]]
    else:
        seed_seqs = seqs
	seed_ids = ids
        non_seed_seqs = []
	non_seed_ids = []

    print 'Num of seed seqs:', len(seed_seqs)
    print 'Num of nonseed seqs:', len(non_seed_seqs)

    # Initialize group and sequence assignment
    group_assign = dict()
    seq_assign = dict()
    for idx in range(len(ids)):
        group_assign[idx] = idx
        seq_assign[idx] = [idx]

    with open(seqs_pkl, 'w') as f:
    	cPickle.dump([seed_seqs, non_seed_seqs, seed_ids, non_seed_ids], f)

    with open(assign_prep, 'w') as f:
        cPickle.dump([group_assign, seq_assign], f)

    with open(raw_seqs_pkl, 'w') as f:
        cPickle.dump([seqs, ids, counts], f)

if 'dist_seed' in args.command:
    with open(seqs_pkl) as f:
	seed_seqs, _, _, _= cPickle.load(f)

    print 'Calculating pairwise distance of the seeds'
    pairwiseDist(dist_seed_batch, seed_seqs, seed_seqs, n_jobs=args.n_jobs, batchsize=args.bs)

if 'dist_nonseed' in args.command:
    with open(seqs_pkl) as f:
	seed_seqs, non_seed_seqs, _, _, = cPickle.load(f)

    print 'Calculating the distance from the nonseeds to the seeds'
    pairwiseDist(dist_nonseed_batch, non_seed_seqs, seed_seqs, n_jobs=args.n_jobs, batchsize=args.bs)

if 'single_seq_dist' in args.command:
    with open(raw_seqs_pkl) as f:
    	seqs, ids, count = cPickle.load(f)

    tmpdir = tempfile.mkdtemp()
    pairwiseDist(tmpdir, seqs, [args.target_seq], n_jobs=args.n_jobs, batchsize=args.bs)
    with open(args.target_outfile, 'w') as fout:
        flag = False
        for idx in range(len([1 for _, _ in enumerate(range(0, len(seqs), args.bs))])):
            with open(join(tmpdir, str(idx)+'.txt')) as f:
                for x in f:
                    if not flag:
                        fout.write(x.strip())
                        flag = True
                    else:
                        fout.write('\t'+x.strip())
    system('rm -r ' + tmpdir)


if 'extract_dist_seed' in args.command:
    with open(seqs_pkl) as f:
        seed_seqs, _, seed_ids, _ = cPickle.load(f)

    target_idx = int(np.where(seed_seqs == args.target_seq)[0])
    batch_num = target_idx / args.bs
    idx_in_batch = target_idx - batch_num*args.bs
    print target_idx, batch_num, idx_in_batch, args.bs
    with open(join(dist_seed_batch, str(batch_num)+'.txt')) as fin, open(args.target_outfile, 'w') as fout:
        for _ in range(idx_in_batch):
            fin.readline()
        fout.write(fin.readline())


if 'find_nb_seed' in args.command:
    with open(seqs_pkl) as f:
    	_, _, seed_ids, _ = cPickle.load(f)

    cut_graph(
       dist_seed_batch,
       seed_ids,
       seed_ids,
       args.cutoff,
       len([1 for _, _ in enumerate(range(0, len(seed_ids), args.bs))]),
       args.n_jobs)

if 'clust_seed' in args.command:
    with open(assign_prep) as f:
    	group_assign, seq_assign = cPickle.load(f)

    with open(seqs_pkl) as f:
    	_, _, seed_ids, _ = cPickle.load(f)

    cutdir = dist_seed_batch + '_cut' + str(args.cutoff)
    print 'Loading nb pickle'
    with open(join(cutdir, 'nb.pkl'), 'rb') as f:
        nb = cPickle.load(f)

    print 'Clustering the seeds'
    merge_cluster_inner(nb, group_assign, seq_assign, [x for x in seed_ids])
    print 'Finish clustering the seeds'

    with open(assign_after_seedclust, 'w') as f:
     	cPickle.dump([group_assign, seq_assign], f)

if 'find_nb_nonseed' in args.command:

    with open(seqs_pkl) as f:
        _, _, seed_ids, non_seed_ids = cPickle.load(f)

    cut_graph(
       dist_nonseed_batch,
       non_seed_ids,
       seed_ids,
       args.cutoff,
       len([1 for _,_ in enumerate(range(0, len(non_seed_ids), args.bs))]),
       args.n_jobs)

if 'clust_nonseed' in args.command:
    with open(assign_after_seedclust) as f:
    	group_assign, seq_assign = cPickle.load(f)

    with open(seqs_pkl) as f:
        _, _, seed_ids, non_seed_ids = cPickle.load(f)


    cutdir = dist_nonseed_batch + '_cut' + str(args.cutoff)
    print 'Loading nb pickle'
    with open(join(cutdir, 'nb.pkl'), 'rb') as f:
        nb = cPickle.load(f)

    print 'Merge the nonseeds with the seeds'
    merge_cluster_outer(nb, group_assign, seq_assign, [x for x in seed_ids], [x for x in non_seed_ids])

    with open(assign_after_nonseedclust, 'w') as f:
    	cPickle.dump([group_assign, seq_assign], f)

if 'output' in args.command:
    with open(assign_after_nonseedclust) as f:
    	group_assign, seq_assign = cPickle.load(f)

    with open(raw_seqs_pkl) as f:
    	seqs, ids, count = cPickle.load(f)

    # Output the clustering result
    with open(join(args.outdir, 'clusters_cut{}'.format(args.cutoff)), 'w') as f:
        f.write('\t'.join(['seq', 'count', 'group', 'group_size'])+'\n')
        for t_group, t_seq_ids in seq_assign.items():
            for seq_id in t_seq_ids:
                f.write('%s\n' % '\t'.join([seqs[seq_id], str(count[seq_id]), str(t_group), str(len(seq_assign[t_group]))]))
    print 'num of groups', len(seq_assign.keys())

