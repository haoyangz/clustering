import multiprocessing as mp, numpy as np, h5py, collections
from itertools import izip
import Levenshtein as lev, cPickle
from os.path import exists, join
from os import makedirs, system

def dist_slave(args):
    '''The slave process to calculate pairwise distance'''
    group1, group2, outfile = args
    with open(outfile, 'w') as fout:
        for s1 in group1:
            fout.write('%s\n' % '\t'.join(map(str, [mydist(s1, s2) for s2 in group2])))

def pairwiseDist(outdir, group1, group2, n_jobs=10, batchsize=1000):
    '''The master process to calculate pairwise distance'''
    if exists(outdir):
        system('rm -r ' + outdir)
    makedirs(outdir)
    keys = [ [group1[idx:min(len(group1), idx+batchsize)], group2, join(outdir, str(batchidx)+'.txt')]   for batchidx, idx in enumerate(range(0, len(group1), batchsize))]
    pool = mp.Pool(processes=n_jobs)
    pool.map(dist_slave, keys)
    pool.close()
    pool.join()

def merge_cluster_outer(nb, seq2group, group2seq, knownseq, seq2update):
    '''Determine cluster assignment of a group of sequences by projecting to sequences whose group assignment
    has been determined
    '''
    for seq in seq2update:
        # The current cluster assignment
        seq_group = seq2group[seq]
        # The new cluster is the most common cluster of known sequences that
        # are neighbors of the current sequence
        nb_group = [seq2group[known] for known in nb[seq]]
	if nb_group:
            known_group = collections.Counter(nb_group).most_common()[0][0]
            # Add the current sequence to the new cluster
            group2seq[known_group].append(seq)
            # Update the cluster assignment
            seq2group[seq] = known_group
            # Remove the original cluster, assuming there was only one member.
            del group2seq[seq_group]

def merge_cluster_inner(nb, seq2group, group2seq, seq2update):
    '''Cluster the sequences by identifying connected components using DFS'''
    # Status vector to record visited nodes
    status = {}
    for x in seq2update:
        status[x] = 1

    rootidx = 0
    newroot = seq2update[rootidx]
    # DFS update stack
    update_stack = [newroot]
    status[newroot] = 0

    while update_stack:
        to_update = update_stack.pop()
        to_update_group = seq2group[to_update]
        # Search for unvisited neighbors
        for _, point in enumerate(nb[to_update]):
            if status[point]==1:
                point_group = seq2group[point]
                # Merge the two cluster
                group2seq[to_update_group] += group2seq[point_group]
                # Update the cluster assignemnt and status of all the nodes in the neighbor's cluster
                for x in group2seq[point_group]:
                    if status[x] == 1:
                        status[x] = 0
                        seq2group[x] = to_update_group
                        update_stack.append(x)

                # Remove the old cluster
                del group2seq[point_group]


        # If the current DFS is ended, start a new one
        if not update_stack:
            if rootidx < len(seq2update)-1:
                found_new_root = False
                while rootidx < len(seq2update)-1:
                    rootidx += 1
                    if status[seq2update[rootidx]] == 1:
                        newroot = seq2update[rootidx]
                        status[newroot] = 0
                        update_stack.append(newroot)
                        found_new_root = True
                        break
                if not found_new_root:
                    break
            else:
                break

def mydist(s1, s2):
    return lev.distance(s1, s2)

def find_neighbor(args):
    '''The slave process to convert a distance row to a neighbor row'''
    oldfile, newfile, cutoff = args
    out = []
    with open(oldfile) as f, open(newfile, 'w') as fout:
        for x in f:
            out.append(np.argwhere(np.asarray(x.split(), dtype=int) <= cutoff).squeeze())
        fout.write('done')
    return out

def cut_graph(dist_dir, row_ids, col_ids, cutoff, n_file, n_jobs):
    '''The master process to convert a distance matrix (in batches on file) to a neighbor matrix'''
    cutdir = dist_dir + '_cut' + str(cutoff)
    if exists(cutdir):
        system('rm -r ' + cutdir)
    makedirs(cutdir)

    keys = [ [join(dist_dir, str(idx)+'.txt'), join(cutdir, str(idx)+'.txt'), cutoff] for idx in range(n_file)]
    pool = mp.Pool(processes=n_jobs)
    result = pool.map(find_neighbor, keys)
    pool.close()
    pool.join()

    nb = {}
    cnt = 0
    for re in result:
        for x in re:
            if x.shape:
                nb[row_ids[cnt]] = col_ids[x]
            else:
                nb[row_ids[cnt]] = [col_ids[x]]
            cnt += 1

    with open(join(cutdir, 'nb.pkl'), 'wb') as f:
        cPickle.dump(nb, f, protocol=cPickle.HIGHEST_PROTOCOL)
