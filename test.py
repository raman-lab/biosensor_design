#!/usr/bin/env python
import gzip
from multiprocessing import Process, Manager
import time
import itertools


def do_work(in_queue, out_list):
    print 'entered'
    item = in_queue.get()
    ident, seq, space, qual = item

    # exit signal
    if item == None:
        return

    out_list[seq] = qual


if __name__ == "__main__":
    num_workers = 4

    manager = Manager()
    results = manager.dict()
    work = manager.Queue(num_workers)

    # start for workers
    pool = []
    for i in xrange(num_workers):
        p = Process(target=do_work, args=(work, results))
        p.start()
        pool.append(p)

    # produce data
    with gzip.open('/Users/nwhoppe/raman_lab/ngs/test2.fastq.gz') as f:
        # iters = itertools.chain(f, (None,)*num_workers)

        for fastq_chunk in itertools.chain(itertools.izip_longest(*[f] * 4, fillvalue=None), [None]*4):
            print fastq_chunk
            work.put(fastq_chunk)

    for p in pool:
        p.join()

    # get the results
    # example:  [(1, "foo"), (10, "bar"), (0, "start")]
    print results
