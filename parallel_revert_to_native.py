#!/usr/bin/env python

""""Run Rosetta revert_design_to_native application on multiple processors"""

import argparse
import multiprocessing as mp
import subprocess


def rosetta_revert(design):
    """Call rosetta revert_design_to_native and wait until completion"""
    database = '/'.join(rosetta_path.split('/')[:-3]) + '/database/'
    log_name = design.rsplit('.', 1)[0] + '.log'
    with open(log_name, 'w') as output:
        command = "{0} -database {1} -extra_res_fa {2} -score:weights enzdes -ex1 -ex2 " \
                  "-revert_app:post_repack {3} -revert_app:threshold {4} -revert_app:ddg_cycles {5} " \
                  "-revert_app:wt {6} -revert_app:design {7}".format(
            rosetta_path, database, params_file, repack, thresh, ddgs, wt_pdb, design)
        split_command = command.split()
        subprocess.call(split_command, stdout=output)


def revert(designs, wt, params, post_repack, threshold, ddg_cycles, cores, rosetta):
    """Create pool of workers and use them to call rosetta"""
    global wt_pdb
    global params_file
    global repack
    global thresh
    global ddgs
    global rosetta_path

    wt_pdb = wt
    params_file = params
    repack = post_repack
    thresh = threshold
    ddgs = ddg_cycles
    rosetta_path = rosetta

    if cores:
        pool = mp.Pool(cores)
    else:
        cores = mp.cpu_count() - 1
        pool = mp.Pool(cores)
    pool.map_async(rosetta_revert, designs)
    pool.close()
    pool.join()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Run Rosetta revert_design_to_native application on multiple processors")
    parser.add_argument("-r", "--rosetta",
                        default='~/rosetta_src_2015.38.58158_bundle/main/source/bin/revert_design_to_native.static.linuxgccrelease',
                        help="""path to rosetta script (default:
                        ~/rosetta_src_2015.38.58158_bundle/main/source/bin/
                        revert_design_to_native.static.linuxgccrelease""")
    parser.add_argument("-c", "--cores", type=int,
                        help="number of processors to use (default is max number available - 1)")
    parser.add_argument("-d", "--ddg_cycles", type=int, default=1,
                        help="ddg iterations (default: 1")
    parser.add_argument("-t", "--threshold", type=float, default=0.5,
                        help="ddg threshold for acceptance of reversion (default: 0.5 )")
    parser.add_argument("-k", "--post_repack", default="false",
                        help="attempt repacking of the reverted structure prior to output. "
                             "true or false. (default: false)")
    parser.add_argument("-p", "--params", default='LG.params',
                        help="name of LG.params file (default: LG.params")
    requiredO = parser.add_argument_group('required arguments')
    requiredO.add_argument("-w", "--wt", required=True,
                           help="pdb for wild type protein")
    requiredO.add_argument("designs", nargs='*',
                           help="one or more designed pdb. reversion is split between the processors")

    args = parser.parse_args()
    revert(args.designs, args.wt, args.params, args.post_repack,
           args.threshold, args.ddg_cycles, args.cores, args.rosetta)
