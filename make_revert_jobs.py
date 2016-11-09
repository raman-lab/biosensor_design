#!/usr/bin/env python

"""make shell scripts and submit files to run Rosetta revert_design_to_native and rescore protocols
on computing cluster. input pdb files are bunched in groups of 100 or less in compliance with
guideline to limit transfer files size to 40-50 MBs.
requires user has a TFs directory in /scratch, which has rescore xml and flags
if 'glob' is passed as an argument, all pdbs in the working dir will be used"""
import glob
import os
import sys
import getpass


def write_sh(name, wt_list):
    lines_pre = ["#!/bin/bash\n", "\n",
                 "tar -xzvf database.tar.gz\n", "\n", "# Revert design to native based on ddg\n",
                 "chmod +x revert_design_to_native.static.linuxgccrelease\n", "shopt -s nullglob\n"]
    loop_lines = []
    for wt in wt_list:
        if wt.startswith('./'):
            wt_base = wt.split('./')[1]
            wt_base = wt_base.split('.pdb')[0]
        else:
            wt_base = wt.split('.pdb')[0]
        loop_lines.append("for fname in {0}_*.pdb\n".format(wt_base))
        loop_lines.append("do\n")
        loop_lines.append("\t./revert_design_to_native.static.linuxgccrelease -database ./database "
                          "-extra_res_fa LG.params -score:weights enzdes -ex1 -ex2 "
                          "-revert_app:ddg_cycles 1 -revert_app:wt {0} -revert_app:design ${{fname}}"
                          " >> REVERT_${{fname%.pdb}}.log\n".format(wt))
        loop_lines.append("done\n")
    rescore_lines = ["\n", "# Rescore reverted designs\n", "chmod +x rosetta_scripts.static.linuxgccrelease\n",
                     "ls *.revert.pdb > infile\n",
                     "./rosetta_scripts.static.linuxgccrelease -database ./database @enzscore_flags "
                     "-parser::protocol enzscore_1.xml -in::file::fullatom -extra_res_fa LG.params "
                     "-out::file::scorefile rescore_{0}.sc -in:file:l infile > RESCORE_{0}.log\n".format(name)]
    lines_post = ["\n", "# Remove unneeded files\n",
                  "find . -type f -name '*.pdb' ! -name '*.revert_0001.pdb' -exec rm {} \;\n",
                  "rm -rf database/\n", "rm infile\n"]
    master_lines = []
    for l in [lines_pre, loop_lines, rescore_lines, lines_post]:
        master_lines.extend(l)
    with open("{0}.sh".format(name), 'w') as f:
        f.writelines(master_lines)


def write_sub(name, designs, wt_list):
    cwd = os.getcwd()
    user = getpass.getuser()
    tfs = ["LG.params", "/scratch/{0}/TFs/enzscore_flags".format(user), "/scratch/{0}/TFs/enzscore_1.xml".format(user),
           "http://proxy.chtc.wisc.edu/SQUID/nwhoppe/database.tar.gz",
           "http://proxy.chtc.wisc.edu/SQUID/nwhoppe/rosetta_scripts.static.linuxgccrelease",
           "http://proxy.chtc.wisc.edu/SQUID/nwhoppe/revert_design_to_native.static.linuxgccrelease"]
    tfs.extend(wt_list)
    tfs.extend(designs)
    tfs_comma = ', '.join(tfs)
    lines = ["universe = vanilla\n", "log = CONDOR_{0}.log\n".format(name), "error = CONDOR_{0}.err\n".format(name),
             "output = CONDOR_{0}.out\n".format(name), "\n", "executable = {0}.sh\n".format(name), "\n",
             'requirements = (OpSys == "LINUX")\n', "should_transfer_files = YES\n",
             "when_to_transfer_output = ON_EXIT\n", "\n",
             "# Leave the queue if job exits w/ 0 (normal exit code), else re-run job\n",
             "on_exit_remove = ExitCode =?= 0\n", "\n", "initialdir = {0}\n".format(cwd),
             "transfer_input_files = {0}\n".format(tfs_comma), "\n", "request_cpus = 1\n", "request_memory = 4GB\n",
             "request_disk = 4GB\n", "\n", "+WantFlocking = true\n", "\n", "queue 1\n"]
    with open("{0}.sub".format(name), 'w') as f:
        f.writelines(lines)


def get_wt(designs):
    wt_list = []
    for design in designs:
        design_split = design.split('_')
        wt = '_'.join(design_split[:5]) + '.pdb'
        if wt not in wt_list:
            wt_list.append(wt)
    return wt_list


def make_files(inputs):
    first = 0
    last = 49
    while last <= len(inputs):
        designs = inputs[first:last]
        name = "revert_{0}-{1}".format(first, last)
        wt_list = get_wt(designs)
        write_sh(name, wt_list)
        write_sub(name, designs, wt_list)
        first += 50
        last += 50
    designs = inputs[first:len(inputs)]
    name = "revert_{0}-{1}".format(first, len(inputs))
    wt_list = get_wt(designs)
    write_sh(name, wt_list)
    write_sub(name,  designs, wt_list)


if __name__ == '__main__':
    files = sys.argv[1:]
    if files == ['glob']:
        files = glob.glob('*.pdb')
    make_files(files)
