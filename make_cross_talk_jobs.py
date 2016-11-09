#!/usr/bin/env python

"""script to create pdb files, submit files, and shell executables to run Rosetta cross talk jobs through HTCondor"""
import argparse
import getpass
import os
import tarfile


def get_lig_name(pose_name):
    name = '_'.join(pose_name.split('_')[1:])
    return name


def make_pdb(design_pose, ligand_pose, ligand_name):
    design_base_name = design_pose.rsplit('.pdb', 1)[0]
    cross_talk_pose = '{0}_w_{1}'.format(design_base_name, ligand_name)
    lines_to_add = []
    with open(ligand_pose, 'r') as lig_f:
        for line in lig_f:
            if line.startswith('HETATM'):
                lines_to_add.append(line)
    with open(design_pose, 'r') as design_f:
        with open(cross_talk_pose, 'w') as new_f:
            lines_to_write = []
            for line in design_f:
                if line.startswith('HETATM'):
                    if lines_to_add:
                        lines_to_write.extend(lines_to_add)
                        lines_to_add = None
                    continue
                lines_to_write.append(line)
            new_f.writelines(lines_to_write)
    return cross_talk_pose


def make_tar_ball(transfer_files, base_name):
    tar_name = '{0}.tar.gz'.format(base_name)
    tar = tarfile.open(tar_name, "w:gz")
    for tf in transfer_files:
        tar.add(tf)
    tar.close()
    return tar_name


def make_executable(tar_name, protocol, nstruct, base_name, flags, resfile, params_file):
    executable_name = '{0}.sh'.format(base_name)
    executable_lines = [
        "#!/bin/bash\n",
        "\n",
        "tar -xzf database.tar.gz\n",
        "tar -xzvf {0} > cross_talk_infile\n".format(tar_name),
        "chmod +x rosetta_scripts.static.linuxgccrelease\n",
        "\n",
        "# Run Rosetta cross talk protocol\n",
        "./rosetta_scripts.static.linuxgccrelease -database ./database -parser:protocol {0} "
        "-extra_res_fa {1} @{2} -parser:script_vars ligchain=X resfile={3} -out:file:o design_pdbs "
        "-in:file:l cross_talk_infile -out:file:scorefile score_{4}.sc -nstruct {5} -out:suffix _ct "
        "> RUN_{4}.log\n".format(protocol, params_file, flags, resfile, base_name, nstruct),
        "\n",
        "# Static rescore\n"
        "ls *ct_00??.pdb > rescore_cross_talk_infile\n",
        "./rosetta_scripts.static.linuxgccrelease -database ./database @enzscore_flags "
        "-parser::protocol enzscore_1.xml -in:file:fullatom -extra_res_fa {0} "
        "-out:file:scorefile rescore_{1}.sc -in:file:l rescore_cross_talk_infile -out:suffix _rs_${{1}} "
        "> RESCORE_{1}.log\n".format(params_file, base_name),
        "\n",
        "# Remove unneeded files\n",
        "find . -maxdepth 1 -type f -name '*.pdb' ! -name '*ct_00??_rs_*_0001.pdb' -exec rm {} \;\n"
        "rm *infile\n"
        "rm -rf database/\n"
    ]
    with open(executable_name, 'w') as executable:
        executable.writelines(executable_lines)
    return executable_name


def make_submit(tar_name, executable_name, protocol, flags, resfile, params_file):
    base_file = executable_name.split('.sh')[0]
    submit_name = '{0}.sub'.format(base_file)
    cwd = os.getcwd()
    user = getpass.getuser()
    home = os.path.expanduser('~{0}'.format(user))
    tfs = ['{0}/{1}'.format(cwd, tar_name),
           '{0}/TFs/{1}, {0}/TFs/{2}, {0}/TFs/{3}, {0}/TFs/{4}'.format(home, protocol, flags, resfile, params_file),
           '{0}/TFs/enzscore_1.xml, {0}/TFs/enzscore_flags'.format(home),
           "http://proxy.chtc.wisc.edu/SQUID/nwhoppe/database.tar.gz",
           "http://proxy.chtc.wisc.edu/SQUID/nwhoppe/rosetta_scripts.static.linuxgccrelease"]
    tfs_comma = ', '.join(tfs)
    lines = [
        'universe = vanilla\n',
        'log = CONDOR_{0}_$(Cluster).log\n'.format(base_file),
        'error = CONDOR_{0}_$(Cluster).err\n'.format(base_file),
        'output = CONDOR_{0}_$(Cluster).out\n'.format(base_file),
        '\n',
        'executable = {0}\n'.format(executable_name),
        'arguments = $(Cluster)\n',
        '\n',
        'requirements = (OpSys == "LINUX")\n',
        'should_transfer_files = YES\n',
        'when_to_transfer_output = ON_EXIT\n',
        'on_exit_remove = ExitCode =?= 0\n',
        '\n',
        'initial_dir = {0}/output\n'.format(cwd),
        'transfer_input_files = {0}\n'.format(tfs_comma),
        '\n',
        'request_cpus = 1\n',
        'request_memory = 2GB\n',
        'request_disk = 2GB\n',
        '\n',
        '+WantFlocking = true\n',
        '+WantGlidein = true \n',
        '\n',
        'queue 1\n'
    ]
    with open(submit_name, 'w') as sub:
        sub.writelines(lines)


def get_file_base_name(first_file, last_file):
    name_list = []
    for file_name in [first_file, last_file]:
        if '_w_' in file_name:
            part_list = file_name.split('.pdb')
            prefix = part_list[0]
            w_index = part_list[1].find('w')
            suffix = part_list[1][w_index:]
            name = '_'.join([prefix, suffix])
            name_list.append(name)
        else:
            name = file_name.split('.pdb', 1)[0]
            name_list.append(name)
    base_name = '{0}_to_{1}'.format(name_list[0], name_list[1])
    return base_name


def make_condor_files(poses, protocol, nstruct, poses_per_executable, flags, resfile, params_file):
    poses.sort()
    first = 0
    last = poses_per_executable
    while last <= len(poses):
        transfer_files = poses[first:last]
        base_name = get_file_base_name(transfer_files[0], transfer_files[-1], )
        tar_name = make_tar_ball(transfer_files, base_name)
        executable_name = make_executable(tar_name, protocol, nstruct, base_name, flags, resfile, params_file)
        make_submit(tar_name, executable_name, protocol, flags, resfile, params_file)
        first += poses_per_executable
        last += poses_per_executable

    transfer_files = poses[first:len(poses) + 1]
    base_name = get_file_base_name(transfer_files[0], transfer_files[-1], )
    tar_name = make_tar_ball(transfer_files, base_name)
    executable_name = make_executable(tar_name, protocol, nstruct, base_name, flags, resfile, params_file)
    make_submit(tar_name, executable_name, protocol, flags, resfile, params_file)


def make_cross_talk_files(designs, ligands, protocol, flags, params_file, nstruct, poses_per_executable, resfile):
    cross_talk_poses = []
    if ligands:
        for design_pose in designs:
            for ligand_pose in ligands:
                ligand_name = get_lig_name(ligand_pose)
                cross_talk_pose = make_pdb(design_pose, ligand_pose, ligand_name)
                cross_talk_poses.append(cross_talk_pose)
        make_condor_files(cross_talk_poses, protocol, nstruct, poses_per_executable, flags, resfile, params_file)
    else:
        make_condor_files(designs, protocol, nstruct, poses_per_executable, flags, resfile, params_file)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="script to create pdb files, submit files, shell executables, and "
                                                 "tar balls to run Rosetta cross talk jobs through HTCondor. "
                                                 "assumes two directories exist: ./output ~/TFs/ "
                                                 "protocol, resfile, and params file are assumed to be in ~/TFs/ "
                                                 "run ./make_cross_talk_jobs.py --help for option descriptions "
                                                 "script creates files and tars combinatorially, so it may take awhile "
                                                 "to run. Omit -l to run cross talk protocol with initial designed "
                                                 "ligand as a control")
    parser.add_argument("-n", "--nstruct", type=int, default=1,
                        help="Rosetta option:Number of times to process each input PDB")
    parser.add_argument("-ppe", "--poses_per_executable", type=int, default=50)
    parser.add_argument("-r", "--resfile", default="empty.resfile",
                        help="Name of resfile used by protocol")
    parser.add_argument("-l", "--ligands", nargs='*',
                        help="Names of single PDB files with alternate ligand to process or file(s) "
                             "containing list(s) of PDB files with alternate ligand to process. "
                             "if no ligand name is detected in PDB file name, 'wt' ligand will be assumed")
    requiredO = parser.add_argument_group('required arguments')
    requiredO.add_argument("-d", "--designs", nargs='*', required=True,
                           help="Names of single design PDB files to process or file(s) containing list(s) of "
                                "design PDB files to process")
    requiredO.add_argument("-x", "--xml_protocol", help="Name of xml file used in Rosetta command line")
    requiredO.add_argument("-f", "--flags", required=True,
                           help="Name of flags file that accompanies protocol")
    requiredO.add_argument("-p", "--params", required=True,
                           help='Name of params file for ligand (can only handle one params file)')
    args = parser.parse_args()
    make_cross_talk_files(args.designs, args.ligands, args.xml_protocol, args.flags, args.params,
                          args.nstruct, args.poses_per_executable, args.resfile)
