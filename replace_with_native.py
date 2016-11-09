"""script to replace current ligand in pdb with its native ligand
native ligand is input with its native coordinates
need input descriptor file
script meant to be invoked with pymol -
    pymol -qc replace_current_ligand.py -- pair_fit_descriptor.txt"""

import re
import sys
import glob


def parse_input_file(in_file):
    with open(in_file, 'r') as f:
        all_lines = f.readlines()
        content = [line.rstrip('\n') for line in all_lines if '#' not in line]
        pdb_to_glob = content[0].split(': ')[1]
        pdb_list = glob.glob(pdb_to_glob)
        native_to_glob = content[1].split(': ')[1]
        native_lig_list = glob.glob(native_to_glob)
    return pdb_list, native_lig_list


def reformat_lig_pdb(ligands):
    """reformat conformer file so that output is Rosetta compatible"""
    for lig_pdb in ligands:
        with open(lig_pdb, 'r') as f:
            write_lines = []
            for line in f:
                if 'HETATM' in line and 'LG1 X' not in line:
                    line_list = re.split(r'(\s+)', line)
                    line_list[6] = 'LG1 X'
                    line_list[8] = '1'
                    write_line = ''.join(line_list)
                    write_lines.append(write_line)
                else:
                    write_lines.append(line)
        with open(lig_pdb, 'w') as f:
            f.writelines(write_lines)


def create_pdb_w_lig(pdb, ligands):
    """create new pdb with input ligand"""
    pdb_base = pdb.rsplit('.pdb', 1)[0]
    base = '//X/LG1`1/'
    lig_identifier = base[4:base.find('`')]
    for lig in ligands:
        lig_base = lig.split('.')[0]
        cmd.reinitialize()
        cmd.load(pdb)
        cmd.select('og_lig', 'resname ' + lig_identifier)
        cmd.load(lig)
        cmd.remove('og_lig')
        cmd.set('pdb_use_ter_records', 0)
        cmd.save('{0}_w_{1}.pdb'.format(pdb_base, lig_base))


def replace_ligand(in_file):
    pdb_list, native_lig_list = parse_input_file(in_file)

    reformat_lig_pdb(native_lig_list)
    for pdb in pdb_list:
        create_pdb_w_lig(pdb, native_lig_list)

input_file = sys.argv[-1]
replace_ligand(input_file)
