"""script to replace current ligand in pdb with input ligand
need input descriptor file
script invoked with pymol -
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
        lig_to_glob = content[1].split(': ')[1]
        lig_list = glob.glob(lig_to_glob)
        current_atoms = content[2].split(': ')[1]
        current_atom_list = current_atoms.split()
        input_atoms = content[3].split(': ')[1]
        input_atom_list = input_atoms.split()
    return pdb_list, lig_list, current_atom_list, input_atom_list


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


def create_pdb_w_lig(pdb, ligands, current_lig_atoms, input_lig_atoms):
    """create new pdb with input ligand"""
    pdb_base = pdb.rsplit('.pdb', 1)[0]
    base = '//X/LG1`1/'
    lig_identifier = base[4:base.find('`')]
    for lig in ligands:
        lig_base = lig.split('.')[0]
        cmd.reinitialize()
        cmd.load(pdb)
        cmd.select('ogLig', 'resname ' + lig_identifier)
        cmd.load(lig)

        atom_count = len(current_lig_atoms)
        if atom_count == 3:
            cmd.pair_fit(pdb_base + base + current_lig_atoms[0], lig_base + base + input_lig_atoms[0],
                         pdb_base + base + current_lig_atoms[1], lig_base + base + input_lig_atoms[1],
                         pdb_base + base + current_lig_atoms[2], lig_base + base + input_lig_atoms[2])
        elif atom_count == 4:
            cmd.pair_fit(pdb_base + base + current_lig_atoms[0], lig_base + base + input_lig_atoms[0],
                         pdb_base + base + current_lig_atoms[1], lig_base + base + input_lig_atoms[1],
                         pdb_base + base + current_lig_atoms[2], lig_base + base + input_lig_atoms[2],
                         pdb_base + base + current_lig_atoms[3], lig_base + base + input_lig_atoms[3])
        elif atom_count == 5:
            cmd.pair_fit(pdb_base + base + current_lig_atoms[0], lig_base + base + input_lig_atoms[0],
                         pdb_base + base + current_lig_atoms[1], lig_base + base + input_lig_atoms[1],
                         pdb_base + base + current_lig_atoms[2], lig_base + base + input_lig_atoms[2],
                         pdb_base + base + current_lig_atoms[3], lig_base + base + input_lig_atoms[3],
                         pdb_base + base + current_lig_atoms[4], lig_base + base + input_lig_atoms[4])
        elif atom_count == 6:
            cmd.pair_fit(pdb_base + base + current_lig_atoms[0], lig_base + base + input_lig_atoms[0],
                         pdb_base + base + current_lig_atoms[1], lig_base + base + input_lig_atoms[1],
                         pdb_base + base + current_lig_atoms[2], lig_base + base + input_lig_atoms[2],
                         pdb_base + base + current_lig_atoms[3], lig_base + base + input_lig_atoms[3],
                         pdb_base + base + current_lig_atoms[4], lig_base + base + input_lig_atoms[4],
                         pdb_base + base + current_lig_atoms[5], lig_base + base + input_lig_atoms[5])

        cmd.remove('ogLig')
        cmd.set('pdb_use_ter_records', 0)
        cmd.save('{0}_w_{1}.pdb'.format(pdb_base, lig_base))


def replace_ligand(in_file):
    pdb_list, lig_list, current_lig_atoms, input_lig_atoms = parse_input_file(in_file)
    assert len(current_lig_atoms) == len(input_lig_atoms), 'Error: must have same number of current and input' \
                                                           'ligand atoms'
    reformat_lig_pdb(lig_list)
    for pdb in pdb_list:
        create_pdb_w_lig(pdb, lig_list, current_lig_atoms, input_lig_atoms)

input_file = sys.argv[-1]
replace_ligand(input_file)
