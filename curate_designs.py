#!/usr/bin/env python
import argparse
import re
import collections
from itertools import izip_longest


def pair_files(unique_sc, filtered_sc, fasta):
    """check if input files meet requirements and extract ligand names"""
    assert len(unique_sc) == len(filtered_sc) == len(fasta), "Number of input fasta files, unique score files, " \
                                                             "and filtered score files must be the same"
    lig_file_d = collections.defaultdict(dict)
    lig_list = []
    for entry in filtered_sc:
        fb_name = entry.split('.')[0].lower()
        f_parts = fb_name.split('_')
        # best way I currently know to dynamically pull out lig names
        # search terms have to be updated as naming conventions changes
        search = ['unique', 'score', 'all', 'rescore', 'revert', 'unique', 'filtered', 'filter', 'cut']
        for part in search:
            try:
                f_parts.remove(part)
            except ValueError:
                pass
        lig = [x for x in f_parts if x.isalpha()]
        lig_list.extend(lig)

    for lig in lig_list:
        file_d = {'unique_sc': [], 'filtered_sc': [], 'fasta': []}
        for filtered in filtered_sc:
            filtered = filtered.lower()
            if lig in filtered:
                file_d['filtered_sc'].append(filtered)
        for unique in unique_sc:
            unique = unique.lower()
            if lig in unique:
                file_d['unique_sc'].append(unique)
        for fas in fasta:
            fas = fas.lower()
            if lig in fas:
                file_d['fasta'].append(fas)

        assert len(file_d['unique_sc']) == len(file_d['filtered_sc']) == len(file_d['fasta']) == 1, \
            "Can only handle one fasta file, one unique score file, and one filtered score file per ligand"
        lig_file_d[lig] = file_d
    return lig_file_d


def parse_combos(sc_file, unique):
    """make dictionary with A<int>_P<int> as key and values are list of descriptions from score files"""
    combo_desc = collections.defaultdict(list)
    with open(sc_file, 'r') as f:
        header = f.readline().split()
        lig_ts_index = header.index('SR_1_total_score')
        for line in f:
            line_list = line.split()
            if (not line_list) or (line_list[0].startswith("#")) or (line_list[0][0].isalpha()):
                continue
            try:
                desc = line_list[-1].split('.')[0]  # .lower() might add if sc files and fasta files don't agree
                ap_combo = re.search('A([0-9]+)_P([0-9]+)', desc).group()
            except AttributeError:
                continue
            if unique:
                lig_ts = float(line_list[lig_ts_index])
                combo_desc[ap_combo].append((lig_ts, desc))
            else:
                combo_desc[ap_combo].append(desc)
    return combo_desc


def make_limits(nested_filtered_d, lig_list, max_seq, min_per_combo):
    """calculate max sequences per ligand and minimum sequences per 'combo' from total input max
     store in hash table"""
    nested_limits = collections.defaultdict(dict)
    total_seqs = 0.0
    for lig in lig_list:
        lig_info = collections.defaultdict(int)
        lig_seqs = 0.0
        for combo in (nested_filtered_d[lig].keys()):
            combo_seqs = len(nested_filtered_d[lig][combo])
            total_seqs += combo_seqs
            lig_seqs += combo_seqs
        max_per_lig = int((lig_seqs / total_seqs) * max_seq)
        lig_info['max_per_lig'] = max_per_lig
        lig_info['min_per_combo'] = min_per_combo
        nested_limits[lig] = lig_info
    return nested_limits


def enforce_limits(dict_to_curate, nested_tuples_to_add, nested_limits, lig_list, blind):
    for lig in lig_list:
        for combo in dict_to_curate[lig].keys():
            remaining_int = len(nested_tuples_to_add[lig][combo]) - len(dict_to_curate[lig][combo])
            if nested_limits[lig]['min_per_combo'] < remaining_int:
                minimum = nested_limits[lig]['min_per_combo']
            else:
                minimum = remaining_int
            sorted_tup = sorted(nested_tuples_to_add[lig][combo])
            sorted_score, sorted_desc = zip(*sorted_tup)
            count = 0
            desc_to_add = []
            if blind:
                while len(desc_to_add) < minimum:
                    if sorted_desc[count] not in dict_to_curate[lig][combo]:
                        desc_to_add.append(sorted_desc[count])
                    count += 1
                dict_to_curate[lig][combo].extend(desc_to_add)
            else:
                combo_len = len(dict_to_curate[lig][combo])
                if combo_len < minimum:
                    diff = minimum - combo_len
                    while len(desc_to_add) < diff:
                        if sorted_desc[count] not in dict_to_curate[lig][combo]:
                            desc_to_add.append(sorted_desc[count])
                        count += 1
                    dict_to_curate[lig][combo].extend(desc_to_add)
    return dict_to_curate


def write_fasta(curated_d, file_dict):
    for lig in file_dict.keys():
        curated_fasta = []
        fasta = file_dict[lig]['fasta'][0]
        with open(fasta, 'r') as f:
            for pdb, seq in izip_longest(f, f, fillvalue=None):
                if '/' in pdb:
                    desc = (pdb.split('/')[-1]).split('.')[0]
                elif pdb.startswith('>'):
                    desc = (pdb.split('>')[1]).split('.')[0]
                else:
                    continue
                combo = re.search('A([0-9]+)_P([0-9]+)', desc).group()
                if desc in curated_d[lig][combo]:
                    curated_fasta.append(pdb)
                    curated_fasta.append(seq)
        with open('curated_{0}.fasta'.format(lig), 'w') as o:
            o.writelines(curated_fasta)


def curate(unique_sc, filtered_sc, fasta, max_seq, min_per_combo, blind):
    """main function"""
    lig_file_d = pair_files(unique_sc, filtered_sc, fasta)
    nested_filtered_d = collections.defaultdict(dict)
    nested_unique_td = collections.defaultdict(dict)
    for lig in lig_file_d.keys():
        sc_f = lig_file_d[lig]['filtered_sc'][0]
        combo_d = parse_combos(sc_f, unique=False)
        nested_filtered_d[lig] = combo_d
        sc_u = lig_file_d[lig]['unique_sc'][0]
        tuple_d = parse_combos(sc_u, unique=True)
        nested_unique_td[lig] = tuple_d
    nested_limits = make_limits(nested_filtered_d, lig_file_d.keys(), max_seq, min_per_combo)
    curated_d = enforce_limits(nested_filtered_d, nested_unique_td, nested_limits, lig_file_d.keys(), blind)
    write_fasta(curated_d, lig_file_d)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="""Script to generate curated fastas for a given scaffold
        from filtered score files for each target ligand.
        the number of seqs per target ligand is scaled by the proportion to all score file entries
        inputs:
            unique and filtered score files following the revert_design_to_native Rosetta application
            max number of sequences for scaffold
        assumes naming convention: <scaffold>_<lig>_A<int>_P<int>_...
        run script with -h for list of options and default settings"""
    )
    parser.add_argument("-min", "--min_per_combo", type=int, default=100,
                        help="min number of seqs for a ligand per A<int>_P<int> combo")
    parser.add_argument("-max", "--max_seq", type=int, default=50000,
                        help="total size of curated set for the scaffold (not by ligand or combo)"
                             "before minimum is added")
    parser.add_argument("-b", "--blind", action="store_true",
                        help="blindly add minimum to each A<int>_P<int> combo and do not consider how many designs"
                             "are already present in the combo")
    requiredO = parser.add_argument_group('required arguments')
    requiredO.add_argument("-a", "--fasta", nargs='*', required=True,
                           help="one or more unique and reverted fasta files. "
                                "seq in fasta files must contain the ligand name")
    requiredO.add_argument("-f", "--filtered_sc", nargs='*', required=True,
                           help="one or more filtered score files with partner unique score file. "
                                "score files are partnered by ligand name, which must be consistent "
                                "with fasta files")
    requiredO.add_argument("-u", "--unique_sc", nargs='*', required=True,
                           help="one or more unique score files with partner filtered score file. "
                                "score files are partnered by ligand name, which must be consistent "
                                "wit fasta files")
    args = parser.parse_args()
    curate(args.unique_sc, args.filtered_sc, args.fasta, args.max_seq, args.min_per_combo, args.blind)
