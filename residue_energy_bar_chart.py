#!/usr/bin/env python
import argparse
import numpy as np
import plotly
import plotly.graph_objs as go


def residue_energy_bar_chart(pdb_files, energies, chains, residues):
    # import glob
    # pdb_files = glob.glob('*_00*.pdb')
    energy_dict = {}

    for pdb_file in pdb_files:
        score_weights = {}
        index_term_dict = {}
        base_name = pdb_file.rsplit('_', 1)[0]
        fi = int(base_name.split('_')[-1])
        if base_name not in energy_dict.keys():
            energy_dict[base_name] = {}
            energy_dict[base_name]['binding_energy'] = []
            energy_dict[base_name]['fold_induction'] = fi
            for residue in residues:
                energy_dict[base_name][residue] = {}
                for chain in chains:
                    energy_dict[base_name][residue][chain] = {}
                    energy_dict[base_name][residue][chain]['residue_energy'] = []
                    for energy in energies:
                        energy_dict[base_name][residue][chain][energy] = []

        with open(pdb_file, 'r') as f:
            for line in f:
                if not line.startswith('REMARK'):
                    break
                if line.startswith('REMARK Binding energy:'):
                    binding_energy = float(line.rstrip().split()[-1])

                if line.startswith('REMARK Non-zero ScoreFunction weights:'):
                    line = f.next()
                    while line != 'REMARK\n':
                        split_line = line.rstrip().split()
                        score_term = split_line[1][0:10]
                        weight = float(split_line[2])
                        score_weights[score_term] = weight
                        line = f.next()
                if line.startswith('REMARK Weighted per-residue energies'):
                    score_terms_str = f.next().rstrip()
                    score_terms_str = score_terms_str.replace('REMARK  pdbtypech', '')
                    score_terms = [score_terms_str[i:i+10].replace(' ', '') for i in range(0, len(score_terms_str), 10)]
                    for t, term in enumerate(score_terms):
                        index_term_dict[t] = term
                    line = f.next()
                    while line != 'REMARK\n':
                        split_line = line.rstrip().split()
                        position = int(split_line[1])
                        chain = split_line[3]
                        if position in residues and chain in chains:
                            total_residue_score = 0.0
                            for s, score in enumerate(split_line[4:]):
                                score = float(score)
                                energy_term = index_term_dict[s]
                                if energy_term in energies:
                                    energy_dict[base_name][position][chain][energy_term].append(score)
                                    # if energy_term == 'fa_rep':
                                    #     total_residue_score += score_weights['fa_atr'] * score
                                    # else:
                                    total_residue_score += score_weights[energy_term] * score
                            energy_dict[base_name][position][chain]['residue_energy'].append(total_residue_score)
                        line = f.next()
            energy_dict[base_name]['binding_energy'].append(binding_energy)

    energies.append('residue_energy')
    for residue in residues:
        for energy in energies:
            data = []
            for pdb in energy_dict.keys():
                means = [np.mean(energy_dict[pdb][residue][c][energy]) for c in chains]
                stdevs = [np.std(energy_dict[pdb][residue][c][energy]) for c in chains]
                medians = [np.median(energy_dict[pdb][residue][c][energy]) for c in chains]
                mads = [np.median(np.abs(energy_dict[pdb][residue][c][energy] - medians[i])) for i, c in enumerate(
                     chains)]
                trace = go.Bar(
                    x=chains,
                    y=medians,
                    name=energy_dict[pdb]['fold_induction'],
                    error_y=dict(
                        type='data',
                        array=mads,
                        visible=True
                    )
                )
                data.append(trace)
            layout = go.Layout(
                barmode='group'
            )
            fig = go.Figure(data=data, layout=layout)
            plotly.offline.plot(fig, filename='{0}_{1}.html'.format(energy, residue), auto_open=False)

    x = []
    y = []
    y_err = []
    for pdb in energy_dict.keys():
        fi = energy_dict[pdb]['fold_induction']
        bind_e_list = energy_dict[pdb]['binding_energy']
        x.append(pdb)
        # y.append(np.mean(bind_e_list))
        # y_err.append(np.std(bind_e_list))
        y.append(np.median(bind_e_list))
        y_err.append(np.median(np.abs(bind_e_list - np.median(bind_e_list))))
    data = [go.Bar(
        x=x,
        y=y,
        error_y=dict(
            type='data',
            array=y_err,
            visible=True
        )
    )]
    plotly.offline.plot(data, filename='binding_energy.html', auto_open=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="""script to plot bar charts of residue specific energies from rosetta.
                                                 input are pdbs with residue specific energies written at the top of
                                                 the file. multiple nstructs will be lumped together""")
    parser.add_argument('-r', '--residues', nargs='*', type=int,
                        default=[5, 9, 12, 34, 35, 36, 40, 45, 47, 48, 50, 51],
                        help='residues to use in making bar charts - one chart for each residue')
    parser.add_argument('-c', '--chains', nargs='*', default=['A', 'B'],
                        help='chains in pdb. multiple chains are on the same bar chart')
    parser.add_argument('-e', '--energies', nargs='*',
                        default=['fa_atr', 'fa_rep', 'hbond_sc', 'fa_dun'],
                        help='energies to make bar charts for. Binding energy is not residue specific and total will '
                             'use the ScoreFunction weights to compute linear combination of terms')
    required = parser.add_argument_group('required')
    required.add_argument('-p', '--pdbs', nargs='*', required=True,
                          help='pdbs with residue specific energies written as REMARKS')
    args = parser.parse_args()
    residue_energy_bar_chart(args.pdbs, args.energies, args.chains, args.residues)
