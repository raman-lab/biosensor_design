#!/usr/bin/env python
import argparse
import numpy as np
import sys


def boltzmann_probability_from_energies(energy_list, temperature):
    k = 0.0083144621    # boltzmann constant in kJ/molK
    energies = np.asarray(energy_list)
    frequencies = np.exp(-energies / (k * temperature))
    probabilities = frequencies / sum(frequencies)
    return probabilities


def conformer_probability(energies_file, temperature):
    energy_list = []
    with open(energies_file, 'r') as f:
        for line in f:
            energy_string = line.split('kJ/mol')[0].split()[-1]
            try:
                energy_list.append(float(energy_string))
            except ValueError:
                pass
    boltzmann_probs = boltzmann_probability_from_energies(energy_list, temperature)
    sys.stdout.write('energy\tprobability\n')
    for energy, probability in zip(energy_list, boltzmann_probs):
        sys.stdout.write('{0}\t{1}\n'.format(round(energy, 2), round(probability, 2)))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="""given input energies, output the probability of observing that 
    energy state from the boltzmanm distribution""")
    parser.add_argument('-t', '--temperature', type=float, default=298)
    required = parser.add_argument_group('required')
    required.add_argument('-e', '--energies', required=True,
                          help='file with one energy listed per line with units of kj/mol. script stores the number '
                               'appearing before "kJ/mol" in each line')
    args = parser.parse_args()
    conformer_probability(args.energies, args.temperature)
