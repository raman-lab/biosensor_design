#!/usr/bin/env python
import argparse
import collections
import datetime as dt
import itertools
import numpy as np
import pandas as pd
import sys


class Cluster(object):
    def __init__(self, cluster_int):
        self._cluster_number = None
        self.cluster_number = cluster_int

    @property
    def cluster_number(self):
        return self._cluster_number

    @cluster_number.setter
    def cluster_number(self, cluster_int):
        self._cluster_number = int(cluster_int)


class Node(Cluster):
    def __init__(self, cluster_int, node_int):
        Cluster.__init__(self, cluster_int=cluster_int)
        self._node_number = None
        self._execute_time = None
        self._times_to_disconnect = []
        self._times_to_evict = []
        self._times_to_shadow = []
        self.node_number = node_int

    @property
    def node_number(self):
        return self._node_number

    @node_number.setter
    def node_number(self, node_int):
        self._node_number = int(node_int)

    @property
    def execute_time(self):
        return self._execute_time

    @execute_time.setter
    def execute_time(self, date_time_execution_began):
        _time = dt.datetime.strptime(date_time_execution_began, '%m/%d %H:%M:%S')
        self._execute_time = _time

    @property
    def times_to_evict(self):
        return self._times_to_evict

    @property
    def times_to_disconnect(self):
        return self._times_to_disconnect

    @property
    def times_to_shadow(self):
        return self._times_to_shadow

    def append_time(self, htcondor_indicator, date_time_str):
        _time = dt.datetime.strptime(date_time_str, '%m/%d %H:%M:%S')
        _time_difference = (_time - self._execute_time).total_seconds()
        if htcondor_indicator == "004":
            self._times_to_evict.append(_time_difference)
        elif htcondor_indicator == "007":
            self._times_to_shadow.append(_time_difference)
        elif htcondor_indicator == "022":
            self._times_to_disconnect.append(_time_difference)


def parse_htcondor_log_files(log_files, indicator_dict):
    node_dict = collections.defaultdict(dict)

    for log_file in log_files:
        with open(log_file, 'r') as f:
            for line in f:
                split_line = line.split()
                if split_line and split_line[0] in indicator_dict.keys():
                    indicator = split_line[0]
                    job_numbers, date, time = split_line[1:4]
                    cluster, node, step = map(int, job_numbers.strip('(').strip(')').split('.'))
                    date_time = '{0} {1}'.format(date, time)

                    if indicator == "000":
                        node_obj = Node(cluster_int=cluster, node_int=node)
                        node_dict[cluster][node] = node_obj
                        node_dict[cluster][node].execute_time = date_time
                    elif indicator == "001":
                        node_dict[cluster][node].execute_time = date_time
                    else:
                        node_dict[cluster][node].append_time(htcondor_indicator=indicator, date_time_str=date_time)
    return node_dict


def htcondor_log_analyzer(log_files):
    log_file_indicators = {
        "000": "submitted",
        "001": "executing",
        "004": "evicted",
        "007": "missing output file",
        "022": "disconnected"
    }

    cluster_node_dict = parse_htcondor_log_files(log_files, log_file_indicators)
    column_tups = list(itertools.product(
        ['evicted', 'disconnected', 'missing output file'],
        ['avg time (min)', 'std time (min)', 'n']
    ))
    indices = pd.MultiIndex.from_tuples(column_tups)
    for cluster, node_dict in cluster_node_dict.iteritems():
        node_matrix = pd.DataFrame(columns=indices)
        sys.stdout.write('cluster: {0}\n'.format(cluster))

        total_evict_array = np.array([])
        total_disconnect_array = np.array([])
        total_shadow_array = np.array([])

        for node, node_obj in node_dict.iteritems():
            evict_array = np.asarray(node_obj.times_to_evict)
            total_evict_array = np.append(total_evict_array, evict_array)

            disconnect_array = np.asarray(node_obj.times_to_disconnect)
            total_disconnect_array = np.append(total_disconnect_array, disconnect_array)

            shadow_array = np.asarray(node_obj.times_to_shadow)
            total_shadow_array = np.append(total_shadow_array, shadow_array)

            array_list = [evict_array, disconnect_array, shadow_array]

            array_to_append = []
            for array in array_list:
                array /= 60
                mean = np.mean(array)
                std = np.std(array)
                n = len(array)
                array_to_append.extend([mean, std, n])
            node_matrix.loc[node] = array_to_append
        sort_tup = zip(['evicted', 'disconnected', 'missing output file'], ['n'] * 3)
        node_matrix.sort_values(by=sort_tup, ascending=[False] * 3, inplace=True)


        total_dict = {
            'evicted': total_evict_array, 'disconnected': total_disconnect_array, 'missing output': total_shadow_array}

        with open('{0}.out'.format(cluster), 'w') as o:
            o.write('{0}\n'.format(node_matrix.to_string()))
            for key, total_array in total_dict.items():
                total_array /= 60
                total_mean = np.mean(total_array)
                total_std = np.std(total_array)
                total_n = len(total_array)
                o.write('cumulative {0}: {1} +/ {2} (n={3})\n'.format(key, total_mean, total_std, total_n))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="""script to calculate statistics about condor job processes. specifically, the number of
         occurrences and times of the following events relative to when the process began executing:
            004 - evicted
            007 - missing output file
            022 - disconnected
            """
    )
    required = parser.add_argument_group('required arguments')
    required.add_argument('-l', '--log', nargs='*', required=True, help='log file(s) of jobs')

    args = parser.parse_args()
    htcondor_log_analyzer(args.log)
