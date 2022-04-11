"""Study the computational costs of frequent patterns finding methods."""


import numpy as np
from random import random
from readfasta import Reader
from sigthrd import thresholds
from freqsubseq import frequent_patterns, Counter, ExceedAllocatedMemory
import json
import logging

logging.basicConfig(format='%(asctime)s - %(message)s', level=logging.INFO)

# Global variables
NUM_SEQ = 118554  # Total number of sequences pre-read from the data
FILEPATH = '../Data/EpiR_all_seq.fasta'
THRD = 100
MIN_LEN = 3
CONF = 0.9


def experiment(prob):
    """
    Sample sequences with a given probability from the data, and obtain the
    number of passes for A-priori and position-based method.

    """
    # Sample sequences and compute its total number of nucleotides
    seq_IDs = [idx for idx in range(NUM_SEQ) if random() < prob]
    data = Reader(FILEPATH, seq_IDs)
    num_nc = sum([len(seq) for seq in data.items()])
    len_thrd = thresholds(num_nc, MIN_LEN, THRD, CONF)

    # Re-sample if the minimum pattern length does not have a threshold of
    # significance
    while MIN_LEN not in len_thrd.keys():
        seq_IDs = [idx for idx in range(NUM_SEQ) if random() < prob]
        data = Reader(FILEPATH, seq_IDs)
        num_nc = sum([len(seq) for seq in data.items()])
        len_thrd = thresholds(num_nc, MIN_LEN, THRD, CONF)

    # Obtain the maximum pattern length of significance that satisfies the
    # pre-specified support lower bound
    max_len = max(len_thrd.keys())

    # logging.info('Sampling is completed.')
    results = list()

    # A-priori
    logging.info('A-priori starts.')
    counter = Counter()
    _ = frequent_patterns(data, MIN_LEN, max_len, THRD, counter=counter,
                          method='apriori', verbose=True)
    results.append(('apriori', num_nc, counter.count))

    # Position-based
    logging.info('Position-based method starts.')
    try:
        counter = Counter()
        _ = frequent_patterns(data, MIN_LEN, max_len, THRD, counter=counter,
                              method='position', verbose=True)
        results.append(('position', num_nc, counter.count))
    except ExceedAllocatedMemory:
        logging.info('The queue size exceeds the amount of allocated memory.')

    return results


def main():
    """Run experiments on sampled sequences."""
    probs = np.logspace(np.log10(5e-5), np.log10(5e-3), num=10)
    output = list()
    for prob in probs:
        output += experiment(prob)
        logging.info('An experiment is done.\n')
    with open('../Results/computational_costs.json', 'w') as fp:
        json.dump(output, fp)

    return


if __name__ == '__main__':
    main()
