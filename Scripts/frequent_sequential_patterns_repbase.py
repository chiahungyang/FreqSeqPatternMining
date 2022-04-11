"""Obtain frequent sequential patterns in the Arabdopsis thiliana data."""


from readfasta import Reader
from sigthrd import thresholds
from freqsubseq import frequent_patterns
import json

# Global variables
MIN_LEN = 4
THRD = 100
CONF = 0.9


def main():
    """
    Obtain sequential patterns with supports not less than the thresholds
    of significance.

    """
    # Compute the support thresholds of significance
    data = Reader('../Data/repbase_arabidopsis_thaliana_TE.fasta')
    num_nc = sum([len(seq) for seq in data.items()])
    len_thrd = thresholds(num_nc, MIN_LEN, THRD, CONF)

    # Search for frequent patterns of significance
    max_len = max(len_thrd.keys())
    results = frequent_patterns(data, MIN_LEN, max_len, THRD,
                                method='hybrid', verbose=True)
    sig_seq = [seq
               for seq, count in results
               if count >= len_thrd[len(seq)]]

    # Output the frequent patterns and their corresponding supports
    with open('../Results/frequent_sequential_patterns_repbase.json', 'w') as fp:
        json.dump(sig_seq, fp)

    return


if __name__ == '__main__':
    main()
