"""Compare the results from genome data to the RepBase TE database."""


import json
import matplotlib.pyplot as plt

plt.style.use('ggplot')


def fraction_rel_comp(A, B):
    """
    Return the fraction of elements in set A that belong to the relative
    complement of set B.

    """
    return len(A - B) / len(A)


def jaccard_similarity(A, B):
    """Return the Jaccard similarity between sets A and B."""
    return len(A & B) / len(A | B)


def main():
    """
    Compare the frequent sequential patterns found in the Arabdopsis thiliana
    genome data to ones found in the RepBase TE database by plotting, for each
    pattern length:
    1. Number of frequent patterns found
    2. Fraction of patterns that belong to the relative complement of the other
    3. Jaccard similarity between the two outcomes

    """
    # Load outcomes of frequent sequential pattern findings
    filepath = '../Results/'
    with open(filepath + 'frequent_sequential_patterns.json', 'r') as fp:
        freq_seq_data = json.load(fp)
    with open(filepath + 'frequent_sequential_patterns_repbase.json', 'r') as fp:
        freq_seq_repb = json.load(fp)

    # Structure the outcomes as sets of sequences according to their length
    len_seqset_data = dict()
    for seq in freq_seq_data:
        m = len(seq)
        len_seqset_data.setdefault(m, set())
        len_seqset_data[m].add(seq)

    len_seqset_repb = dict()
    for seq in freq_seq_repb:
        m = len(seq)
        len_seqset_repb.setdefault(m, set())
        len_seqset_repb[m].add(seq)

    # Plot the number of frequent patterns found for each length
    fig, ax = plt.subplots()
    ax.set_xlabel('Pattern Length', fontsize=16)
    ax.set_ylabel('Number of Patterns Found', fontsize=16)
    _len, num = zip(*[(l - 0.2, len(s)) for l, s in len_seqset_repb.items()])
    ax.bar(_len, num, width=0.3, color='C1', label='RepBase')
    _len, num = zip(*[(l + 0.2, len(s)) for l, s in len_seqset_data.items()])
    ax.bar(_len, num, width=0.3, color='C0', label='Genome Data')
    ax.set_yscale('log')
    ax.legend()
    fig.savefig('../Results/number_of_patterns_found.png', dpi=300)

    # Extract range of common pattern lengths in the two outcomes
    min_len = max(min(len_seqset_data.keys()), min(len_seqset_repb.keys()))
    max_len = min(max(len_seqset_data.keys()), max(len_seqset_repb.keys()))
    lengths = range(min_len, max_len + 1)

    # Obtain the fraction of relative complement for each pattern length
    frac_comp_data = [fraction_rel_comp(len_seqset_data[_len],
                                        len_seqset_repb[_len])
                      for _len in lengths]
    frac_comp_repb = [fraction_rel_comp(len_seqset_repb[_len],
                                        len_seqset_data[_len])
                      for _len in lengths]

    # Plot the fractionn of relative complement for each length
    fig, ax = plt.subplots()
    ax.set_xlabel('Pattern Length', fontsize=16)
    ax.set_ylabel('Fraction of Relative Complement', fontsize=16)
    ax.plot(lengths, frac_comp_repb, marker='o', color='C1', label='RepBase')
    ax.plot(lengths, frac_comp_data, marker='o', color='C0', label='Genome Data')
    ax.set_xticks(lengths)
    ax.set_xticklabels(lengths)
    ax.legend()
    fig.savefig('../Results/fraction_of_relative_complement.png', dpi=300)

    # Plot the Jaccard similarity between the two outcomes
    jaccard = [jaccard_similarity(len_seqset_data[_len], len_seqset_repb[_len])
               for _len in lengths]
    fig, ax = plt.subplots()
    ax.set_xlabel('Pattern Length', fontsize=16)
    ax.set_ylabel('Jaccard Similarity', fontsize=16)
    ax.plot(lengths, jaccard, marker='o', color='C2')
    ax.set_xticks(lengths)
    ax.set_xticklabels(lengths)
    fig.savefig('../Results/jaccard_similarity.png', dpi=300)

    return


if __name__ == '__main__':
    main()
