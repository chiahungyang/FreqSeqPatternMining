"""Summarize the Arabdopsis thiliana nucleotide sequence data."""


from readfasta import sequences
import numpy as np
import matplotlib.pyplot as plt
import logging

plt.style.use('ggplot')
logging.basicConfig(format='%(asctime)s - %(message)s', level=logging.INFO)


def main():
    """
    Obtain the sequence length distribution of Arabdopsis thiliana data, total
    number of sequences and number of base pairs.

    """
    # Read data in to sequence length
    logging.info('Start.')
    seq_len = [len(seq) for seq in sequences('../Data/EpiR_all_seq.fasta')]
    logging.info('Data loaded.')

    # Output statistic of sequences
    num_seq = len(seq_len)
    min_len = min(seq_len)
    max_len = max(seq_len)
    num_nc = sum(seq_len)
    with open('../Results/data_summary.txt', 'w') as fp:
        fp.write(f'Total number of sequences = {num_seq}\n')
        fp.write(f'Minimum length of sequences = {min_len}\n')
        fp.write(f'Maximum length of sequences = {max_len}\n')
        fp.write(f'Total number of base pairs = {num_nc}')
    logging.info('Statistic output is done.')

    # Plot the distribution of sequence length
    fig, ax = plt.subplots()
    ax.set_title('Sequence Length Distribution', fontsize=18)
    ax.set_xlabel('Length', fontsize=16)
    ax.set_ylabel('Fraction of Sequences', fontsize=16)
    ax.set_xscale('log')

    bins = np.logspace(np.log10(min_len), np.log10(max_len), num=20)
    ax.hist(seq_len, bins=bins, density=True, color='C5')

    ax.set_xlim(left=1)
    fig.subplots_adjust(left=0.14, bottom=0.13)
    fig.savefig('../Results/sequence_length_distribution.png', dpi=300)

    return


if __name__ == '__main__':
    main()
