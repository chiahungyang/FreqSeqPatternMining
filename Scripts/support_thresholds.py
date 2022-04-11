"""Obtain support thresholds of significance for Arabdopsis thiliana data."""


from readfasta import sequences
from sigthrd import thresholds
import matplotlib.pyplot as plt
import logging

plt.style.use('ggplot')
logging.basicConfig(format='%(asctime)s - %(message)s', level=logging.INFO)

# Global variables
MIN_LEN = 4  # Minimum pattern length of interests
LWR_BD = 100  # Lower bound of support thresholds
CONF = 0.9  # Confidence level of significance


def main():
    """Obtain and plot the support thresholds. of significance."""
    # Extract total number of nucleotides from the data
    logging.info('Start.')
    num_nc = sum([len(seq) for seq in sequences('../Data/EpiR_all_seq.fasta')])
    logging.info('Data loaded.')

    # Obtain support thresholds of significance
    len_thrd = thresholds(num_nc, MIN_LEN, LWR_BD, CONF)
    _len, thrd = zip(*len_thrd.items())
    logging.info('Support thresholds computed.')

    # Plot the support thresholds
    fig, ax = plt.subplots()
    ax.set_title('Support Thresholds of Significance', fontsize=18)
    ax.set_xlabel('Sequential Pattern Length', fontsize=16)
    ax.set_ylabel('Support Threshold', fontsize=16)
    ax.set_yscale('log')
    ax.bar(_len, thrd, color='C2', width=0.6)
    ax.set_ylim(bottom=LWR_BD)
    fig.savefig('../Results/support_thresholds.png', dpi=300)

    return


if __name__ == '__main__':
    main()
