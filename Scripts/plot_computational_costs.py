"""Plot the computational costs of frequent pattern finding methods."""


import numpy as np
import json
import matplotlib.pyplot as plt

plt.style.use('ggplot')


def binning(x, y, bins):
    """
    Return the binned function of the data.

    Parameters
    ----------
    x: list
        Data to be binned.

    y: list
        Values of a function of the data, which must have the same length as x.

    bins: numpy array
        Pre-specified bins to be used.

    Returns
    -------
    centers: list
        Centers of each bins

    means: list
        Mean of the binned function values.

    stderrs: list
        Standard deviation of the binned function values.

    """
    data = np.array(x)
    values = np.array(y)
    num = len(bins) - 1
    binned = [list() for _ in range(num)]

    # Categorized the data into bins
    # B_0 = [b_0, b_1]; B_i = (b_i, b_i+1]
    for j in range(len(data)):
        i = np.argmax(bins[1:] >= data[j])
        binned[i].append(values[j])

    # Obtain the center of each bin
    centers = bins[:-1] + 0.5 * np.diff(bins)

    # Obtain the means and standard errors of binned function values
    means = np.array([np.mean(b) for b in binned if b])
    stderrs = np.array([np.std(b) for b in binned if b])

    return centers, means, stderrs


def main():
    """
    Plot the number of passes along with the number of base pairs in the input.

    """
    with open('../Results/computational_costs.json', 'r') as fp:
        data = json.load(fp)
        _min = min([n for _, n, _ in data])
        _max = max([n for _, n, _ in data])
        bins = np.logspace(np.log10(_min), np.log10(_max), 10 + 1)

    fig, ax = plt.subplots()
    ax.set_xlabel('Number of Base Pairs', fontsize=16)
    ax.set_ylabel('Number of Passes', fontsize=16)

    # A-priori
    pairs = [(n, c) for method, n, c in data if method == 'apriori']
    num_nc, _pass = zip(*pairs)
    cent, mean, stdr = binning(num_nc, _pass, bins)
    ax.errorbar(cent, mean, yerr=stdr, fmt='o', color='C1', ecolor='C1',
                label='A-priori')

    # Position-based
    pairs = [(n, c) for method, n, c in data if method == 'position']
    num_nc, _pass = zip(*pairs)
    cent, mean, stdr = binning(num_nc, _pass, bins)
    ax.errorbar(cent, mean, yerr=stdr, fmt='o', color='C0', ecolor='C0',
                label='Position-based')

    ax.set_xscale('log')
    ax.legend()
    fig.savefig('../Results/computational_costs.png', dpi=300)

    return


if __name__ == '__main__':
    main()
