"""To read nucleotide sequence data from fasta file."""


from Bio import SeqIO


def sequences(filepath):
    """
    Return a generator of sequences (in strings) in the data.

    Parameter
    ---------
    filepath: str
        File path to the data, which must be a fasta file.

    """
    for record in SeqIO.parse(filepath, 'fasta'):
        yield str(record.seq).upper()


def samples(filepath, seq_IDs):
    """
    Return a generator of sampled sequences (in strings) from the data.

    Parameters
    ----------
    filepath: str
        File path to the data, which must be a fasta file.

    seq_IDs: list
        IDs of sequences in the data to be sampled, which has been sorted.

    """
    k = len(seq_IDs)
    idx = 0
    for i, record in enumerate(SeqIO.parse(filepath, 'fasta')):
        if i == seq_IDs[idx]:
            yield str(record.seq).upper()

            if idx < k - 1:
                idx += 1
            else:
                break


class Reader:
    """Data reader of the fasta file."""

    def __init__(self, filepath, sample=None):
        """
        Initiate reader with file path to the data and additional information.

        Parameters
        ----------
        filepath: str
            File path to the data, which must be a fasta file.

        sample: list of None
            A list of sampled sequences' ID to be read. Default to None.

        """
        self.filepath = filepath
        self.sample = sample

    def items(self):
        """Return a generator of items in the data."""
        if self.sample is None:
            return sequences(self.filepath)
        else:
            return samples(self.filepath, self.sample)


def main():
    """Empty main function."""
    return


if __name__ == '__main__':
    main()
