"""Algorithms for mining frequent contiguous sequential patterns in genomes."""


import logging

logging.basicConfig(format='%(asctime)s - %(message)s', level=logging.INFO)

# Global variable
NUCLEOTIDES = ('A', 'C', 'T', 'G')
MEM_QUEUE = 8 * 10**9  # Amount of memory for queue
MEM_PER_ELEM = 80  # Amount of memory used for every element in the queue
MAX_SIZE = int(round(MEM_QUEUE / MEM_PER_ELEM))  # Maximum size of the queue


class PatternNotFound(Exception):
    """Exception that a sequential pattern has not been stored."""
    pass


class ExceedAllocatedMemory(Exception):
    """Exception that the queue exceeds the amount of allocated memory."""
    pass


class Tree:
    """
    Root of a (sub)tree that store information of sequential patterns, where
    a node representing a length m pattern has the parent as its length (m-1)
    prefix.

    """

    def __init__(self, seq=''):
        """
        Initiate the root of a (sub)tree.

        Parameter
        ---------
        seq: str
            Sequential pattern representation of the node. Default to ''.
        """
        self.seq = seq

    def has_child(self, nc):
        """
        Return whether the node has a child that append the pattern by a given
        nucleotide.

        Parameter
        ---------
        nc: string
            Last nucleotide of child pattern, if any.

        """
        return hasattr(self, nc)

    def child(self, nc):
        """
        Return the child node that append the pattern by a given nucleotide. A
        PatternNotFound exception is raised if no such a node exists.

        Parameter
        ---------
        nc: string
            Last nucleotide of child pattern, if any.

        """
        if self.has_child(nc):
            return getattr(self, nc)
        else:
            raise PatternNotFound

    def add(self, seq, count=None, attr_val=None):
        """
        Extend the (sub)tree if a non-exsisting sequential pattern suffix is
        added.

        Parameters
        ----------
        seq: str
            Sequential pattern suffix to be added.

        count: int or None
            If not None, increment the count of the target pattern by given
            integer or initiate its attribute of count.

        attr_val: dict or None
            Dictionary mapping from attribute names to corresponding values. If
            not None, initiate the attributes at the target node with given
            values.

        """
        # If the pattern suffix is empty, update the current node
        if len(seq) == 0:
            if count is not None:
                self.count = self.count + count if hasattr(self, 'count') else count
            if attr_val is not None:
                for k, v in attr_val.items():
                    setattr(self, k, v)
        # Otherwise, add suffix to the appropriate child
        else:
            nc = seq[0]
            # Build the appropriate if one does not exist
            if not hasattr(self, nc):
                setattr(self, nc, Tree(self.seq + nc))
            self.child(nc).add(seq[1:], count, attr_val)

    def pattern(self, seq):
        """
        Return the node corresponding to the given pattern suffix. A
        PatternNotFound exception is raised if no such a node exists.

        seq: str
            Sequential pattern suffix to search for.

        """
        # If the pattern suffix is empty, return current node.
        if len(seq) == 0:
            return self
        # Otherwise, search for an appropriate child
        else:
            nc = seq[0]
            if not self.has_child(nc):
                raise PatternNotFound
            else:
                return self.child(nc).pattern(seq[1:])

    def increment_count(self, seq, count):
        """
        Increment the count of a given pattern suffix, if any exists in the
        (sub)tree.

        Parameters
        ----------
        seq: str
            Sequential pattern suffix to increment count.

        count: int
            Amount of count increment.

        """
        try:
            node = self.pattern(seq)
            node.count = node.count + count if hasattr(node, 'count') else count
        except PatternNotFound:
            pass

    def has_count_at_least(self, thrd):
        """Retrun whether the node has count not less than the threshold."""
        return hasattr(self, 'count') and (self.count >= thrd)

    def all_children(self):
        """Return a generator of all children of the current node."""
        for nc in NUCLEOTIDES:
            if self.has_child(nc):
                yield self.child(nc)

    def all_nodes(self):
        """Return a generator of all nodes in the (sub)tree."""
        # Yield itself
        yield self
        # and traverse all nodes in its children's subtree
        for child in self.all_children():
            for node in child.all_nodes():
                yield node

    def all_patterns_at_least(self, thrd):
        """
        Return a generator of sequential patterns and corresponding supports
        such that the supports are not less than a given threshold.

        """
        for node in self.all_nodes():
            if node.has_count_at_least(thrd):
                yield (node.seq, node.count)

    def nodes_at_level(self, level):
        """
        Return a generator of all nodes in the (sub)tree at a given level.

        Parameter
        ---------
        level: int
            Level at which nodes are extracted.

        """
        # If at level 0, yield the current node
        if level == 0:
            yield self
        # Otherwise, search in its children's subtree
        else:
            for child in self.all_children():
                for node in child.nodes_at_level(level - 1):
                    yield node


class Counter:
    """Counter of number of passes through the data."""

    def __init__(self, count=0):
        """Initiate counter."""
        self.count = count

    def increment(self, count=1):
        """Increment the recorded number of passes."""
        self.count += count


def increment(counter):
    """Increment the counter if any."""
    if counter is not None:
        counter.increment()


def subsequences(seq, _len):
    """
    Return a generator of all subsequences with given length in the given
    sequence.

    """
    n = len(seq)
    for i in range(n - _len + 1):
        yield seq[i:i+_len]


def initialized_tree(data, _len):
    """
    Return a Tree that stores sequential patterns of given length and their
    support after running the first pass through the data.

    Parameters
    ----------
    data: readfasta Reader
        Data reader along with a items method to read elements in the data.

    _len: int
        Length of target sequential patterns in the first pass.

    """
    tree = Tree()
    for seq in data.items():
        for subseq in subsequences(seq, _len):
            tree.add(subseq, count=1)

    return tree


def candidates(tree, _len, thrd):
    """
    Return a generator of candidate sequential patterns of given length that
    can be deduced from the prefix and suffix patterns with support not less
    than the given threshold.

    """
    # For all possible prefix patterns with support not less than the threshold
    for node_p in tree.nodes_at_level(_len - 1):
        if node_p.has_count_at_least(thrd):
            prefix = node_p.seq
            # Examine existence of the node representing intersection between
            # the prefix and suffix
            try:
                node_i = tree.pattern(prefix[1:])
            except PatternNotFound:
                continue
            # Seach for suffixes also with support not less than the threshold
            for node_s in node_i.all_children():
                if node_s.has_count_at_least(thrd):
                    suffix = node_s.seq
                    yield prefix + suffix[-1:]


def run_apriori(data, tree, _len, thrd):
    """
    Run a single iteration of the A-priori algorithm and update the tree with
    sequential patterns of given length with supports not less than the
    threshold.

    Parameters
    ----------
    data: readfasta Reader
        Data reader along with an items method to read elements in the data.

    tree: Tree
        Tree object storing sequential patterns.

    _len: int
        Length of target patterns.

    thrd: int
        Support threshold of sequential patterns.

    """
    # Add candidate patterns to the tree
    for seq in candidates(tree, _len, thrd):
        tree.add(seq)
    # Pass the data once and obtain supports of the candidates
    for seq in data.items():
        for subseq in subsequences(seq, _len):
            tree.increment_count(subseq, count=1)


def assign_pattern_index(tree, _len, thrd):
    """
    Add an index attribute to patterns of given length in the tree, whose
    supports are not less than the threshold.

    """
    idx = 0
    for node in tree.nodes_at_level(_len):
        if node.has_count_at_least(thrd):
            node.idx = idx
            idx += 1


def initialized_queue(data, tree, _len, thrd):
    """
    Return a queue that stores position information of frequent subsequences
    of given length. Frequent sequential patterns in the tree must have the
    attribute 'idx'.

    Parameters
    ----------
    data: readfasta Reader
        Data reader along with an items method to read elements in the data.

    tree: Tree
        Tree object storing sequential patterns.

    _len: int
        Length of target patterns.

    thrd: int
        Support threshold of sequential patterns.

    Return
    ------
    queue: list
        List of frequent subsequences of gien length that have beendiscovered,
        each of whose elements is a tuple of sequence ID, first position of the
        subsequence and the sequential pattern index. The ordering of
        subsequences remains the same as they were read from the data.

    """
    queue = list()
    for _id, seq in enumerate(data.items()):
        for pos, subseq in enumerate(subsequences(seq, _len)):
            try:
                node = tree.pattern(subseq)
                if node.has_count_at_least(thrd):
                    queue.append((_id, pos, node.idx))
            except PatternNotFound:
                continue

    return queue


def candidate_mappings(tree, _len, thrd):
    """
    Return a generator of candidate patterns and corresponding mappings from
    indices of prefix and suffix patterns to the candidate's index. All prefix
    (suffix) patterns must have the attribute 'idx'.

    """
    idx = 0  # Index of candidate pattern
    # For all possible prefix patterns with support not less than the threshold
    for node_p in tree.nodes_at_level(_len - 1):
        if node_p.has_count_at_least(thrd):
            prefix = node_p.seq
            # Examine existence of the node representing intersection between
            # the prefix and suffix
            try:
                node_i = tree.pattern(prefix[1:])
            except PatternNotFound:
                continue
            # Seach for suffixes also with support not less than the threshold
            for node_s in node_i.all_children():
                if node_s.has_count_at_least(thrd):
                    suffix = node_s.seq
                    seq = prefix + suffix[-1:]
                    yield (seq, (node_p.idx, node_s.idx), idx)
                    idx += 1


def run_position(queue, tree, _len, thrd):
    """
    Run a single iteration of the position-based algorithm and update the tree
    with sequential patterns of given length with supports not less than the
    threshold.

    Parameters
    ----------
    queue: list
        List of frequent subsequences discovered in the previous iteration,
        each of whose elements is a tuple of sequence ID, first position of the
        subsequence and the sequential pattern index.

    tree: Tree
        Tree object storing sequential patterns.

    _len: int
        Length of target patterns.

    thrd: int
        Support threshold of sequential patterns.

    Note
    ----
    Invariance: The ordering of subsequences in the queue remains the same as
                they were read from the data.

    """
    # If the queue is empty, do nothing
    if len(queue) == 0:
        return

    # Pre-computation to:
    # 1. Add candidate patterns to the tree with indices
    # 2. Build a mapping from prefix and suffix indices to the candidate index
    # 3. Build a dictionary mapping from candidate index to sequential pattern
    index_map = dict()
    pattern = dict()
    for seq, indices, cand_idx in candidate_mappings(tree, _len, thrd):
        tree.add(seq, attr_val={'idx': cand_idx})
        index_map[indices] = cand_idx
        pattern[cand_idx] = seq

    # Search for occurence of candidate patterns over adjacent subsequences in
    # the queue from the beginning, and increment supports in the tree
    new_queue = list()
    current_ID, current_pos, current_idx = queue.pop(0)
    while len(queue) > 0:
        next_ID, next_pos, next_idx = queue.pop(0)
        if (next_ID == current_ID) and (next_pos - current_pos == 1):
            # If the adjacent patterns are prefix and suffix of a candidate
            if (current_idx, next_idx) in index_map:
                cand_idx = index_map[(current_idx, next_idx)]
                new_queue.append((current_ID, current_pos, cand_idx))
                tree.increment_count(pattern[cand_idx], count=1)
        current_ID, current_pos, current_idx = next_ID, next_pos, next_idx

    # Update the queue to frequent subsequences discovered in the current
    # iteration
    queue.extend(new_queue)


def number_of_freq_subseq(tree, _len, thrd):
    """
    Return the total number of frequent subsequences dicovered in the tree,
    which is of the given length.

    """
    total = 0
    for node in tree.nodes_at_level(_len):
        if node.has_count_at_least(thrd):
            total += node.count

    return total


def frequent_patterns(data, min_len, max_len, thrd, counter=None,
                      method='hybrid', verbose=False):
    """
    Return a generator of sequential patterns and supports, which have length
    not less than a given minimum and support not less than a given threshold,
    using an A-priori algorithm.

    Parameters
    ----------
    data: readfasta Reader
        Data reader along with an items method to read elements in the data.

    min_len: int
        Minimum length of sequential patterns of interests.

    max_len: int
        Maximum length of sequential patterns of interests.

    thrd: int
        Support threshold of sequential patterns.

    counter: Counter
        Counter object to keep track of the number of passes through the data.

    method: str
        Algorithm to be used to obtain frequent sequential patterns. Must be
        one of {'apriori', 'position', 'hybrid'}. Default to 'hybrid'.

    verbose: bool
        Whether to log the process.

    """
    # Build a Tree to store sequential patterns
    _len = min_len
    tree = initialized_tree(data, _len)
    increment(counter)

    if method == 'position':
        # If the size of the queue to be constructed does not fit in the
        # allocated memory, raise an ExceedAllocatedMemory exception.
        num = number_of_freq_subseq(tree, _len, thrd)
        if num > MAX_SIZE:
            raise ExceedAllocatedMemory
        # Otherwise build the queue to store position information
        else:
            assign_pattern_index(tree, _len, thrd)
            queue = initialized_queue(data, tree, _len, thrd)
            increment(counter)

    elif method == 'hybrid':
        use_apriori = True

    if verbose:
        num = number_of_freq_subseq(tree, _len, thrd)
        logging.info('Initialization is done. {} instances found.'.format(num))

    # Iterate over all pattern length of interests
    while _len < max_len:
        if method == 'hybrid':
            # If a queue can fit in the allocated memory, use the
            # position-based; otherwise, use the A-priori algorithm
            if use_apriori:
                if number_of_freq_subseq(tree, _len, thrd) <= MAX_SIZE:
                    use_apriori = False
                    assign_pattern_index(tree, _len, thrd)
                    queue = initialized_queue(data, tree, _len, thrd)
                    increment(counter)
        _len += 1

        if method == 'apriori':
            run_apriori(data, tree, _len, thrd)
            increment(counter)
        elif method == 'position':
            run_position(queue, tree, _len, thrd)
        elif method == 'hybrid':
            if use_apriori:
                run_apriori(data, tree, _len, thrd)
                increment(counter)
            else:
                run_position(queue, tree, _len, thrd)

        if verbose:
            num = number_of_freq_subseq(tree, _len, thrd)
            logging.info('An iteration is done. {} instances found.'.format(num))

    return tree.all_patterns_at_least(thrd)


def main():
    """Empty main function."""
    return


if __name__ == '__main__':
    main()
