"""Test the A-priori, position-based and hybrid method in freqsubseq.py. """


import unittest
from freqsubseq import initialized_tree, candidates, run_apriori
from freqsubseq import assign_pattern_index, initialized_queue
from freqsubseq import candidate_mappings, run_position
from freqsubseq import number_of_freq_subseq
from freqsubseq import frequent_patterns, Counter


class AuxReader:
    """Auxiliary reader to set up test cases."""

    def __init__(self, seqs):
        """Initiate reader with a collection of sequences."""
        self.seqs = seqs

    def items(self):
        """Return a generator of sequences that are stored."""
        for seq in self.seqs:
            yield seq


class AlgorithmTestCase(unittest.TestCase):
    def setUp(self):
        seq = 'ACGATTCGATCG'
        num = 2
        self.data = AuxReader([seq for _ in range(num)])
        self.thrd = 4
        self.min_len = 1
        self.max_len = 3
        self.ground_truth = {('A', 6), ('C', 6), ('G', 6), ('T', 6), ('AT', 4),
                             ('CG', 6), ('GA', 4), ('TC', 4), ('CGA', 4),
                             ('GAT', 4), ('TCG', 4)}

    def test_initialized_tree(self):
        m = 2
        tree = initialized_tree(self.data, m)

        # Check the edgelist of the tree
        el_true = {('', 'A'), ('', 'C'), ('', 'G'), ('', 'T'), ('A', 'AC'),
                   ('A', 'AT'), ('C', 'CG'), ('G', 'GA'), ('T', 'TT'),
                   ('T', 'TC')}
        el = {(parent.seq, child.seq)
              for parent in tree.all_nodes()
              for child in parent.all_children()}
        self.assertEqual(el, el_true)

        # Check the counts for leaves
        counts_true = {('AC', 2), ('AT', 4), ('CG', 6), ('GA', 4), ('TT', 2),
                       ('TC', 4)}
        counts = {(node.seq, node.count) for node in tree.nodes_at_level(m)}
        self.assertEqual(counts, counts_true)

    def test_candidates(self):
        m = 2
        tree = initialized_tree(self.data, m)
        cand_true = {'ATC', 'CGA', 'GAT', 'TCG'}
        cand = set(candidates(tree, m+1, self.thrd))
        self.assertEqual(cand, cand_true)

    def test_run_apriori(self):
        m = 2
        tree = initialized_tree(self.data, m)
        run_apriori(self.data, tree, m+1, self.thrd)

        # Check the edgelist of the tree
        el_true = {('', 'A'), ('', 'C'), ('', 'G'), ('', 'T'), ('A', 'AC'),
                   ('A', 'AT'), ('C', 'CG'), ('G', 'GA'), ('T', 'TT'),
                   ('T', 'TC'), ('AT', 'ATC'), ('CG', 'CGA'), ('GA', 'GAT'),
                   ('TC', 'TCG')}
        el = {(parent.seq, child.seq)
              for parent in tree.all_nodes()
              for child in parent.all_children()}
        self.assertEqual(el, el_true)

        # Check the counts for the candidate patterns
        counts_true = {('ATC', 2), ('CGA', 4), ('GAT', 4), ('TCG', 4)}
        counts = {(node.seq, node.count) for node in tree.nodes_at_level(m+1)}
        self.assertEqual(counts, counts_true)

    def test_assign_pattern_index(self):
        m = 2
        tree = initialized_tree(self.data, m)
        assign_pattern_index(tree, m, self.thrd)
        has_idx_ture = {('AC', False), ('AT', True), ('CG', True),
                        ('GA', True), ('TT', False), ('TC', True)}
        has_idx = {(node.seq, hasattr(node, 'idx'))
                   for node in tree.nodes_at_level(m)}
        self.assertEqual(has_idx, has_idx_ture)

    def test_initialized_queue(self):
        m = 2
        tree = initialized_tree(self.data, m)
        assign_pattern_index(tree, m, self.thrd)
        idx = {seq: tree.pattern(seq).idx
               for seq in ['AT', 'CG', 'GA', 'TC']}
        queue_true = [(0, 1, idx['CG']), (0, 2, idx['GA']), (0, 3, idx['AT']),
                      (0, 5, idx['TC']), (0, 6, idx['CG']), (0, 7, idx['GA']),
                      (0, 8, idx['AT']), (0, 9, idx['TC']), (0, 10, idx['CG'])]
        queue_true += [(1, 1, idx['CG']), (1, 2, idx['GA']), (1, 3, idx['AT']),
                       (1, 5, idx['TC']), (1, 6, idx['CG']), (1, 7, idx['GA']),
                       (1, 8, idx['AT']), (1, 9, idx['TC']), (1, 10, idx['CG'])]
        queue = initialized_queue(self.data, tree, m, self.thrd)
        self.assertEqual(queue, queue_true)

    def test_cadidate_mappings(self):
        m = 2
        tree = initialized_tree(self.data, m)
        assign_pattern_index(tree, m, self.thrd)
        idx = {seq: tree.pattern(seq).idx
               for seq in ['AT', 'CG', 'GA', 'TC']}

        idx_set = {0, 1, 2, 3}
        mapping_true = {('ATC', (idx['AT'], idx['TC'])),
                        ('CGA', (idx['CG'], idx['GA'])),
                        ('GAT', (idx['GA'], idx['AT'])),
                        ('TCG', (idx['TC'], idx['CG']))}
        mapping = set()
        for seq, indices, cand_idx in candidate_mappings(tree, m+1, self.thrd):
            idx_set.remove(cand_idx)
            mapping.add((seq, indices))
        self.assertEqual(len(idx_set), 0)
        self.assertEqual(mapping, mapping_true)

    def test_run_position(self):
        m = 2
        tree = initialized_tree(self.data, m)
        assign_pattern_index(tree, m, self.thrd)
        queue = initialized_queue(self.data, tree, m, self.thrd)
        run_position(queue, tree, m+1, self.thrd)
        idx = {seq: tree.pattern(seq).idx
               for seq in ['ATC', 'CGA', 'GAT', 'TCG']}

        # Check queue
        queue_true = [(0, 1, idx['CGA']), (0, 2, idx['GAT']),
                      (0, 5, idx['TCG']), (0, 6, idx['CGA']),
                      (0, 7, idx['GAT']), (0, 8, idx['ATC']), (0, 9, idx['TCG'])]
        queue_true += [(1, 1, idx['CGA']), (1, 2, idx['GAT']),
                       (1, 5, idx['TCG']), (1, 6, idx['CGA']),
                       (1, 7, idx['GAT']), (1, 8, idx['ATC']), (1, 9, idx['TCG'])]
        self.assertEqual(queue, queue_true)

        # Check counts in the tree
        counts_true = {('ATC', 2), ('CGA', 4), ('GAT', 4), ('TCG', 4)}
        counts = {(node.seq, node.count) for node in tree.nodes_at_level(m+1)}
        self.assertEqual(counts, counts_true)

    def test_number_of_freq_subseq(self):
        m = 1
        tree = initialized_tree(self.data, m)
        num = number_of_freq_subseq(tree, m, self.thrd)
        self.assertEqual(num, 24)

    def test_frequent_patterns(self):
        # A-priori
        counter = Counter()
        results = frequent_patterns(self.data, self.min_len, self.max_len,
                                    self.thrd, method='apriori',
                                    counter=counter)
        self.assertEqual(set(results), self.ground_truth)
        self.assertEqual(counter.count, 3)

        # Position-based method
        counter = Counter()
        results = frequent_patterns(self.data, self.min_len, self.max_len,
                                    self.thrd, method='position',
                                    counter=counter)
        self.assertEqual(set(results), self.ground_truth)
        self.assertEqual(counter.count, 2)

        # Hybrid method
        counter = Counter()
        results = frequent_patterns(self.data, self.min_len, self.max_len,
                                    self.thrd, method='hybrid',
                                    counter=counter)
        self.assertEqual(set(results), self.ground_truth)
        self.assertEqual(counter.count, 2)


if __name__ == '__main__':
    unittest.main()
