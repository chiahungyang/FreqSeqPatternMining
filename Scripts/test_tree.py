"""Test the Tree class used in algorithms find frequent sequential patterns."""


import unittest
from freqsubseq import Tree, PatternNotFound


class TreeTestCase(unittest.TestCase):
    def setUp(self):
        # Manually built Tree
        self.root = Tree('This is a root.')
        self.root.A = Tree('This is a child.')

        # Tree built by add method
        seq = 'ACGATTCGATCG'
        m = 3
        self.tree = Tree()
        for i in range(len(seq)-m+1):
            self.tree.add(seq[i:i+m], count=1)
        self.tree.add('ACG', attr_val={'idx': 0})

    def test_init(self):
        self.assertIsInstance(self.root, Tree)
        self.assertEqual(self.root.seq, 'This is a root.')

    def test_has_child(self):
        self.assertTrue(self.root.has_child('A'))

    def test_child(self):
        self.assertIsInstance(self.root.child('A'), Tree)
        self.assertEqual(self.root.child('A').seq, 'This is a child.')
        self.assertRaises(PatternNotFound, self.root.child, 'B')

    def test_pattern(self):
        self.assertIsInstance(self.tree.pattern('ACG'), Tree)
        self.assertEqual(self.tree.pattern('ACG').seq, 'ACG')
        self.assertRaises(PatternNotFound, self.tree.pattern, 'ACA')

    def test_increment_count(self):
        old_count = self.tree.pattern('ACG').count
        self.tree.increment_count('ACG', 1)
        new_count = self.tree.pattern('ACG').count
        self.assertEqual(new_count, old_count + 1)
        self.tree.increment_count('ACG', -1)
        new_count = self.tree.pattern('ACG').count
        self.assertEqual(new_count, old_count)

    def test_all_children(self):
        are_trees = all([isinstance(node, Tree)
                         for node in self.tree.all_children()])
        self.assertTrue(are_trees)
        children_true = {'A', 'T', 'C', 'G'}
        children = {node.seq for node in self.tree.all_children()}
        self.assertEqual(children, children_true)

    def test_all_nodes(self):
        are_trees = all([isinstance(node, Tree)
                         for node in self.tree.all_nodes()])
        self.assertTrue(are_trees)
        nodes_true = {'', 'A', 'C', 'G', 'T', 'AC', 'AT', 'CG', 'GA', 'TT',
                      'TC', 'ACG', 'ATT', 'ATC', 'CGA', 'GAT', 'TTC', 'TCG'}
        nodes = {node.seq for node in self.tree.all_nodes()}
        self.assertEqual(nodes, nodes_true)

    def test_nodes_at_level(self):
        m = 3
        are_trees = all([isinstance(node, Tree)
                         for node in self.tree.nodes_at_level(m)])
        self.assertTrue(are_trees)
        leaves_true = {'ACG', 'ATT', 'ATC', 'CGA', 'GAT', 'TTC', 'TCG'}
        leaves = {node.seq for node in self.tree.nodes_at_level(m)}
        self.assertEqual(leaves, leaves_true)

    def test_add(self):
        # Check the edgelist of tree
        el_true = {('', 'A'), ('', 'C'), ('', 'G'), ('', 'T'), ('A', 'AC'),
                   ('A', 'AT'), ('C', 'CG'), ('G', 'GA'), ('T', 'TT'),
                   ('T', 'TC'), ('AC', 'ACG'), ('AT', 'ATT'), ('AT', 'ATC'),
                   ('CG', 'CGA'), ('GA', 'GAT'), ('TT', 'TTC'), ('TC', 'TCG')}
        el = {(parent.seq, child.seq)
              for parent in self.tree.all_nodes()
              for child in parent.all_children()}
        self.assertEqual(el, el_true)

        # Check whether additional information was properly added
        node = self.tree.pattern('ACG')
        self.assertTrue(hasattr(node, 'idx'))
        self.assertEqual(node.idx, 0)

        # Check whether the counts are correct
        m = 3
        has_count = all([hasattr(node, 'count')
                         for node in self.tree.nodes_at_level(m)])
        self.assertTrue(has_count)
        counts_true = {('ACG', 1), ('ATT', 1), ('ATC', 1), ('CGA', 2),
                       ('GAT', 2), ('TTC', 1), ('TCG', 2)}
        counts = {(node.seq, node.count)
                  for node in self.tree.nodes_at_level(m)}
        self.assertEqual(counts, counts_true)

    def test_has_count_at_least(self):
        self.assertTrue(self.tree.pattern('ACG').has_count_at_least(0))
        self.assertTrue(self.tree.pattern('ACG').has_count_at_least(1))
        self.assertFalse(self.tree.pattern('ACG').has_count_at_least(2))
        self.assertFalse(self.tree.has_count_at_least(0))

    def test_all_patterns_at_least(self):
        thrd = 2
        pairs_true = {('CGA', 2), ('GAT', 2), ('TCG', 2)}
        pairs = set(self.tree.all_patterns_at_least(thrd))
        self.assertEqual(pairs, pairs_true)


if __name__ == '__main__':
    unittest.main()
