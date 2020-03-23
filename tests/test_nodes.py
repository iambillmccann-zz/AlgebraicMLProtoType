from AlgebraicML.nodes import Node

import pytest

@pytest.mark.usefixtures('nodes_for_testing')
class TestNodes:
    def test_the_test_node_framework(self):
        assert True == True

    def test_the_name_is_set(self):
        assert '{}'.format(self.root) == 'root'

    def test_the_initial_name(self):
        expected = 'root'
        actual = self.root.name
        assert actual == expected

    def test_the_initial_parents(self):
        assert len(self.root.parents) == 0

    def test_the_initial_children(self):
        assert len(self.root.children) == 0

    def test_the_initial_dual(self):
        assert len(self.root.dual) == 0

    def test_add_a_child(self):

        number_of_parents = len(self.root.parents)

        self.root.addChild(self.child1)
        assert len(self.root.children) == 1

        self.root.addChild(self.child2)
        assert len(self.root.children) == 2

        # Test that duplicates will not get added!
        self.root.addChild(self.child1)
        assert len(self.root.children) == 2

        # Test that adding children did not affect parents
        assert len(self.root.parents) == number_of_parents

    def test_add_a_parent(self):

        number_of_children = len(self.root.children)

        self.root.addParent(self.parent1)
        assert len(self.root.parents) == 1

        self.root.addParent(self.parent2)
        assert len(self.root.parents) == 2

        # Test that duplicates will not get added!
        self.root.addParent(self.parent2)
        assert len(self.root.parents) == 2

        # Test that adding parents did not affect children
        assert len(self.root.children) == number_of_children

    def test_remove_child(self):

        number_of_parents = len(self.root.parents)

        # Test that you cannot remove a parent from the children
        self.root.removeChild(self.parent1)
        assert len(self.root.children) == 2

        self.root.removeChild(self.child1)
        assert len(self.root.children) == 1

        # Test that you cannot remove a child twice
        self.root.removeChild(self.child1)
        assert len(self.root.children) == 1

        self.root.removeChild(self.child2)
        assert len(self.root.children) == 0

        # Test that you do not crash if there are no children
        self.root.removeChild(self.child2)
        assert len(self.root.children) == 0

        # Test that adding children did not affect parents
        assert len(self.root.parents) == number_of_parents

    def test_remove_parent(self):
    
        number_of_children = len(self.root.children)

        # Test that you cannot remove a parent from the children
        self.root.removeParent(self.child1)
        expected = 2
        actual = len(self.root.parents)
        assert actual == expected

        self.root.removeParent(self.parent1)
        assert len(self.root.parents) == 1

        # Test that you cannot remove a parent twice
        self.root.removeParent(self.parent1)
        assert len(self.root.parents) == 1

        self.root.removeParent(self.parent2)
        assert len(self.root.parents) == 0

        # Test that you do not crash if there are no children
        self.root.removeParent(self.parent2)
        assert len(self.root.parents) == 0

        # Test that adding children did not affect parents
        assert len(self.root.children) == number_of_children

    def test_add_a_dual(self):
    
        self.root.addDual(self.child1)
        assert len(self.root.dual) == 1

        self.root.addDual(self.parent1)
        assert len(self.root.dual) == 2

    def test_add_root_to_itself(self):

        number_of_children = len(self.root.children)
        number_of_parents = len(self.root.parents)
        number_of_duals = len(self.root.dual)

        self.root.addChild(self.root)
        self.root.addParent(self.root)
        self.root.addDual(self.root)

        assert len(self.root.children) == number_of_children + 1
        assert len(self.root.parents) == number_of_parents + 1
        assert len(self.root.dual) == number_of_duals + 1
