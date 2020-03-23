from AlgebraicML.nodes import Node

import pytest

@pytest.mark.usefixtures('training_data')
class TestClass:
    """
    def __init__(self):
        self.test_this = 'Test This'
    """

    def test_nodes_setup(self):
        assert self.test_this == 'Test this!'

    def test_the_test_frameword(self):
        assert True == True
