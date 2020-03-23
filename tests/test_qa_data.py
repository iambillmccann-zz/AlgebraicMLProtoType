from AlgebraicML.nodes import Node

import pytest

@pytest.mark.usefixtures('training_data')
class TestClass:
    """
    def __init__(self):
        self.test_this = 'Test This'
    """
    def test_the_test_frameword(self):
        assert True == True

    def test_number_of_labels(self):
        actual = len(self.labels)
        expected = 4
        assert actual == 4

    def test_content_of_labels(self):
        expected = ["Name", "Eye_Color", "Hair_Color", "Label"]
        assert self.labels == expected

    def test_number_of_records(self):
        expected = 4
        assert len(self.records) == expected

    def test_number_of_columns(self):
        expected = 4
        assert len(self.records[1]) == expected