from AlgebraicML.utilities import *
from AlgebraicML.nodes import Node

import pytest
"""
When pytest discovers a conftest.py, it modifies sys.path so it can import stuff from the conftest module. 
So, pytest will be forced to append it to sys.path. The fixture is used to setup a global test environment
for all the training_data classes.
"""
@pytest.fixture(scope = 'class')
def training_data(request):

    labels, records = readData('./data/people.data')
    request.cls.labels = labels
    request.cls.records = records
    yield

@pytest.fixture(scope = 'class')
def nodes_for_testing(request):

    root = Node(name = "root")
    child1 = Node(name = "child1")
    child2 = Node(name = "child2")
    parent1 = Node(name = "parent1")
    parent2 = Node(name = "parent2")

    request.cls.root = root
    request.cls.child1 = child1
    request.cls.child2 = child2
    request.cls.parent1 = parent1
    request.cls.parent2 = parent2
    yield