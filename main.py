from AlgebraicML.utilities import *
from AlgebraicML.algorithms import *
from AlgebraicML.graph import Graph

#
# Main Program
#

_DEFAULT_TRAINING_DATA_ = "./data/people.data"
_DEFAULT_TEST_DATA_ = "./data/test.data"

labels, records = readData(_DEFAULT_TRAINING_DATA_)
graph = Graph(labels, records)
enforceNegativeTraceConstraints(graph)
enforcePositiveTraceContraints(graph)
cross(graph)
#reduceAtoms(graph)

with open(_DEFAULT_TEST_DATA_, 'r') as file:
    for line in file:
        data = list(map(lambda x: x.strip(), line.split(",")))
        print(data[0], ":", classify(data, labels, graph))

#print(classify(["Steve","blue","brown"], labels, layers))          
    