from AlgebraicML.utilities import *
from AlgebraicML.algorithms import *
from AlgebraicML.graph import Graph

#
# Main Program
#

_DEFAULT_TRAINING_DATA_ = "./data/grids.data"
_DEFAULT_TEST_DATA_ = "./data/gridTest.data"

if sys.argv.count == 2:
    inputFile = sys.argv[1]
    testFile = sys.argv[2]
else:
    inputFile = _DEFAULT_TRAINING_DATA_
    testFile = _DEFAULT_TEST_DATA_

labels, records = readData(inputFile)

graph = Graph(labels, records)
enforceNegativeTraceConstraints(graph)
enforcePositiveTraceContraints(graph)
print(verifyTraceConstraints(graph))
cross(graph)
# reduceAtoms(graph)

correct = 0
total = 0

with open(testFile, 'r') as file:
    for line in file:
        data = list(map(lambda x: x.strip(), line.split(",")))
        classified = int(classify(data[:-1], labels, graph))
        trueLabel = int(data[~0])
        total += 1
        correct += int(trueLabel == classified)
        print(data[0], ":", classified, trueLabel == classified)

print("Error Margin:", correct / total) 
obeserveAtoms(graph)          
