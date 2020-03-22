from AlgebraicML import nodes

def readData(file):
    """
    """
    with open(file, 'r') as data:
        labels = list(map(lambda x: x.strip(), data.readline().split(",")))
        records = []

        for line in data:
            records.append(list(map(lambda x: x.strip(),line.split(","))))
    
    return labels, records
