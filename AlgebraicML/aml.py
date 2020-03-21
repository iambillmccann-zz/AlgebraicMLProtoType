

class Node(object):   
    def __init__(self, name=""):
        self.children = set()
        self.parents = set()
        self.name = name
        self.dual = None  
    def addChild(self, node):
        self.children.add(node) 
    def removeChild(self, node):
        self.children.discard(node)  
    def addParent(self, node):
        self.parents.add(node)    
    def removeParent(self, node):
        self.parents.discard(node)   
    def addDual(self, node):
        self.dual.add(node)
    def __str__(self):
        return self.name

def linkNodes(child, parent):
    child.addParent(parent)
    parent.addChild(child)

def linkDuals(node, dual):
    node.dual = dual
    dual.dual = node

def fullLink(child, parent):
    linkNodes(child, parent)
    linkNodes(parent.dual, child.dual)

def deleteNode(node, layers, layer, dlayer):
    for parent in node.parents:
        parent.removeChild(node)
        parent.dual.removeParent(node)
    for child in node.children:
        child.removeParent(node)
        child.dual.removeChild(node)
    
    layers[layer].discard(node)
    if dlayer:
        layers[dlayer].discard(node.dual)

    for child in node.children:
        for parent in node.parents:
            fullLink(child, parent)
    
    return layers

def getLower(node):
    lowerset = set([node])

    for child in node.children:
        lowerset = lowerset.union(getLower(child))
    
    return lowerset

def getConstrainedLower(node, layer):
    return getLower(node).intersection(layer)

def getUpper(node):
    upperset = set()

    for parent in node.parents:
        upperset = upperset.union(getLower(parent))
    
    return upperset

def getConstrainedUpper(node, layer):
    return getUpper(node).intersection(layer)

def getTrace(node, layers):
    atoms = getConstrainedLower(node, layers["atoms"])
    trace = layers["dualAtoms"]
    for atom in atoms:
        trace = trace.intersection(getLower(atom.dual))
    return trace

def traceConstraint(a, b, layers):
    return getTrace(b, layers).issubset(getTrace(a, layers))

def findStronglyDiscriminantCoefficient(a, b, graph):
    Omega = set(map(lambda c: c.dual, getConstrainedLower(a, graph["constants"])))
    U = getTrace(b, layers)
    while U:
        zeta = U.pop()
        T = Omega.difference(getUpper(zeta))
        if T:
            c = T.pop()
            return c.dual
    return None

# Algorithm 1

def enforceNegativeTraceConstraints(graph):
    # Preprocessing Step
    count = 1
    for neg in graph["nterms"]:
        zeta = Node(name="zeta"+str(count))
        linkNodes(zeta, neg.dual)
        graph["dualAtoms"].add(zeta)
        count += 1

    # Main Algorithm
    for neg in graph["nterms"]:
        if not traceConstraint(graph["$"], neg, graph):
            continue
        
        c = None
        while not c:
            c = findStronglyDiscriminantCoefficient(graph["$"], neg, graph)
            if not c:
                h = getConstrainedLower(neg.dual, graph["terms*"]).difference(getLower(graph["$"].dual)).pop()
                zeta = Node(name="zeta"+str(count))
                linkNodes(zeta, h.dual)
                graph["dualAtoms"].add(zeta)
                count += 1
        
        phi = Node(name="phi")
        dphi = Node(name="[phi]")
        linkDuals(phi, dphi)
        fullLink(phi, c)
        graph["atoms"].add(phi)
        graph["atoms*"].add(dphi)

# Algorithm 2 

def enforcePositiveTraceContraints(graph):
    phiCount = 1
    for pos in graph["pterms"]:
        while not traceConstraint(graph["$"], pos, graph):
            zeta = getTrace(pos, graph).difference(getTrace(graph["$"], graph)).pop()
            Gamma = set()
            for c in getConstrainedLower(pos, graph["constants"]):
                if zeta not in getLower(c.dual):
                    Gamma.add(c)
            if not Gamma:
                linkNodes(zeta, graph["$"])
            else:
                c = Gamma.pop()
                phi = Node(name="phi"+str(phiCount))
                dphi = Node(name="[phi"+str(phiCount)+"]")
                phiCount += 1
                linkDuals(phi, dphi)
                fullLink(phi, c)
                graph["atoms"].add(phi)
                graph["atoms*"].add(dphi)
                
    
    return
    
# Algorithm 3

def sparseCrossing(a, b, graph):
    psiCount = 1
    A = getConstrainedLower(a, graph["atoms"]).difference(getLower(b))
    U = set()
    for phi in A:
        U = set()
        temp = getLower(b)
        B = temp.intersection(graph["atoms"])
        Delta = graph["dualAtoms"].difference(getLower(phi.dual))
        flag = True
        while Delta or flag:
            epsilon = B.pop()
            DeltaP = Delta.intersection(getLower(epsilon.dual))
            if not Delta or any([x not in DeltaP for x in Delta]) or any([x not in Delta for x in DeltaP]):
                psi = Node(name="psi" + str(psiCount))
                dpsi = Node(name="[psi"+ str(psiCount)+"]")
                linkDuals(psi, dpsi)
                fullLink(psi, phi)
                fullLink(psi, epsilon)
                graph["atoms"].add(psi)
                graph["atoms*"].add(dpsi)
                Delta = DeltaP
                U.add(epsilon)
                psiCount += 1
            
            flag = False
    ecount = 1
    for epsilon in U:
        epsilonp = Node("espilon'"+str(ecount))
        depsilonp = Node("[espilon'"+str(ecount)+"]")
        linkDuals(epsilonp, depsilonp)
        fullLink(epsilonp, epsilon)
        graph["atoms"].add(epsilonp)
        graph["atoms*"].add(depsilonp)
        ecount += 1

    for node in A.union(U):
        graph = deleteNode(node, graph, "atoms", "atoms*")
    
    return graph

def cross(graph):
    for pos in graph["pterms"]:
        sparseCrossing(graph["$"], pos, graph)
    
    return

# Algorithm 4

def reduceAtoms(graph):
    Q = set()
    Lambda = graph["constants"]
    while Lambda:
        c = Lambda.pop()
        Sc = Q.intersection(getLower(c))
        Wc = graph["dualAtoms"]
        if Sc:
            for phi in Sc:
                Wc = Wc.intersection(getConstrainedLower(phi, graph["dualAtoms"]))
        
        Phic = set([x.dual for x in getConstrainedLower(c, graph["atoms"])])
        T = getTrace(c, graph)
        while (not T.issubset(Wc) or not Wc.issubset(T)) and Wc:
            eta = Wc.difference(T).pop()
            temp = getUpper(eta).union(set([graph["Base"].dual]))
            t2 = Phic.difference(temp)
            phi = t2.pop().dual
            #Phic.remove(phi.dual)
            Q.add(phi)
            Wc = Wc.intersection(getConstrainedLower(phi.dual, graph["dualAtoms"]))
    
    for atom in graph["atoms"].difference(Q.union([graph["Base"]])):
        deleteNode(atom, graph, "atoms", "atoms*")

def readData(file):
    with open(file, 'r') as data:
        labels = list(map(lambda x: x.strip(), data.readline().split(",")))
        records = []

        for line in data:
            records.append(list(map(lambda x: x.strip(),line.split(","))))
    
    return labels, records

def classify(dataList, labels, graph):
    consts = {}
    for const in graph["constants"]:
        consts[const.name] = const.children
    
    atoms = set()
    for i in range(1, len(dataList)):
        attr = labels[i] + "_" + dataList[i]
        if attr in consts:
            atoms = atoms.union(consts[attr])
    
    return graph["$"].children.issubset(atoms)

def createGraph(labels, records):
    consts = {}
    layers = {"pterms":set(), "nterms":set(), "constants":set(), 
    "atoms":set(), "terms*":set(), "constants*":set(), 
    "atoms*":set(), "dualAtoms":set(), "$":None, "Base":None}

    consts["$"] = Node(name="$")
    cdual = Node(name="[$]")
    linkDuals(cdual, consts["$"])
    zero = Node(name="0")
    zdual = Node(name="[0]")
    linkDuals(zero, zdual)
    dualBase = Node(name="0*")

    for record in records:
        term = Node(name=record[0])
        dualTerm = Node(name="["+record[0]+"]")
        linkDuals(term, dualTerm)
        for i in range(1, len(record)-1):
            name = labels[i] + "_" + record[i]
            if name in consts:
                fullLink(consts[name], term)
            else:
                consts[name] = Node(name=name)
                dualConst = Node(name="["+name+"]")
                linkDuals(consts[name], dualConst)
                fullLink(consts[name], term)
        
        if record[~0] == "1":
            linkNodes(term.dual, consts["$"].dual)
            layers["pterms"].add(term)
            
        else:
            layers["nterms"].add(term)
        linkNodes(dualBase, term.dual)
        layers["terms*"].add(dualTerm)
    
    for const in consts:
        fullLink(zero, consts[const])
        layers["constants"].add(consts[const])
        layers["constants*"].add(consts[const].dual)
        

    layers["dualAtoms"].add(dualBase)
    layers["atoms"].add(zero)
    layers["atoms*"].add(zdual)
    layers["$"] = consts["$"]
    layers["Base"] = zero
    return layers

labels, records = readData("people.data")
layers = createGraph(labels, records)
enforceNegativeTraceConstraints(layers)
enforcePositiveTraceContraints(layers)
cross(layers)
#reduceAtoms(layers)

with open("test.data", 'r') as file:
    for line in file:
        data = list(map(lambda x: x.strip(), line.split(",")))
        print(data[0], ":", classify(data, labels, layers))

#print(classify(["Steve","blue","brown"], labels, layers))          
    