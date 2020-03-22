from AlgebraicML.nodes import Node
from AlgebraicML.graph import Graph

import sys

# Algorithm 1

def enforceNegativeTraceConstraints(graph):
    """
    """
    # Preprocessing Step
    count = 1
    for neg in graph.layers["nterms"]:
        zeta = Node(name="zeta"+str(count))
        graph.linkNodes(zeta, neg.dual)
        graph.layers["dualAtoms"].add(zeta)
        count += 1

    # Main Algorithm
    for neg in graph.layers["nterms"]:
        if not graph.traceConstraint(graph.layers["$"], neg):
            continue
        
        c = None
        while not c:
            c = graph.findStronglyDiscriminantCoefficient(graph.layers["$"], neg)
            if not c:
                h = graph.getConstrainedLower(neg.dual, graph.layers["terms*"]).difference(graph.getLower(graph.layers["$"].dual)).pop()
                zeta = Node(name="zeta"+str(count))
                graph.linkNodes(zeta, h.dual)
                graph.layers["dualAtoms"].add(zeta)
                count += 1
        
        phi = Node(name="phi")
        dphi = Node(name="[phi]")
        graph.linkDuals(phi, dphi)
        graph.fullLink(phi, c)
        graph.layers["atoms"].add(phi)
        graph.layers["atoms*"].add(dphi)

# Algorithm 2 

def enforcePositiveTraceContraints(graph):
    """
    """
    phiCount = 1
    for pos in graph.layers["pterms"]:
        while not graph.traceConstraint(graph.layers["$"], pos):
            zeta = graph.getTrace(pos).difference(graph.getTrace(graph.layers["$"])).pop()
            Gamma = set()
            for c in graph.getConstrainedLower(pos, graph.layers["constants"]):
                if zeta not in graph.getLower(c.dual):
                    Gamma.add(c)
            if not Gamma:
                graph.linkNodes(zeta, graph.layers["$"])
            else:
                c = Gamma.pop()
                phi = Node(name="phi"+str(phiCount))
                dphi = Node(name="[phi"+str(phiCount)+"]")
                phiCount += 1
                graph.linkDuals(phi, dphi)
                graph.fullLink(phi, c)
                graph.layers["atoms"].add(phi)
                graph.layers["atoms*"].add(dphi)
                
    
    return
    
# Algorithm 3

def sparseCrossing(a, b, graph, count):
    """
    """
    psiCount = 1
    A = graph.getConstrainedLower(a, graph.layers["atoms"]).difference(graph.getLower(b))
    U = set()
    for phi in A:
        U = set()
        B = graph.getConstrainedLower(b, graph.layers["atoms"])
        Delta = graph.layers["dualAtoms"].difference(graph.getLower(phi.dual))
        flag = True
        while Delta or flag:
            epsilon = B.pop()
            DeltaP = Delta.intersection(graph.getLower(epsilon.dual))
            if not Delta or (not Delta.issubset(DeltaP) or not DeltaP.issubset(Delta)):
                psi = Node(name="psi" + str(psiCount))
                dpsi = Node(name="[psi"+ str(psiCount)+"]")
                graph.linkDuals(psi, dpsi)
                graph.fullLink(psi, phi)
                graph.fullLink(psi, epsilon)
                graph.layers["atoms"].add(psi)
                graph.layers["atoms*"].add(dpsi)
                Delta = DeltaP
                U.add(epsilon)
                psiCount += 1
            
            flag = False

    ecount = 1
    for epsilon in U:
        epsilonp = Node("epsilon'"+str(ecount))
        depsilonp = Node("[epsilon'"+str(ecount)+"]")
        graph.linkDuals(epsilonp, depsilonp)
        graph.fullLink(epsilonp, epsilon)
        graph.layers["atoms"].add(epsilonp)
        graph.layers["atoms*"].add(depsilonp)
        ecount += 1

    #for node in A.union(U):
    #    deleteNode(atom, graph, "atoms", "atoms*")
    Joe = set()
    for atom in graph.layers["atoms"]:
        if atom.children:
            Joe.add(atom)
    
    for atom in Joe:
        graph.deleteNode(atom, "atoms", "atoms*")

    return graph

def cross(graph):
    """
    """
    count = 1
    for pos in graph.layers["pterms"]:
        sparseCrossing(graph.layers["$"], pos, graph, count)
        count += 1
    
    return

# Algorithm 4

def reduceAtoms(graph):
    """
    """
    Q = set()
    Lambda = graph.layers["constants"]
    while Lambda:
        c = Lambda.pop()
        Sc = Q.intersection(graph.getLower(c))
        Wc = graph.layers["dualAtoms"]
        if Sc:
            for phi in Sc:
                Wc = Wc.intersection(graph.getConstrainedLower(phi, graph.layers["dualAtoms"]))
        
        Phic = set([x.dual for x in graph.getConstrainedLower(c, graph.layers["atoms"])])
        T = graph.getTrace(c)
        count = 0
        while (not T.issubset(Wc) or not Wc.issubset(T)) and Wc:
            eta = Wc.difference(T).pop()
            temp = Phic.difference(graph.getUpper(eta))
            phi = list(temp)[count%len(temp)].dual
            count += 1
            Q.add(phi)
            Wc = Wc.intersection(graph.getConstrainedLower(phi.dual, graph.layers["dualAtoms"]))
    
    for atom in graph.layers["atoms"].difference(Q.union([graph.layers["Base"]])):
        graph.deleteNode(atom, "atoms", "atoms*")

# Observe Atoms

def obeserveAtoms(graph):
    """
    """
    importantAtoms = graph.layers["$"].children
    importantAtoms.discard(graph.layers["Base"])
    for atom in importantAtoms:
        pset = atom.parents
        print(" or ".join(map(lambda x: x.name, pset.difference(set([graph.layers["$"]])))))

def verifyTraceConstraints(graph):
    """
    """
    for neg in graph.layers["nterms"]:
        if graph.traceConstraint(graph.layers["$"], neg):
            return False
    for pos in graph.layers["pterms"]:
        if not graph.traceConstraint(graph.layers["$"], pos):
            return False
    return True


def classify(dataList, labels, graph):
    consts = {}
    for const in graph.layers["constants"]:
        consts[const.name] = const.children
    
    atoms = set()
    for i in range(1, len(dataList)):
        attr = labels[i] + "_" + dataList[i]
        if attr in consts:
            atoms = atoms.union(consts[attr])
    
    return graph.layers["$"].children.issubset(atoms)
