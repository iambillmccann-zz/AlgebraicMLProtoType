from AlgebraicML.nodes import Node

class Graph(object):
    """
    Hold and manage the Graph
    NOTE. The Graph is a data structure of vertices and edges (nodes and links).
    """
    def __init__(self, labels, records):
        consts = {}
        self.layers = { "pterms" : set(), 
            "nterms" : set(), 
            "constants" : set(), 
            "atoms" : set(), 
            "terms*" : set(), 
            "constants*" : set(), 
            "atoms*" : set(), 
            "dualAtoms" : set(), 
            "$" : None, 
            "Base" : None }

        consts["$"] = Node(name="$")
        cdual = Node(name="[$]")
        self.linkDuals(cdual, consts["$"])
        zero = Node(name="0")
        zdual = Node(name="[0]")
        self.linkDuals(zero, zdual)
        dualBase = Node(name="0*")

        for record in records:

            # for each record in the data, create a vertice (node) and set up initial links (edges)
            term = Node(name=record[0])
            dualTerm = Node(name="["+record[0]+"]")
            self.linkDuals(term, dualTerm)
            for i in range(1, len(record)-1):
                name = labels[i] + "_" + record[i]
                if name in consts:
                    self.fullLink(consts[name], term)
                else:
                    consts[name] = Node(name=name)
                    dualConst = Node(name="["+name+"]")
                    self.linkDuals(consts[name], dualConst)
                    self.fullLink(consts[name], term)

            # Set up the prediction variable
            if record[~0] == "1":
                self.linkNodes(term.dual, consts["$"].dual)
                self.layers["pterms"].add(term)
            else:
                self.layers["nterms"].add(term)

            self.linkNodes(dualBase, term.dual)
            self.layers["terms*"].add(dualTerm)

        # Set up the Graph constants
        for const in consts:
            self.fullLink(zero, consts[const])
            self.layers["constants"].add(consts[const])
            self.layers["constants*"].add(consts[const].dual)

        # Set the layers
        self.layers["dualAtoms"].add(dualBase)
        self.layers["atoms"].add(zero)
        self.layers["atoms*"].add(zdual)
        self.layers["$"] = consts["$"]
        self.layers["Base"] = zero

    def linkNodes(self, child, parent):
        child.addParent(parent)
        parent.addChild(child)

    def linkDuals(self, node, dual):
        node.dual = dual
        dual.dual = node

    def fullLink(self, child, parent):
        self.linkNodes(child, parent)
        self.linkNodes(parent.dual, child.dual)

    def deleteNode(self, node, layer, dlayer):
        """ Delete a node from the Graph
        Loop through the parent vertices and remove those edges, then
        loop through the child vertices and remove those edges too.
        Discard the node, and if necessary, discard the dual.
        Lastly establish edges between all premutations of parents and children.

        Args:
            node        The node to remove
            layer       The layer to remove the node from
            dlayer      The dual layer to remove the node from

        Returns:
            N/A
        """
        for parent in node.parents:
            parent.removeChild(node)
            parent.dual.removeParent(node)
        for child in node.children:
            child.removeParent(node)
            child.dual.removeChild(node)

        self.layers[layer].discard(node)
        if dlayer:
            self.layers[dlayer].discard(node.dual)

        for child in node.children:
            for parent in node.parents:
                self.fullLink(child, parent)

    def getLower(self, node):
        """ Make a new Graph where node is the root, and descendents comprise the remainder
            of the Graph.

        Args:
            node           The node to make the root of the new graph

        Returns:
            A set representing a Graph where node is the root
        """
        lowerset = set([node])

        for child in node.children:
            lowerset = lowerset.union(self.getLower(child))
    
        return lowerset

    def getConstrainedLower(self, node, layer):
        return self.getLower(node).intersection(layer)

    def getUpper(self, node):
        """ Make a new Graph where node is the root, and parents comprise the remainder
            of the Graph.

        Args:
            node           The node to make the root of the new graph

        Returns:
            A set representing a Graph where node is the root
        """
        upperset = set()

        for parent in node.parents:
            upperset = upperset.union(self.getLower(parent))
    
        return upperset

    def getConstrainedUpper(self, node, layer):
        return self.getUpper(node).intersection(layer)

    def getTrace(self, node):
        """

        Args:
            node

        Returns:
            a new set representing the trace
        """
        atoms = self.getConstrainedLower(node, self.layers["atoms"])
        trace = self.layers["dualAtoms"]
        for atom in atoms:
            trace = trace.intersection(self.getLower(atom.dual))
        return trace

    def traceConstraint(self, a, b):
        return self.getTrace(b).issubset(self.getTrace(a))

    def findStronglyDiscriminantCoefficient(self, a, b):
        """ Find the discriminant coefficient

        Args:
            a       Graph a
            b       Graph b

        Returns:
            Either the dual or none
        """
        Omega = set(map(lambda c: c.dual, self.getConstrainedLower(a, self.layers["constants"])))
        U = self.getTrace(b)
        while U:
            zeta = U.pop()
            T = Omega.difference(self.getUpper(zeta))
            if T:
                c = T.pop()
                return c.dual
        return None
