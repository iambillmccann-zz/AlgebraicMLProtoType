class Node(object):
    """
    The nodes are the verticies of the graph.
    """

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
