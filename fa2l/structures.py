"""
Data structures
"""

class Node:
    def __init__(self):
        self.mass = 0
        self.old_dx = 0
        self.old_dy = 0
        self.dx = 0
        self.dy = 0
        self.x = 0
        self.y = 0
        self.size = 10
        self.weight = 1


class Edge:
    def __init__(self):
        self.node1 = -1
        self.node2 = -1
        self.weight = 0
