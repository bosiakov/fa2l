import math

"""
The crucial idea in speeding up the brute force n-body algorithm is to group nearby bodies and approximate them as a 
single body. If the group is sufficiently far away, we can approximate its gravitational effects by using its 
center of mass. The center of mass of a group of bodies is the average position of a body in that group, 
weighted by mass.
"""


class Quadtree:
    """
    A quad-tree is similar to a binary tree, except that each node has 4 children (some of which may be empty).
    Each node represents a region of the two dimensional space.
    The topmost node represents the whole space, and its four children represent the four quadrants of the space.
    The space is recursively subdivided into quadrants until each subdivision contains 0 or 1 bodies
    (some regions do not have bodies in all of their quadrants)
    """

    nodes = []

    sum_mass = 0.0
    center_x = 0.0
    center_y = 0.0
    size = 0

    # Cardinal direction

    NW = None  # left top
    NE = None  # right top
    SW = None  # left bottom
    SE = None  # right bottom

    def __init__(self, nodes):
        self.nodes = nodes
        self.compute()

    def compute(self):
        if len(self.nodes):
            self.sum_mass = 0
            sum_x = 0
            sum_y = 0

            for node in self.nodes:
                # print(node.mass)
                self.sum_mass += node.mass
                sum_x += node.x * node.mass
                sum_y += node.y * node.mass

            self.center_x = sum_x / self.sum_mass
            self.center_y = sum_y / self.sum_mass

            size = float('-inf')
            for n in self.nodes:
                distance = math.sqrt(
                    (n.x - self.center_x) * (n.x - self.center_x) + (n.y - self.center_y) * (
                        n.y - self.center_y))
                size = max(size, 2 * distance)

            self.size = size

    def build(self):
        if len(self.nodes) > 1:
            NW = []
            NE = []
            SW = []
            SE = []

            for node in self.nodes:
                if node.y > self.center_y:
                    if node.x > self.center_x:
                        NE.append(node)
                    if node.x < self.center_x:
                        NW.append(node)
                else:
                    if node.x > self.center_x:
                        SW.append(node)
                    if node.x < self.center_x:
                        SE.append(node)

            self.NW = Quadtree(NW)
            self.NE = Quadtree(NE)
            self.SW = Quadtree(SW)
            self.SE = Quadtree(SE)

            for reg in (self.NW, self.NE, self.SW, self.SE):
                reg.build()

    def apply_force(self, node, force, theta):
        if len(self.nodes) > 1:
            distance = math.sqrt(
                (node.x - self.center_x) * (node.x - self.center_x) + (node.y - self.center_y) * (
                    node.y - self.center_y))
            if distance * theta > self.size:
                force.apply_approximation(node, self)
            else:
                for reg in (self.NW, self.NE, self.SW, self.SE):
                    reg.apply_force(node, force, theta)
        if len(self.nodes) == 1:
            force.apply_node_to_node(node, self.nodes[0])
