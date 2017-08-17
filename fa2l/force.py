import math


def apply_repulsion(repulsion, nodes, barnes_hut_optimize=False, region=None, barnes_hut_theta=1.2):
    """
    Iterate through the nodes or edges and apply the forces directly to the node objects.
    """
    if not barnes_hut_optimize:
        for i in range(0, len(nodes)):
            for j in range(0, i):
                repulsion.apply_node_to_node(nodes[i], nodes[j])
    else:
        for i in range(0, len(nodes)):
            region.apply_force(nodes[i], repulsion, barnes_hut_theta)


def apply_gravity(repulsion, nodes, gravity, scaling_ratio):
    """
    Iterate through the nodes or edges and apply the gravity directly to the node objects.
    """
    for i in range(0, len(nodes)):
        repulsion.apply_gravitation(nodes[i], gravity / scaling_ratio)


def apply_attraction(attraction, nodes, edges, edge_weight_influence):
    # Optimization, since usually edgeWeightInfluence is 0 or 1, and pow is slow
    if edge_weight_influence == 0:
        for edge in edges:
            attraction.apply(nodes[edge.node1], nodes[edge.node2], 1)
    elif edge_weight_influence == 1:
        for edge in edges:
            attraction.apply(nodes[edge.node1], nodes[edge.node2], edge.weight)
    else:
        for edge in edges:
            attraction.apply(nodes[edge.node1], nodes[edge.node2], pow(edge.weight, edge_weight_influence))


def get_repulsion(adjust_by_size, coefficient):
    if adjust_by_size:
        return LinRepulsionAntiCollision(coefficient)
    else:
        return LinearRepulsion(coefficient)


def get_strong_gravity(coefficient):
    return StrongGravity(coefficient)


def get_attraction(log_attraction, distributed_attraction, adjust_by_size, coefficient):
    if adjust_by_size:
        if log_attraction:
            if distributed_attraction:
                return LogAttractionDegreeDistributedAntiCollision(coefficient)
            else:
                return LogAttractionAntiCollision(coefficient)

        else:
            if distributed_attraction:
                return LinAttractionDegreeDistributedAntiCollision(coefficient)
            else:
                return LinAttractionAntiCollision(coefficient)
    else:
        if log_attraction:
            if distributed_attraction:
                return LogAttractionDegreeDistributed(coefficient)
            else:
                return LogAttraction(coefficient)

        else:
            if distributed_attraction:
                return LinAttractionMassDistributed(coefficient)
            else:
                return LinAttraction(coefficient)


class AttractionForce:
    def __init__(self, coefficient):
        self.coefficient = coefficient

    def __str__(self):
        return str(self.__class__)

    def apply(self, node1, node2, edge_weight):
        """
        Model for node-node attraction.
        This will adjust the dx and dy values of `node1` (and optionally `node2`).
        """
        raise NotImplementedError


class RepulsionForce:
    """
    Here are all the formulas for attraction and repulsion.
    """

    def __init__(self, coefficient):
        self.coefficient = coefficient

    def __str__(self):
        return str(self.__class__)

    def apply_node_to_node(self, node1, node2):
        """
        Model for node-node repulsion
        """
        raise NotImplementedError

    @staticmethod
    def apply_approximation(self, node, region):
        """
        Model for Barnes Hut approximation
        """
        raise NotImplementedError

    def apply_gravitation(self, node, gravity):
        """
        Model for gravitation (anti-repulsion)
        """
        raise NotImplementedError


class LinearRepulsion(RepulsionForce):
    """
    Repulsion force: Linear
    """

    def __init__(self, *args, **kwargs):
        RepulsionForce.__init__(self, *args, **kwargs)

    def __str__(self):
        return RepulsionForce.__str__(self)

    def apply_node_to_node(self, node1, node2):
        x_dist = node1.x - node2.x
        y_dist = node1.y - node2.y
        distance2 = x_dist * x_dist + y_dist * y_dist  # Distance squared

        if distance2 > 0:
            factor = self.coefficient * node1.mass * node2.mass / distance2
            node1.dx += x_dist * factor
            node1.dy += y_dist * factor
            node2.dx -= x_dist * factor
            node2.dy -= y_dist * factor

    def apply_approximation(self, node, region):
        x_dist = node.x - region.center_x
        y_dist = node.y - region.center_y
        distance = math.sqrt(x_dist ** 2 + y_dist ** 2)

        if distance > 0:
            # NB: factor = force / distance
            factor = self.coefficient * node.mass * region.sum_mass / distance / distance

            node.dx += x_dist * factor
            node.dy += y_dist * factor

    def apply_gravitation(self, node, gravity):
        x_dist = node.x
        y_dist = node.y
        distance = math.sqrt(x_dist ** 2 + y_dist ** 2)

        if distance > 0:
            factor = self.coefficient * node.mass * gravity / distance
            node.dx -= x_dist * factor
            node.dy -= y_dist * factor


class LinRepulsionAntiCollision(RepulsionForce):
    """
    Repulsion force: Strong Gravity (as a Repulsion Force because it is easier)
    """

    def __init__(self, *args, **kwargs):
        RepulsionForce.__init__(self, *args, **kwargs)

    def __str__(self):
        return RepulsionForce.__str__(self)

    def apply_node_to_node(self, node1, node2):
        x_dist = node1.x - node2.x
        y_dist = node1.y - node2.y

        distance = math.sqrt(x_dist ** 2 + y_dist ** 2) - node1.size - node2.size

        if distance > 0:
            factor = self.coefficient * node1.mass * node2.mass / distance / distance

            node1.dx += x_dist * factor
            node1.dy += y_dist * factor
            node2.dx -= x_dist * factor
            node2.dy -= y_dist * factor

        if distance < 0:
            factor = 100 * self.coefficient * node1.mass * node2.mass
            node1.dx += x_dist * factor
            node1.dy += y_dist * factor

            node2.dx -= x_dist * factor
            node2.dy -= y_dist * factor

    def apply_approximation(self, node, region):
        x_dist = node.x - region.center_x
        y_dist = node.y - region.center_y

        distance = math.sqrt(x_dist ** 2 + y_dist ** 2)
        if distance > 0:
            factor = self.coefficient * node.mass * region.sum_mass / distance / distance

            node.dx += x_dist * factor
            node.dy += y_dist * factor
        if distance < 0:
            factor = -self.coefficient * node.mass * region.sum_mass / distance

            node.dx += x_dist * factor
            node.dy += y_dist * factor

    def apply_gravitation(self, node, gravity):
        x_dist = node.x
        y_dist = node.y
        distance = math.sqrt(x_dist ** 2 + y_dist ** 2)

        if distance > 0:
            factor = self.coefficient * node.mass * gravity / distance
            node.dx -= x_dist * factor
            node.dy -= y_dist * factor


class StrongGravity(RepulsionForce):
    """
    Strong gravity force function
    """

    def __init__(self, *args, **kwargs):
        RepulsionForce.__init__(self, *args, **kwargs)

    def __str__(self):
        return RepulsionForce.__str__(self)

    def apply_node_to_node(self, node1, node2):
        """
        Not Relevant
        """
        pass

    @staticmethod
    def apply_approximation(self, node, region):
        """
        Not Relevant
        """
        pass

    def apply_gravitation(self, node, gravity):
        x_dist = node.x
        y_dist = node.y

        distance = math.sqrt(x_dist * x_dist + y_dist * y_dist)
        if distance > 0:
            # NB: factor = force / distance
            factor = self.coefficient * node.mass * gravity
            node.dx -= x_dist * factor
            node.dy -= y_dist * factor


class LinAttraction(AttractionForce):
    def __init__(self, *args, **kwargs):
        AttractionForce.__init__(self, *args, **kwargs)

    def __str__(self):
        return AttractionForce.__str__(self)

    def apply(self, node1, node2, edge_weight):
        x_dist = node1.x - node2.x
        y_dist = node1.y - node2.y
        factor = -self.coefficient * edge_weight

        node1.dx += x_dist * factor
        node1.dy += y_dist * factor

        node2.dx -= x_dist * factor
        node2.dy -= y_dist * factor


class LinAttractionMassDistributed(AttractionForce):
    """
    Attraction force: Linear, distributed by mass (typically, degree)
    """

    def __init__(self, *args, **kwargs):
        AttractionForce.__init__(self, *args, **kwargs)

    def __str__(self):
        return AttractionForce.__str__(self)

    def apply(self, node1, node2, edge_weight):
        x_dist = node1.x - node2.x
        y_dist = node1.y - node2.y

        factor = -self.coefficient * edge_weight / node1.weight

        node1.dx += x_dist * factor
        node1.dy += y_dist * factor

        node2.dx -= x_dist * factor
        node2.dy -= y_dist * factor


class LinAttractionAntiCollision(AttractionForce):
    """
    Attraction force: Linear, with Anti-Collision
    """

    def __init__(self, *args, **kwargs):
        AttractionForce.__init__(self, *args, **kwargs)

    def __str__(self):
        return AttractionForce.__str__(self)

    def apply(self, node1, node2, edge_weight):
        x_dist = node1.x - node2.x
        y_dist = node1.y - node2.y

        distance = math.sqrt(x_dist ** 2 + y_dist ** 2) - node1.size - node2.size

        if distance > 0:
            factor = -self.coefficient * edge_weight

            node1.dx += x_dist * factor
            node1.dy += y_dist * factor

            node2.dx -= x_dist * factor
            node2.dy -= y_dist * factor


class LinAttractionDegreeDistributedAntiCollision(AttractionForce):
    """
    Attraction force: Linear, distributed by Degree, with Anti-Collision
    """

    def __init__(self, *args, **kwargs):
        AttractionForce.__init__(self, *args, **kwargs)

    def __str__(self):
        return AttractionForce.__str__(self)

    def apply(self, node1, node2, edge_weight):
        x_dist = node1.x - node2.x
        y_dist = node1.y - node2.y

        distance = math.sqrt(x_dist ** 2 + y_dist ** 2) - node1.size - node2.size

        if distance > 0:
            factor = -self.coefficient * edge_weight / node1.mass

            node1.dx += x_dist * factor
            node1.dy += y_dist * factor

            node2.dx -= x_dist * factor
            node2.dy -= y_dist * factor


class LogAttraction(AttractionForce):
    """
    Attraction force: Logarithmic
    """

    def __init__(self, *args, **kwargs):
        AttractionForce.__init__(self, *args, **kwargs)

    def __str__(self):
        return AttractionForce.__str__(self)

    def apply(self, node1, node2, edge_weight):
        x_dist = node1.x - node2.x
        y_dist = node1.y - node2.y

        distance = math.sqrt(x_dist ** 2 + y_dist ** 2)
        if distance > 0:
            factor = -self.coefficient * edge_weight * math.log(1 + distance) / distance

            node1.dx += x_dist * factor
            node1.dy += y_dist * factor

            node2.dx -= x_dist * factor
            node2.dy -= y_dist * factor


class LogAttractionDegreeDistributed(AttractionForce):
    """
    Attraction force: Linear, distributed by Degree
    """

    def __init__(self, *args, **kwargs):
        AttractionForce.__init__(self, *args, **kwargs)

    def __str__(self):
        return AttractionForce.__str__(self)

    def apply(self, node1, node2, edge_weight):
        x_dist = node1.x - node2.x
        y_dist = node1.y - node2.y

        distance = math.sqrt(x_dist ** 2 + y_dist ** 2)

        if distance > 0:
            factor = -self.coefficient * edge_weight * math.log(1 + distance) / distance / node1.mass

            node1.dx += x_dist * factor
            node1.dy += y_dist * factor

            node2.dx -= x_dist * factor
            node2.dy -= y_dist * factor


class LogAttractionAntiCollision(AttractionForce):
    """
    Attraction force: Linear, distributed by Degree
    """

    def __init__(self, *args, **kwargs):
        AttractionForce.__init__(self, *args, **kwargs)

    def __str__(self):
        return AttractionForce.__str__(self)

    def apply(self, node1, node2, edge_weight):
        x_dist = node1.x - node2.x
        y_dist = node1.y - node2.y

        distance = math.sqrt(x_dist ** 2 + y_dist ** 2) - node1.size - node2.size

        if distance > 0:
            factor = -self.coefficient * edge_weight * math.log(1 + distance) / distance

            node1.dx += x_dist * factor
            node1.dy += y_dist * factor

            node2.dx -= x_dist * factor
            node2.dy -= y_dist * factor


class LogAttractionDegreeDistributedAntiCollision(AttractionForce):
    """
    Attraction force: Linear, distributed by Degree, with Anti-Collision
    """

    def __init__(self, *args, **kwargs):
        AttractionForce.__init__(self, *args, **kwargs)

    def __str__(self):
        return AttractionForce.__str__(self)

    def apply(self, node1, node2, edge_weight):
        x_dist = node1.x - node2.x
        y_dist = node1.y - node2.y

        distance = math.sqrt(x_dist ** 2 + y_dist ** 2) - node1.size - node2.size

        if distance > 0:
            factor = -self.coefficient * edge_weight * math.log(1 + distance) / distance / node1.mass
            node1.dx += x_dist * factor
            node1.dy += y_dist * factor

            node2.dx -= x_dist * factor
            node2.dy -= y_dist * factor
