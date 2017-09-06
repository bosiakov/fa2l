import numpy
import random
import math
import networkx
from .structures import Node, Edge
from .approximation import Quadtree
from .force import apply_repulsion, apply_gravity, apply_attraction, get_repulsion, get_strong_gravity, get_attraction


def force_atlas2_layout(graph,
                        pos_list=None,
                        node_masses=None,
                        iterations=100,

                        outbound_attraction_distribution=False,
                        lin_log_mode=False,
                        prevent_overlapping=False,
                        edge_weight_influence=1.0,

                        jitter_tolerance=1.0,
                        barnes_hut_optimize=False,
                        barnes_hut_theta=1.2,

                        scaling_ratio=2.0,
                        strong_gravity_mode=False,
                        multithread=False,
                        gravity=1.0):
    """
    Position nodes using ForceAtlas2 force-directed algorithm

    Parameters
    ----------
    graph: NetworkX graph
        A position will be assigned to every node in G.

    pos_list : dict or None  optional (default=None)
        Initial positions for nodes as a dictionary with node as keys
        and values as a coordinate list or tuple.  If None, then use
        random initial positions.

    node_masses : dict or None  optional (default=None)
        Predefined masses for nodes with node as keys and masses as values.
        If None, then use degree of nodes.

    iterations : int  optional (default=50)
        Number of iterations

    outbound_attraction_distribution : boolean
        Distributes attraction along outbound edges. Hubs attract less and thus are pushed to the borders.
        This mode is meant to grant authorities (nodes with a high indegree) a more central position than hubs (nodes with a high outdegree).
        This is useful for social networks and web networks, where authorities are sometimes considered more important than hubs

    lin_log_mode: boolean
        Switch ForceAtlas model from lin-lin to lin-log (tribute to Andreas Noack). Makes clusters more tight

    prevent_overlapping: boolean
        With this mode enabled, the repulsion is modified so that the nodes do not overlap.
        The goal is to produce a more readable and aesthetically pleasing image.

    edge_weight_influence: float
        How much influence you give to the edges weight. 0 is “no influence” and 1 is “normal”.

    jitter_tolerance: float
        How much swinging you allow. Above 1 discouraged. Lower gives less speed and more precision

    barnes_hut_optimize: boolean
        Barnes Hut optimization: n² complexity to n.ln(n) ; allows larger graphs.

    barnes_hut_theta: float
        Theta of the Barnes Hut optimization

    scaling_ratio: float
        How much repulsion you want. More makes a more sparse graph.

    strong_gravity_mode: boolean
        The “Strong gravity” option sets a force that attracts the nodes that are distant from the center more ( is this distance).
        This force has the drawback of being so strong that it is sometimes stronger than the other forces.
        It may result in a biased placement of the nodes.
        However, its advantage is to force a very compact layout, which may be useful for certain purposes.

    multithread: boolean

    gravity: float
        Attracts nodes to the center. Prevents islands from drifting away.

    Returns
    -------
    pos : dict
        A dictionary of positions keyed by node
    """

    assert isinstance(graph, networkx.classes.graph.Graph), "Not a networkx graph"
    assert isinstance(pos_list, dict) or (pos_list is None), "pos must be specified as a dictionary, as in networkx"

    assert multithread is False, "Not implemented yet"

    G = numpy.asarray(networkx.to_numpy_matrix(graph))

    pos = None

    if pos_list is not None:
        pos = numpy.asarray([pos_list[i] for i in graph.nodes()])

    masses = None

    if node_masses is not None:
        masses = numpy.asarray([node_masses[node] for node in graph.nodes()])

    assert G.shape == (G.shape[0], G.shape[0]), "G is not 2D square"
    assert numpy.all(G.T == G), "G is not symmetric."

    # speed and speed efficiency describe a scaling factor of dx and dy
    # before x and y are adjusted.  These are modified as the
    # algorithm runs to help ensure convergence.
    speed = 1
    speed_efficiency = 1

    nodes = []
    for i in range(0, G.shape[0]):
        n = Node()
        if node_masses is None:
            n.mass = 1 + numpy.sum(G[i])
        else:
            n.mass = masses[i]
        n.old_dx = 0
        n.old_dy = 0
        n.dx = 0
        n.dy = 0
        if pos is None:
            n.x = random.random()
            n.y = random.random()
        else:
            n.x = pos[i][0]
            n.y = pos[i][1]
        nodes.append(n)

    edges = []
    es = numpy.asarray(G.nonzero()).T
    for e in es:
        if e[1] <= e[0]: continue  # Avoid duplicate edges
        edge = Edge()
        edge.node1 = e[0]  # The index of the first node in `nodes`
        edge.node2 = e[1]  # The index of the second node in `nodes`
        edge.weight = G[tuple(e)]
        edges.append(edge)

    repulsion = get_repulsion(prevent_overlapping, scaling_ratio)

    if strong_gravity_mode:
        gravity_force = get_strong_gravity(scaling_ratio)
    else:
        gravity_force = repulsion

    if outbound_attraction_distribution:
        outbound_att_compensation = numpy.mean([n.mass for n in nodes])

    attraction_coef = outbound_att_compensation if outbound_attraction_distribution else 1
    attraction = get_attraction(lin_log_mode, outbound_attraction_distribution, prevent_overlapping,
                                attraction_coef)
    # Main loop

    for _i in range(0, iterations):
        for n in nodes:
            n.old_dx = n.dx
            n.old_dy = n.dy
            n.dx = 0
            n.dy = 0

        # Barnes Hut optimization
        root_region = None

        if barnes_hut_optimize:
            root_region = Quadtree(nodes)
            root_region.build()

        apply_repulsion(repulsion, nodes, barnes_hut_optimize=barnes_hut_optimize, barnes_hut_theta=barnes_hut_theta,
                        region=root_region)
        apply_gravity(gravity_force, nodes, gravity, scaling_ratio)

        apply_attraction(attraction, nodes, edges, edge_weight_influence)

        # Auto adjust speed.
        total_swinging = 0.0  # How much irregular movement
        total_effective_traction = 0.0  # How much useful movement
        for n in nodes:
            swinging = math.sqrt((n.old_dx - n.dx) * (n.old_dx - n.dx) + (n.old_dy - n.dy) * (n.old_dy - n.dy))
            total_swinging += n.mass * swinging
            total_effective_traction += .5 * n.mass * math.sqrt(
                (n.old_dx + n.dx) * (n.old_dx + n.dx) + (n.old_dy + n.dy) * (n.old_dy + n.dy))

        # Optimize jitter tolerance.
        # The 'right' jitter tolerance for this network.
        # Bigger networks need more tolerance. Denser networks need less tolerance.
        # Totally empiric.

        estimated_optimal_jitter_tolerance = .05 * math.sqrt(len(nodes))
        min_jt = math.sqrt(estimated_optimal_jitter_tolerance)
        max_jt = 10
        jt = jitter_tolerance * max(min_jt, min(max_jt, estimated_optimal_jitter_tolerance * total_effective_traction /
                                                (len(nodes) ** 2)))

        min_speed_efficiency = 0.05

        # Protective against erratic behavior
        if total_swinging / total_effective_traction > 2.0:
            if speed_efficiency > min_speed_efficiency:
                speed_efficiency *= .5
            jt = max(jt, jitter_tolerance)

        target_speed = jt * speed_efficiency * total_effective_traction / total_swinging

        if total_swinging > jt * total_effective_traction:
            if speed_efficiency > min_speed_efficiency:
                speed_efficiency *= .7
        elif speed < 1000:
            speed_efficiency *= 1.3

        # But the speed shoudn't rise too much too quickly, since it would
        # make the convergence drop dramatically.
        max_rise = .5
        speed = speed + min(target_speed - speed, max_rise * speed)

        # Apply forces.
        if prevent_overlapping:
            for n in nodes:
                swinging = n.mass * math.sqrt(
                    (n.old_dx - n.dx) * (n.old_dx - n.dx) + (n.old_dy - n.dy) * (n.old_dy - n.dy))
                factor = 0.1 * speed / (1 + math.sqrt(speed * swinging))

                df = math.sqrt(math.pow(n.dx, 2) + n.dy ** 2)
                factor = min(factor * df, 10.) / df

                x = n.dx * factor
                y = n.dy * factor
        else:
            for n in nodes:
                swinging = n.mass * math.sqrt(
                    (n.old_dx - n.dx) * (n.old_dx - n.dx) + (n.old_dy - n.dy) * (n.old_dy - n.dy))
                factor = speed / (1.0 + math.sqrt(speed * swinging))
                n.x = n.x + (n.dx * factor)
                n.y = n.y + (n.dy * factor)
    positions = [(n.x, n.y) for n in nodes]
    return dict(zip(graph.nodes(), positions))
