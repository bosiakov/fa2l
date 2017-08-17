Force Atlas 2 Layout
===========================

ForceAtlas2 is a continuous graph layout algorithm for handy network visualization.

This implementation is based on this `paper <http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0098679>`_.

**Warning:** Some features (especially *Prevent Overlapping*) are not completely implemented. I'm waiting for your pull-requests.

Installing
----------

Supports Python 3.3+

Install from pip:

.. code-block:: bash

     pip install fa2l


To build and install run from source:

.. code-block:: bash

     python setup.py install

Usage
-----

.. code-block:: python

    import networkx as nx
    from fa2l import force_atlas2_layout
    import matplotlib.pyplot as plt

    G = nx.erdos_renyi_graph(100, 0.15, directed=False)

    positions = force_atlas2_layout(G,
                                    None,
                                    iterations=1000,

                                    outbound_attraction_distribution=False,
                                    lin_log_mode=False,
                                    prevent_overlapping=False,
                                    edge_weight_influence=1.0,

                                    jitter_tolerance=1.0,
                                    barnes_hut_optimize=True,
                                    barnes_hut_theta=0.5,

                                    scaling_ratio=2.0,
                                    strong_gravity_mode=False,
                                    multithread=False,
                                    gravity=1.0)

    nx.draw_networkx(G, positions, cmap=plt.get_cmap('jet'), node_size=50, with_labels=False)
    plt.show()

Example of social graph rendered with force atlas 2 layout:

.. image:: https://raw.githubusercontent.com/bosiakov/fa2l/master/_static/result.jpg

Features
--------

Force Atlas 2 features these settings:

- *Approximate Repulsion*: Barnes Hut optimization: n² complexity to n.ln(n).
- *Gravity*: Attracts nodes to the center. Prevents islands from drifting away.
- *Dissuade Hubs*: Distributes attraction along outbound edges. Hubs attract less and thus are pushed to the borders.
- *LinLog mode*: Switch ForceAtlas model from lin-lin to lin-log. Makes clusters more tight.
- *Prevent Overlap*. WARNING! Does not work very well.
- *Tolerance*: How much swinging you allow. Above 1 discouraged. Lower gives less speed and more precision.
- *Edge Weight Influence*: How much influence you give to the edges weight. 0 is "no influence" and 1 is "normal".

Documentation
-------------

You will find all the documentation in the source code

Copyright
---------

Copyright Eugene Bosiakov. Licensed under the GNU GPLv3.

This files are based on the java files included in Gephi (Copyright 2011 Gephi Consortium).

Also thanks to Max Shinn.