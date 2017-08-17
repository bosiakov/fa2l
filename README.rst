Force Atlas 2 Layout
===========================

ForceAtlas2 is a continuous graph layout algorithm for handy network visualization.

This implementation is based on this `paper <http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0098679>`_.

**Warning:** Some features (especially *Prevent Overlapping*) are not completely implemented. I'm waiting for your pull-requests.

Installing
----------

Supports Python 3.3+

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


Documentation
-------------

You will find all the documentation in the source code

Copyright
---------

Copyright Eugene Bosiakov. Licensed under the GNU GPLv3.

This files are based on the java files included in Gephi (Copyright 2011 Gephi Consortium).

Also thanks to Max Shinn.