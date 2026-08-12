==============================
Generate arcs from centrosomes
==============================

1. Tune the settings
====================

The arcs will be generated centered on the centrioles and oriented using the direction
vector from one centriole to the other. The angle of the arcs is evenly splitted on both 
sides of this direction vector. This way, the arc always points from one centriole 
to the other, and the kymograph will always be built in the same direction.

On the following figures:

- R = Radius of the arc
- A = 2×α = Angle of the arc

+------------------------------------+-----------------------------------+------------------------------------+
| .. image:: _images/vectors.png     | .. image:: _images/arc.png        | .. image:: _images/distance.png    |
|   :align: center                   |   :align: center                  |   :align: center                   |
+------------------------------------+-----------------------------------+------------------------------------+

+-----------------------+-----------------------------------------------------------------+
| Name                  | Description                                                     |
+=======================+=================================================================+
| Radius                | The radius of the arc centered on the centriole.                |
+-----------------------+-----------------------------------------------------------------+
| Angle                 | The angle of the arc centered on the centriole.                 |
+-----------------------+-----------------------------------------------------------------+

2. Launch the process
=====================

If you click on "Build arcs", the process should be much quicker than the detection and tracking of centrioles. 
There will be one new layer per arc (so per centriole) in the Napari viewer.
There layers are named "_Arc XX" with XX being the ID of the centriole.

.. figure:: _images/arcs_napari.png
  :align: center
  :width: 60%

  The two layers associated to the two arcs.