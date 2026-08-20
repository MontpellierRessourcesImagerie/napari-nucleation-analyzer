===================
Nucleation Analyzer
===================

.. figure:: _images/screen-all.png
  :align: center
  :width: 80%

1. Introduction
===============

"Nucleation Analyzer" (NA) from :code:`napari-nucleation-analyzer` is a plugin 
for `Napari <https://napari.org/>`_ that allows to:

- Detect and track centrosomes over time on 2D+t images.
- Bind them in pairs.
- Build arcs at a user-defined distance and angle from each centrosome.
- Build kymographs along the arcs.
- Detect spots on the kymographs.
- Report the number of events per point of time per kymograph.

The images expected by this plugin are 2D+t fluorescence images.

.. figure:: _images/demo-img.png
  :align: center
  
  Starting from the original image, we bind and track a pair of centrosomes, 
  then we build arcs and kymographs to detect events (spots) over time.

2. Principle
============

We work on 2D+t fluorescence images. These images can't be registered because centrosomes
naturally move away from each other over time. Doing so would artificially change the distance between 
centrosomes and would bias the analysis. 
Instead we detect centrosomes on each frames, track them over time, and bind them in pairs. 
Then we build arcs at a user-defined distance and angle from each centrosome, and build kymographs along the arcs.
Finally we detect spots on the kymographs to report the number of events per point of time per kymograph.

.. toctree::
   :maxdepth: 2
   :caption: Content
   
   install_guide
   user_guide
   workflow