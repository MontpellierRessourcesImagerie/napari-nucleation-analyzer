===============================
Sample arcs to build kymographs
===============================

A kymograph consists in taking a line of pixels on several time points and stack them to
form a new 2D image. If the line of pixels is P pixels long and the image has T time points, 
the resulting kymograph should be a 2D image of size P×T.

In our case, the line that we use is the arc that we build from the centrioles.

If you click on the "Build kymographs" button, the kymographs will be built for each arc. 
All the other layers will be hidden in the viewer, and there will be one new layer per kymograph in the Napari viewer.
Each kymograph should have a colored outline indicating from which centrosome it comes from.
Above each kymograph, you should see something like "CXX -> cYY" with XX being the ID of the 
centrosome and YY being the ID of the centriole.

**Note:** From this step on, you will certainly have the message: :code:`WARNING: Inconsistent units across layers; units will not be used for rendering.`. 
It is a normal behavior, it is due to the fact that kymographs cannot be calibrated. Indeed, the Y axis of the 
kymograph is actually time while the X axis is a distance. So it is impossible to have an actual physical size 
as for all other layers.

.. figure:: _images/kymos.png
  :align: center
  :width: 60%

  The two orange kymos come from the same centrosome and so do the the two cyan ones.
  The two first ones come from the centrosome of ID 1 and the two last ones come from the 
  centrosome of ID 2. The source centrioles are in order: 1, 5, 3 and 4.