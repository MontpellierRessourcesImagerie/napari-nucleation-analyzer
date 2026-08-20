===============================
Sample arcs to build kymographs
===============================

A kymograph consists in taking a line of pixels on several time points and stack them to
form a new 2D image. If the line of pixels is P pixels long and the image has T time points, 
the resulting kymograph should be a 2D image of size P×T.

In our case, the line that we use is the arc that we build from the centrosomes.

If you click on the "Build kymographs" button, the kymographs will be built for each arc. 
They should show up on the right side of your original image in the viewer.
Each kymograph should have a colored outline indicating from which centrosome it comes from.
Above each kymograph, you should see something like "PXX -> cYY" with XX being the ID of the 
pair and YY being the ID of the centrosome.

**Note 01:** From this step on, you will certainly have the message: :code:`WARNING: Inconsistent units across layers; units will not be used for rendering.`. 
It is a normal behavior, it is due to the fact that kymographs cannot be calibrated. Indeed, the Y axis of the 
kymograph is actually time while the X axis is a distance. So it is impossible to have an actual physical size 
as for all other layers.

.. figure:: _images/kymos.png
  :align: center
  :width: 60%

  The two orange kymos come from the same pair and so do the the two cyan ones.
  The two first ones come from the pair of ID 1 and the two last ones come from the 
  centrosome of ID 2. The source centrosomes are in order: 4, 5, 1, 6.

**Note 02:** The IDs of centrosomes are not necessarily in order nor consecutive. Indeed, the IDs 
are created when centrosomes are detected and only the centrosomes bound to your hints are kept. 
So if there are objects detected as centrosomes that are not bound to any hint, they will be discarded 
but they keep their ID. You only see the IDs of the centrosomes that are bound to your hints.