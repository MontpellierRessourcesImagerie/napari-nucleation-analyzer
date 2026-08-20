=================
Detailed workflow
=================


1. Detection and tracking of centrosomes
========================================

a. Detect centrosomes
---------------------

To make the detection of the centrosomes an easier task, several filters are applied time 
point per time point to the whole image.

1. A median filter with a square kernel of 3×3 pixels is applied to reduce noise.
2. The background is estimated and subtracted using morphological opening bigger than a centrosome.
3. A Laplacian of Gaussian filter is applied to enhance the centrosomes.
4. The result is normalized between 0.0 and 1.0 and inverted.

.. figure:: _images/prefilter.png
  :align: center
  
  a) original image, b) after median filter, c) isolated background, d) background subtracted 
  from median filtered image, e) after Laplacian of Gaussian filter, f) normalized and inverted image.

- On this filtered image, we look for local maxima that are above a user-defined prominence, the 
  prominence being the difference between the local maximum and its most shalow local minimum.

.. math::
  prominence_{final} = prominence_{user} \cdot StdDev_{filtered}

.. figure:: _images/prominence.png
  :align: center
  :width: 50%
  
  P1 is the prominence of the first local maximum, P2 is the prominence of the second local maximum.

This step produces a collection of 2D points for each frame that are candidates for being centrosome. 
There can still be false positives or centrosomes you are not interested in.

b. Track the centrosomes
------------------------

The tracking of centrosomes relies on [TrackPy]_, a Python library for particle tracking.
The `nearest velocity algorithm`_ is used to link the detected centrosomes over time.
The "memory" parameter determines for how many frames a centrosome can disappear without changing 
of identity. The "searching range" parameter determines how far a centrosome can move between two consecutive
time points before it is considered as a new centrosome. 

The higher these two parameters, the better centrosomes will be tracked but the more likely it is to make mistakes.

At the end of this step, we have a collection of trajectories (== a 2D coordinate for each time point) 
that are candidates for being centrosomes.


c. From detected centrosomes to centrosomes
-------------------------------------------

The only thing allowing to know which centrosomes are bound together to form a pair  
and which one are worth analyzing are the hints that you provided.

For each hint:

- We take the centrosome candidates at the designated starting point as well as the hints for 
  this same time point.
- For each hint point, we look at the closest candidate centrosome. If it is closer than the 
  "max binding distance" parameter, we link this centrosome to the hint and we consider it as a 
  centrosome of interest. By finding the candidate at this time point, we find it for the whole 
  trajectory thanks to the previous tracking step.
- To have a dense and uniform sampling, missing points on trajectories are interpolated using a 
  linear interpolation. This is done for each centrosome of interest.
- The trajectories are filtered to keep only the time points from the requested time range.

At this moment, we have a trajectory for each centrosome of interest. Each centrosome has a unique ID 
and each pair has a unique ID. Thanks to the hints, each centrosome ID is associated to a pair ID. 
The two centrosomes of a same pair are linked together by a line in the Napari viewer.

The result is stored in a :code:`pandas.DataFrame` (CSV-like structure) having one line per point per time. 
The columns present are:

- **T**: the time point index
- **pair_id**: the pair ID
- **centrosome_id**: the centrosome ID
- **X**: the X coordinate of the centrosome in pixels
- **Y**: the Y coordinate of the centrosome in pixels
- **track_id**: a duplication of the centrosome ID, used for compatibility with Napari track layers.

2. Build arcs from centrosomes
==============================

a. Process vectors
------------------

At this point, we have for each pair, the two points corresponding to the centrosomes. 
To know towards which direction the arcs should be built, we need to compute a direction vector
for each of them. If we say that C1 is the position of a centrosome and C2 the position of the other, 
the vectors are equal to:

.. math::

   \overrightarrow{v_{1}} &= \left| C_{2} - C_{1} \right| \\
   \overrightarrow{v_{2}} &= -\overrightarrow{v_{1}}


b. Build the arcs
-----------------

The number of points on the arc is processed as being a fraction of the perimeter of the whole circle. 
In the following equation, the arc angle is noted 𝜃 and r is its radius.

.. math::

   L_{\text{arc}} &= \pi \cdot r_{\text{px}} \cdot \frac{\theta_{deg}}{360} \\
   n_{\text{points}} &= \max\left(2, \left\lceil L_{\text{arc}} \right\rceil\right)

Once we have the number of points, we can generate a "neutral" arc that is centered on the 
origin and oriented along the X axis (like the X axis is the direction vector). 
The points are generated using a parametric equation of a circle:

.. math::

   x &= r_{\text{px}} \cdot \cos(\alpha) \\
   y &= r_{\text{px}} \cdot \sin(\alpha) \\
   \alpha &\in \left[-\frac{\theta_{rad}}{2}, \frac{\theta_{rad}}{2}\right]

Since the arcs are oriented along the X axis, the next step consists in processing the angle 
between the X-axis' vector and the direction vector of the centrosomes. We name this angle 𝛽.

.. math::

   \beta = \operatorname{atan2}\left(v_{y}, v_{x}\right)

This angle is used to build a rotation matrix, which, applied to the neutral arc, will rotate 
it around the origin to be oriented along the direction vector of the centrosomes.

.. math::

   R = \begin{bmatrix} \cos\beta & \sin\beta \\ -\sin\beta & \cos\beta \end{bmatrix}

The last step consists in translating the arc to be centered on the centrosomes. This is a simple
addition of the coordinates of the centrosomes to the coordinates of the points on the arc.

We just produced for each centrosome of interest and for each point of time, a collection of points 
that are the coordinates of the arc.


3. Sample arcs to build kymographs
==================================

We have arcs center on centrosomes and oriented along the direction vector of the centrosomes.
We now need to sample the pixel values along these arcs and stack them to build kymographs. 
However, points coordinates are in sub-pixel precision, which mean that each coordinate doesn't
correspond exactly to a pixel, but rather a combination of 4 pixels.

.. figure:: _images/sampling.png
  :align: center
  :width: 30%
  
  The points of the arc are not exactly on the pixels, we need to make a linear combination of the 4 neighbors.

The coefficients of the linear interpolation are computed using the distance between the point and the 4 neighbors.


4. Searching for spots on kymographs
====================================

The same prominence-based detection algorithm that was used to detect centrosomes is used to detect spots on the kymographs.


---

.. _nearest velocity algorithm: https://soft-matter.github.io/trackpy/dev/tutorial/prediction.html
.. [TrackPy] Allan, D. B., Caswell, T., Keim, N. C., van der Wel, C. M.& Verweij, R. W. (2025). soft-matter/trackpy: v0.7 (Version v0.7) [Computer software]. Zenodo. https://doi.org/10.5281/zenodo.16089574