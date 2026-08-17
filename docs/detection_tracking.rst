==================================================
Detect and track centrioles to rebuild centrosomes
==================================================

1. Setup your workspace
=======================

a. Open an image
----------------

- In Napari, open the image you want to analyze. It should be a 2D+t fluorescence image. 
  You should be able to simply drag and drop it into the viewer. 
- The slider beneath the viewer should allow you to scroll through the time points of your image.

.. figure:: _images/drag_image.png
  :align: center
  :width: 30%

b. Calibrate your image
-----------------------

- In the top-bar menu of Napari, go for the "Set Scale" tool in the "Calibration Tool" plugin.
- It should open a new widget in the right panel of Napari.
- You will be asked the physical size of your pixels on the X, Y and Z axes. You can provide 
  the actual size for the X and Y axes (ex: 0.103) but you should set the Z axis to 1 (since 
  you don't actually have one).
- Now, you have to tell the axes order of this image. You only have two spacial dimensions 
  and time so your axes should be set to "TYX".
- If you click on the "Apply" button, this configuration will be applied to your current image.
- The default pixel size is 1.0 for all axes so after applying your calibration, your
  image will become very small in the viewer and move in the corner. To address this problem, you
  can click on the |reset_view| "Reset view" button in the lower-left part of Napari's window.

2. Provide hints for the centrosomes
====================================

- At this point, you can open the "Nucleation analyzer" plugin in the top-bar menu of Napari. 
  It should open a new widget in the right panel of Napari.
- A hint consists in:

    - a hand-drawn line between two centrioles to indicate that they belong to the same centrosome
    - the time point at which we should start analyzing this centrosome
    - the time point at which we should stop analyzing this centrosome

- Here is the sequence of actions to do to give a hint for a centrosome:

    1. Identify the first centrosome that you want to study and navigate to the first time point at which
       you want to study it. Zoom on it using the mouse wheel and pan using the left mouse button.
    2. In the "Centrosome tracks" sub-section of the widget, click on the "Add centrosome" button. 
       A new row should appear with its own color and two buttons to set the start and stop time points.
    3. In the list of layers, a new layer named "_Centrosome XX" should have appeared and be selected.
       In the upper-left region of Napari's window, make sure the the |move| "Move" button is selected.
       While it is selected, you can still navigate in the viewer.
    4. Activate the |line| "Line" tool still in the upper-left region of Napari's window. 
       You can now draw a line going from one centriole to the other of the same centrosome.
       For safety reasons, you should go back to the |move| "Move" tool after drawing the line 
       to avoid accidentally drawing another one. If you need to adjust the line, you have to
       use the |move_verts| "Move vertices" tool.
    5. Get back in the "Centrosome tracks" sub-section and click on the "Start" button. It will
       set the start time point of this centrosome to the current time point.
    6. Navigate to the last time point at which you want to study this centrosome. The line that you
       drew won't follow, it is normal. Click on the "End" button to set the stop time point of 
       this centrosome to the current time point.

.. tabs::

   .. tab:: Step 1

      +----------------------------------------------------+
      | .. image:: _images/def_centro/01-centro.png        |
      |   :align: center                                   |
      +----------------------------------------------------+

   .. tab:: Step 2

      +----------------------------------------------------+
      | .. image:: _images/def_centro/02-centro.png        |
      |   :align: center                                   |
      +----------------------------------------------------+

   .. tab:: Step 3

      +----------------------------------------------------+
      | .. image:: _images/def_centro/03-centro.png        |
      |   :align: center                                   |
      +----------------------------------------------------+
   
   .. tab:: Step 4

      +----------------------------------------------------+
      | .. image:: _images/def_centro/04-centro.png        |
      |   :align: center                                   |
      +----------------------------------------------------+

   .. tab:: Step 5

      +----------------------------------------------------+
      | .. image:: _images/def_centro/05-centro.png        |
      |   :align: center                                   |
      +----------------------------------------------------+

3. Tune the settings
====================

+-----------------------+--------------------------------------------------------------------+------------+
| Name                  | Description                                                        | Default    |
+=======================+====================================================================+============+
| Image                 | This dropdown menu allows you to select the image you want to      |            |
|                       | analyze. This is simply the list of image layers currently open in |            |
|                       | Napari.                                                            |            |
+-----------------------+--------------------------------------------------------------------+------------+
| Prominence            | This value must be strictly positive. It is a factor that          | ×10.0      |
|                       | represents how much a local maximum on the preprocessed image must |            |
|                       | be above the filtered image's standard deviation to be considered  |            |
|                       | as a centriole. A local maximum on the preprocessed image must be  |            |
|                       | above its most shallow local minimum to be considered as a         |            |
|                       | centriole. The higher this value, the less centrioles will be      |            |
|                       | detected. More detailed explanation are available in               |            |
|                       | :doc:`workflow`                                                    |            |
+-----------------------+--------------------------------------------------------------------+------------+
| Searching range       | In the tracking process, how much a centriole is allowed to move   | 3.25 µm    |
|                       | between two consecutive time points. The higher this value, the    |            |
|                       | more centrioles will be tracked but the more likely it is to make  |            |
|                       | mistakes.                                                          |            |
+-----------------------+--------------------------------------------------------------------+------------+
| Memory                | In the tracking process, how many frames a centriole is allowed to | 10 f       |
|                       | disappear for and still be linked to a trajectory. The higher this |            |
|                       | value, the more robust the tracking will be, but the more likely   |            |
|                       | it is to make mistakes. This is useful when some centrioles move   |            |
|                       | on another Z.                                                      |            |
+-----------------------+--------------------------------------------------------------------+------------+
| Max binding distance  | When the actual detected centrioles have to be linked to your      | 0.6 µm     |
|                       | manually drawn hints, this value is the maximum distance allowed   |            |
|                       | between a detected centriole and a hint. If a detected centriole   |            |
|                       | is further than this distance from a hint, it won't be linked to   |            |
|                       | it.                                                                |            |
+-----------------------+--------------------------------------------------------------------+------------+

4. Launch the calculation
=========================

You can now click on the "Find centrosomes" button to launch the calculation. This is a heavy
process, it can take a minute or two depending on the size of your image. It performs the detection, 
the tracking and the linking of the centrioles in one go.

At the end of the calculation, if you navigate through the time points of your image, you should see 
the detected and tracked centrioles and the link between them to form centrosomes.

.. figure:: _images/tracked.gif
  :align: center
  :width: 60%

  After the operation, the line points from one centriole to the other of the same centrosome 
  at all time points.


.. |reset_view| image:: _images/reset_view.png
   :width: 20px
   :height: 20px
.. |line| image:: _images/line.png
   :width: 20px
   :height: 20px
.. |move| image:: _images/move.png
   :width: 20px
   :height: 20px
.. |move_verts| image:: _images/move_vertices.png
   :width: 20px
   :height: 20px