====================================
Detect spots over time on kymographs
====================================

1. Detect spots
===============

The detection of spots on kymographs is done using the same principle as the detection of 
centrioles: we search for local maxima on a prefiltered image. The only difference is the 
way the prominence is computed.

When we were detecting centrioles, you were providing the prominence value and it was used 
as such to determine local maxima.

Here, the **spot prominence** parameter is actually a factor. The base prominence is individually 
computed for each kymograph and is equal to the standard deviation of the pixel values on the kymograph. 
The final prominence is this base prominence multiplied by the **spot prominence** factor. 
The base value is 1.0 (since it is a factor) and the higher the factor, the less spots will be detected.

Spots will be displayed as gray circles on each kymograph.

.. figure:: _images/detected_spots.png
  :align: center
  :width: 60%

  The gray circles are the detected spots on the kymographs.

2. Export the results
=====================

a. Export summary
-----------------

The first thing that you can export is a summary of the number of events per point of time per kymograph.
The produced file is a CSV file. All the columns start with either **T** or **Count**.

- **T** is the index of the time point at which the event occured.
- **Count** is the number of events that occured at this time point on this kymograph.

The second part of each column name is an expression of the form **[CXX -> cYY]** with XX being 
the ID of the centrosome and YY being the ID of the centriole.

.. figure:: _images/spreadsheet.png
  :align: center
  :width: 60%

  This is a screenshot from LibreOffice Calc in which the CSV was opened.

b. Export archive
-----------------

In order to keep track of all the results, you can export an archive containing all the kymographs
by clicking on the "Export archive" button. You will be asked a folder into which the ZIP 
archive will be created.

You can re-open this archive whenever you want. To do so, you have to:

- Open Napari and open the "Nucleation analyzer" plugin.
- Open again the image that you analyzed
- Drag and drop the ZIP archive into the Napari viewer.

It should restore the whole process **for visualization only**, you cannot use that as a 
starting point to continue the analysis. 