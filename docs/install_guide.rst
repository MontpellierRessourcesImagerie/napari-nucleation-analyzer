=========================
Quick start: A user guide
=========================

1. Installation
===============

You have three choices to install NA, but all of them require you to install something to deal 
with virtual Python environments. In the following instruction, we will use Miniconda (https://repo.anaconda.com/miniconda/)
but you are free to use the one you are the most comfortable with.

Create a new environment with Python 3.10 or higher, and activate it. Then you can install Napari
with the command: :code:`pip install napari[all]`. You can check that Napari is correctly installed 
by running :code:`napari` in your terminal.

a. Install a frozen dev version
-------------------------------

- Start by going on the GitHub repository of :code:`napari-nucleation-analyzer` (https://github.com/MontpellierRessourcesImagerie/napari-nucleation-analyzer).
- Click on the green button :code:`Code` and click on :code:`Download ZIP`.
- Decompress the ZIP file and move its content to somewhere you won't delete it by mistake.
- In your terminal with the environment activated, run the command :code:`pip install -e path/to/napari-nucleation-analyzer` 
  where :code:`path/to/napari-nucleation-analyzer` is the path to the folder you just downloaded. Depending on your system,
  you can simply drag and drop the folder in your terminal to get the path.

b. Install an updatable dev version
-----------------------------------

- If you wish to go with this option, you need to have :code:`git` installed on your system. 
  If you don't have it, you can download it from https://git-scm.com/downloads. To check if you have it, 
  run the command :code:`git --version` in any terminal.
- In your terminal with the environment activated, run the command :code:`pip install git+https://github.com/MontpellierRessourcesImagerie/napari-nucleation-analyzer.git`

c. Install the stable version from PyPI
---------------------------------------

- In your terminal with the environment activated, run the command :code:`pip install napari-nucleation-analyzer`


2. Final check
==============

- To make sure that NA is correctly installed, you can run :code:`napari` in your terminal. 
- Then click on the "Plugins" menu in the top-bar manu of Napari. In there you should see:

  - Calibration Tool > Scale Tool
  - Nucleation analyzer

- If this is the case, everything is correctly installed and you can start using NA.
