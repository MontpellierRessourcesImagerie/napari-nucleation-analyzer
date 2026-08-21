# napari-nucleation-analyzer

![GitHub License](https://img.shields.io/github/license/MontpellierRessourcesImagerie/napari-nucleation-analyzer)
![Python Version](https://img.shields.io/badge/Python-3.10|3.11|3.12-blue?logo=python)
![Unit tests](https://img.shields.io/github/actions/workflow/status/MontpellierRessourcesImagerie/napari-nucleation-analyzer/test_and_deploy.yml?logo=pytest&label=tests)
![PyPI version](https://img.shields.io/pypi/v/napari-nucleation-analyzer)

This Napari plugin allows you to:
- Detect pairs of centrioles and track them over time to form centrosomes.
- Build arcs of desired radius and angle aligned with the vector pointing towards the 
  other centriole of the pair.
- Build kymographs of these arcs.
- Find spots on the kymographs.
- Produce a CSV file containing a column per centrosome and a line per time point 
  giving the number of spots at each T.

This README is very short, for further information, refer yourself to the [documentation](https://montpellierressourcesimagerie.github.io/napari-nucleation-analyzer/index.html).

If you run in any trouble, you can [❗open a ticket](https://github.com/MontpellierRessourcesImagerie/napari-nucleation-analyzer/issues).
