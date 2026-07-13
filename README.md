# napari-nucleation-analyzer

This Napari plugin allows you to:
- Detect pairs of centrosomes and track them over time.
- Build arcs of desired radius and angle aligned with the vector pointing towards the 
  other centrosome of the pair.
- Build kymographs of these arcs.
- Find spots on the kymographs.
- Produce a CSV file containing a column per centrosome and a line per time point 
  giving the number of spots at each T.