import numpy as np
import tifffile as tiff
from skimage.feature import peak_local_max
from scipy.ndimage import gaussian_laplace
import pandas as pd

class FindSpotsOperator:

    def __init__(self):
        self.kymographs = None
        self.coordinates = None

    def set_kymographs(self, kymographs):
        self.kymographs = kymographs

    def _find_spots(self, kymographs):
        """
        Takes the dict of kymographs.
        For Each kymograph, find the maximas corresponding to spots.
        Write them in a buffer having the same size as the kymo and make a sum projection of it to keep only the time axis.
        Make a new dataframe with each column corresponding to a track id.
        """
        summed_tracks = {}
        buffer = {}
        for track_id, kymo in kymographs.items():
            kymo = gaussian_laplace(kymo, sigma=2)
            kymo *= -1
            coordinates = peak_local_max(
                kymo,
                min_distance=2,
                threshold_abs=np.std(kymo) * 1.5
            )
            buffer[track_id] = coordinates
        return buffer
    
    def get_coordinates(self):
        if self.coordinates is None:
            raise ValueError("Coordinates not computed. Please run the operator first.")
        return self.coordinates
    
    def run(self):
        if self.kymographs is None:
            raise ValueError("Kymographs not set. Please set kymographs before running the operator.")
        self.coordinates = self._find_spots(self.kymographs)


if __name__ == "__main__":
    from pathlib import Path

    kymographs = {}
    kymos_path = Path("/home/clement/Documents/projects/nucleation/draft/implementation/dump/kymos")
    for kymo_path in kymos_path.glob("kymo_*.tif"):
        track_id = int(kymo_path.name.replace("kymo_", "").replace(".tif", ""))
        kymograph = tiff.imread(kymo_path)
        kymographs[track_id] = kymograph

    op = FindSpotsOperator()
    op.set_kymographs(kymographs)
    op.run()