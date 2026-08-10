import numpy as np
from skimage.feature import peak_local_max
from scipy.ndimage import gaussian_laplace
import pandas as pd

class FindSpotsOperator:

    def __init__(self):
        self.kymographs = None
        self.coordinates = None
        self.std_coeff = self.default_std_coeff()

    @staticmethod
    def default_std_coeff():
        return 1.5

    def set_kymographs(self, kymographs):
        self.kymographs = kymographs

    def set_std_coeff(self, coeff):
        if coeff <= 0:
            raise ValueError("Standard deviation coefficient must be positive.")
        self.std_coeff = coeff

    def get_coordinates(self):
        if self.coordinates is None:
            raise ValueError("Coordinates not computed. Please run the operator first.")
        return self.coordinates

    def _find_spots(self, kymographs):
        """
        Takes the dict of kymographs.
        For Each kymograph, find the maximas corresponding to spots.
        Write them in a buffer having the same size as the kymo and make a sum projection of it to keep only the time axis.
        Make a new dataframe with each column corresponding to a track id.
        """
        buffer = {}
        for centriole_id, kymo in kymographs.items():
            kymo = gaussian_laplace(kymo, sigma=2)
            kymo *= -1
            coordinates = peak_local_max(
                kymo,
                min_distance=2,
                threshold_abs=np.std(kymo) * self.std_coeff
            )
            buffer[centriole_id] = coordinates
        return buffer
    
    def run(self):
        if self.kymographs is None:
            raise ValueError("Kymographs not set. Please set kymographs before running the operator.")
        self.coordinates = self._find_spots(self.kymographs)


if __name__ == "__main__":
    from pathlib import Path
    import tifffile as tiff
    import napari

    # Loading the produced kymographs
    kymographs = {}
    kymos_path = Path("/home/clement/Documents/projects/nucleation/draft/implementation/dump/kymos")

    for kymo_path in kymos_path.glob("kymo_*.tif"):
        centriole_id = int(kymo_path.name.replace("kymo_", "").replace(".tif", ""))
        kymograph = tiff.imread(kymo_path)
        kymographs[centriole_id] = kymograph

    # Running the operator
    op = FindSpotsOperator()
    op.set_kymographs(kymographs)
    op.run()

    # Showing the result in Napari
    coordinates = op.get_coordinates()
    viewer = napari.Viewer()
    translation = 0

    for centriole_id, coords in coordinates.items():
        kymo = kymographs[centriole_id]
        viewer.add_image(
            kymo, 
            name=f"Kymograph {centriole_id}",
            translate=[0, translation]
        )
        viewer.add_points(
            coords, 
            name=f"Spots {centriole_id}", 
            size=5, 
            face_color='transparent',
            border_color='red',
            translate=[0, translation]
        )
        translation += (kymo.shape[1] + 10)
        np.save(kymos_path / f"spots_{centriole_id}.npy", coords)

    napari.run()