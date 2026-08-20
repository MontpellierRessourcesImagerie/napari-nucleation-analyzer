import numpy as np
import pandas as pd
import xarray as xr
from scipy.ndimage import median_filter, map_coordinates
from concurrent.futures import ThreadPoolExecutor


class MakeKymographOperator:
    """
    Takes a collection of arcs and build kymographs for each of them.
    To attenuate the importance of noise, the mean intensity is processed for each point of the arc
    along the local normal vector. The mean is computed over a range of pixels defined by the 'mean_range' parameter.
    The kymographs are 2D images having a shape of (time, arc_points). The values are the mean intensity along the arc.
    """
    def __init__(self):
        self.arcs        = None
        self.pairs = None
        self.input_image = None

        self.kymographs  = {}

    def set_arcs(self, arcs: dict):
        for centrosome_id, arc in arcs.items():
            if arc.ndim != 3 or arc.shape[2] != 3:
                raise ValueError(f"Arc for centrosome {centrosome_id} must have shape (N, 3).")
        self.arcs = arcs

    def set_pairs(self, pairs: pd.DataFrame):
        required_columns = {"centrosome_id", "T", "Y", "X", "pair_id"}
        if not required_columns.issubset(pairs.columns):
            raise ValueError(f"Pairs dataframe must contain columns: {required_columns}")
        self.pairs = pairs

    def set_input_image(self, image_arr: np.ndarray, calibration: dict, units: dict):
        if image_arr.ndim != 3:
            raise ValueError("Input image must be a 3D array (T, Y, X).")
        self.input_image = xr.DataArray(
            image_arr,
            dims=["T", "Y", "X"],
            attrs={"scale": calibration, "units": units},
        )

    def get_input_image(self) -> xr.DataArray:
        if self.input_image is None:
            raise ValueError("Input image has not been set.")
        return self.input_image

    def get_kymographs(self):
        if self.kymographs is None:
            raise ValueError("Kymographs have not been generated yet. Please run the operator first.")
        return self.kymographs

    def _make_kymograph(self, tracked_points, arcs, intensities):
        """
        This function creates a kymograph for each track ID.
        The values are measured along the arc from the 'args' parameter.
        For a better smoothing, the vector made from (VY, VX) is used to take a few pixels before and after the point of the arc.
        The kymographs are 2D images having a shape of (time, arc_points). The values are the mean intensity along the arc.
        """
        if len(arcs) == 0:
            raise ValueError("No arcs provided for kymograph generation.")

        img = self.get_input_image()
        kymographs = {}
        intensities = np.zeros_like(img.values)
        img_axes = [str(ax) for ax in img.dims]

        def median_t(t):
            intensities[t] = median_filter(img.values[t], size=(3, 3))

        with ThreadPoolExecutor() as executor:
            executor.map(median_t, range(intensities.shape[0]))

        for centrosome_id in arcs.keys():
            arc = arcs[centrosome_id]
            centrosome = tracked_points[tracked_points["centrosome_id"] == centrosome_id].sort_values("T")[img_axes].values

            # Arc: (T, num_points, 3)
            # Centrosome: (T, 3)
            dir_vector = arc - centrosome[:, np.newaxis, :] # new axis for broadcasting, shape: (T, num_points, 3)
            dir_vector /= np.linalg.norm(dir_vector, axis=2, keepdims=True)

            coefs = np.array([-1.0, -0.5, 0.0, 0.5, 1.0])
            # Samples: (T, num_coefs, num_points, 3)
            sample_pts = arc[:, np.newaxis, :, :] + coefs[np.newaxis, :, np.newaxis, np.newaxis] * dir_vector[:, np.newaxis, :, :]

            # Linearize the coordinates to have a long list to feed to map_coordinates
            original_shape = sample_pts.shape
            stack = sample_pts.transpose(3, 0, 1, 2).reshape(3, -1)  # shape (3, T * Num coefs * Num points)
            # Reading the interpolated values
            values = map_coordinates(intensities, stack, order=1, prefilter=False) # shape: (T * Num coefs * Num points,)
            pxl_values = values.reshape(*original_shape[:-1]) # shape: (T, num_coefs, num_points)

            kymograph = np.mean(pxl_values, axis=1)
            kymographs[centrosome_id] = kymograph

        return kymographs

    def run(self):
        if (self.arcs is None) or (self.pairs is None) or (self.input_image is None):
            raise ValueError("Arcs, pairs, and input image must be set before running the operator.")
        
        self.kymographs = self._make_kymograph(
            self.pairs, 
            self.arcs, 
            self.input_image
        )


if __name__ == "__main__":
    import tifffile as tiff
    from pathlib import Path
    import pandas as pd

    # Importing the image to work on
    folder_img = Path("/home/clement/Documents/projects/nucleation/3VPCs")
    filename = "251119_#4_30_001_016.vsi - C561.tif"
    path_in = folder_img / filename
    image = tiff.imread(path_in)
    calib = {'T': 1, 'Y': 0.1083333, 'X': 0.1083333}
    units = {'T': 's', 'Y': 'um', 'X': 'um'}

    # Retrieving the pairs dataframe
    centros_folder = Path("/home/clement/Documents/projects/nucleation/draft/implementation/dump")
    centro_path = centros_folder / "pairs_arcs.csv"
    pairs = pd.read_csv(centro_path)

    # Retrieving the processed arcs
    arcs = {}
    arcs_folder = Path("/home/clement/Documents/projects/nucleation/draft/implementation/dump")
    for centrosome_id in pairs["centrosome_id"].unique():
        arc_path = arcs_folder / f"arc_{centrosome_id}.npy"
        if arc_path.exists():
            arcs[centrosome_id] = np.load(arc_path)
        else:
            print(f"Warning: Arc file for centrosome {centrosome_id} not found at {arc_path}")

    # Running the operator to generate kymographs
    op = MakeKymographOperator()
    op.set_input_image(image, calib, units)
    op.set_arcs(arcs)
    op.set_pairs(pairs)
    op.run()

    # Saving the kymographs to disk
    kymos = op.get_kymographs()
    out_folder = Path("/home/clement/Documents/projects/nucleation/draft/implementation/dump/kymos")
    out_folder.mkdir(exist_ok=True)

    for centrosome_id, kymo in kymos.items():
        out_path = out_folder / f"kymo_{centrosome_id}.tif"
        print(centrosome_id, kymo.shape, out_path)
        tiff.imwrite(out_path, kymo.astype(np.float32))