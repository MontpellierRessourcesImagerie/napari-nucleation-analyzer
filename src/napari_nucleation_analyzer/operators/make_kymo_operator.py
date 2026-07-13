import numpy as np
from scipy.ndimage import median_filter

class MakeKymographOperator:
    
    def __init__(self):
        self.kymographs  = None
        self.arcs        = None
        self.centrosomes = None
        self.input_image = None
        self.shift       = 0.5

    def set_arcs(self, arcs):
        self.arcs = arcs

    def set_centrosomes(self, centrosomes):
        self.centrosomes = centrosomes

    def set_input_image(self, input_image):
        self.input_image = input_image

    def set_shift(self, shift):
        self.shift = shift

    def make_kymograph(self, tracked_points, arcs, intensities, shift=0.5):
        """
        This function creates a kymograph for each track ID.
        The values are measured along the arc from the 'args' parameter.
        For a better smoothing, the vector made from (VY, VX) is used to take a few pixels before and after the point of the arc.
        The kymographs are 2D images having a shape of (time, arc_points). The values are the mean intensity along the arc.
        """
        if len(arcs) == 0:
            raise ValueError("No arcs provided for kymograph generation.")
        
        kymographs = {}
        intensities = median_filter(intensities, size=(1, 3, 3))

        for track_id in arcs.keys():
            arc = arcs[track_id]
            centrosome = tracked_points[tracked_points["track_id"] == track_id].sort_values("T")[["Y", "X"]].values
            
            dir_vector = arc - centrosome[:, np.newaxis, :]
            dir_vector /= np.linalg.norm(dir_vector, axis=2, keepdims=True)

            coefs = np.arange(1-shift, 1+shift + 0.5, 0.5)  # coefs for sampling along the direction vector
            # (T, num_coefs, num_points, 2)
            sample_pts = arc[:, np.newaxis, :, :] + coefs[np.newaxis, :, np.newaxis, np.newaxis] * dir_vector[:, np.newaxis, :, :]

            sample_pxls = np.round(sample_pts).astype(int)  # (T, num_coefs, num_points, 2)

            # Split into separate y/x index arrays and clip to stay inside the image bounds
            T, Y, X = intensities.shape
            ys = np.clip(sample_pxls[..., 0], 0, Y - 1)  # (T, num_coefs, num_points)
            xs = np.clip(sample_pxls[..., 1], 0, X - 1)  # (T, num_coefs, num_points)
            t_idx = np.arange(T)[:, np.newaxis, np.newaxis]  # (T, 1, 1) -> broadcasts against ys/xs

            pxl_values = intensities[t_idx, ys, xs]
            kymograph = np.mean(pxl_values, axis=1)
            kymographs[track_id] = kymograph

        return kymographs
    
    def get_kymographs(self):
        if self.kymographs is None:
            raise ValueError("Kymographs have not been generated yet. Please run the operator first.")
        return self.kymographs

    def run(self):
        if self.arcs is None or self.centrosomes is None or self.input_image is None:
            raise ValueError("Arcs, centrosomes, and input image must be set before running the operator.")
        
        self.kymographs = self.make_kymograph(
            self.centrosomes, 
            self.arcs, 
            self.input_image, 
            self.shift
        )


if __name__ == "__main__":
    import tifffile as tiff
    from pathlib import Path
    import pandas as pd

    folder_in = Path("/home/clement/Documents/projects/nucleation/3VPCs")
    filename  = "251119_#4_30_001_016.vsi - C561.tif"
    path_in   = folder_in / filename
    image     = tiff.imread(path_in)
    t_start   = 17
    t_end     = 463

    arcs = {}
    arcs_folder = Path("/home/clement/Documents/projects/nucleation/draft/implementation/dump/arcs")
    arcs_content = [f for f in arcs_folder.iterdir() if f.is_file() and f.suffix == ".npy"]
    for arc_file in arcs_content:
        track_id = int(arc_file.name.replace(".npy", ""))
        arc = np.load(arc_file)
        arcs[track_id] = arc

    centros_folder = Path("/home/clement/Documents/projects/nucleation/draft/implementation/dump")
    centro_path = centros_folder / "centrosome_points_track.csv"
    centrosomes = pd.read_csv(centro_path)

    op = MakeKymographOperator()
    op.set_input_image(image[t_start:t_end+1])
    op.set_arcs(arcs)
    op.set_centrosomes(centrosomes)
    op.run()

    kymos = op.get_kymographs()
    out_folder = Path("/home/clement/Documents/projects/nucleation/draft/implementation/dump/kymos")
    out_folder.mkdir(exist_ok=True)

    for track_id, kymo in kymos.items():
        out_path = out_folder / f"kymo_{track_id}.tif"
        print(track_id, kymo.shape, out_path)
        tiff.imwrite(out_path, kymo.astype(np.float32))