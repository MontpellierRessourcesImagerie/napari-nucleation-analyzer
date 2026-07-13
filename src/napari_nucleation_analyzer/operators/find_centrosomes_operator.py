import xarray as xr
import pandas as pd
from tqdm import tqdm
import numpy as np
from scipy.ndimage import gaussian_filter
from skimage.feature import peak_local_max
import trackpy as tp


class FindCentrosomesOperator:

    def __init__(self):
        self.input_image = None
        self.calibration = None
        self.time_range  = None
        self.hint_points = None
        self.centrosomes = None

    def set_input_image(self, image_arr, calibration):
        if image_arr.ndim != 3:
            raise ValueError("Input image must be a 3D array (T, Y, X).")
        self.input_image = xr.DataArray(image_arr, dims=["T", "Y", "X"])
        self.calibration = calibration

    def set_time_range(self, start, end):
        if self.input_image is None:
            raise ValueError("Input image must be set before setting time range.")
        if start < 0:
            raise ValueError("Start time must be non-negative.")
        if end >= self.input_image.sizes["T"]:
            raise ValueError(f"End time must be less than {self.input_image.sizes['T']}.")
        if start > end:
            raise ValueError("Start time must be less than or equal to end time.")
        self.time_range = (start, end)

    def set_hint_points(self, hint_points):
        if not isinstance(hint_points, np.ndarray):
            raise ValueError("Hint points must be a numpy array.")
        if hint_points.ndim != 2 or hint_points.shape[1] not in [2, 3]:
            raise ValueError("Hint points must be a 2D array with shape (N, 2) or (N, 3) for (T, Y, X).")
        if hint_points.shape[1] == 2:
            self.hint_points = hint_points
        else:
            self.hint_points = hint_points[:, 1:]

    def get_centrosomes(self):
        if self.centrosomes is None:
            raise ValueError("Centrosomes have not been computed yet. Run the operator first.")
        return self.centrosomes

    def _tracking_pre_processing(self, img):
        if self.time_range is None:
            raise ValueError("Time range must be set before running the operator.")
        img = img.astype(np.float32)
        wide_gaussian = gaussian_filter(img, sigma=(0, 8, 8))
        small_gaussian = gaussian_filter(img, sigma=(0, 3, 3))
        img = small_gaussian - wide_gaussian
        print("Preprocessing done.")
        return img[self.time_range[0]:self.time_range[1]+1]
    
    def _find_maximas(self, img, num_peaks):
        print("Searching for local maxima...")
        norm = (img - img.min()) / (img.max() - img.min())
        rows = []
        for t in range(norm.shape[0]):
            coords = peak_local_max(
                norm[t],
                min_distance=2,
                num_peaks=num_peaks
            )
            for coord in coords:
                rows.append({
                    "Y": coord[0],
                    "X": coord[1],
                    "T": int(t)
                })
        return pd.DataFrame(rows)
    
    def _track_spots(self, spots_df):
        search_range = 30
        memory = 4
        predictor = tp.predict.NearestVelocityPredict()
        linked = predictor.link_df(
                spots_df,
                search_range=search_range,
                memory=memory,
                pos_columns=["Y", "X"],
                t_column='T',
                adaptive_stop=0.01,
                adaptive_step=0.75
            )
        spots_df["track_id"] = linked["particle"].astype(int) + 1
        return spots_df
    
    def _filter_track_length(self, tracked, time_thr=10):
        track_lengths = tracked.groupby("track_id").size()
        valid_tracks = track_lengths[track_lengths >= time_thr].index
        filtered_tracked = tracked[tracked["track_id"].isin(valid_tracks)].copy()
        return filtered_tracked
    
    def _interpolate_missing_time_points(self, tracked):
        """
        All the tracks must be complete in terms of locations.
        This function filss the gap by creating new rows for missing time points.
        The nex coordinates are interpolated linearly between the two closest known points.
        """
        complete_tracks = []
        minT = tracked["T"].min()
        maxT = tracked["T"].max()
        for track_id, group in tracked.groupby("track_id"):
            group = group.set_index("T").reindex(range(minT, maxT + 1))
            group["track_id"] = track_id
            group[["Y", "X"]] = group[["Y", "X"]].interpolate(method='linear')
            complete_tracks.append(group.reset_index())
        df = pd.concat(complete_tracks, ignore_index=True)
        df["track_id"] = df["track_id"].astype(int)
        return df

    def run(self):
        if self.input_image is None:
            raise ValueError("Input image must be set before running the operator.")
        if self.time_range is None:
            raise ValueError("Time range must be set before running the operator.")
        if self.hint_points is None:
            raise ValueError("Hint points must be set before running the operator.")
        
        preprocessed = self._tracking_pre_processing(self.input_image.transpose("T", "Y", "X").values)
        
        maximas = self._find_maximas(
            preprocessed, 
            num_peaks=len(self.hint_points)
        )

        tracked = self._track_spots(maximas)
        tracked = self._filter_track_length(tracked)
        tracked = self._interpolate_missing_time_points(tracked)
        tracked['T'] += self.time_range[0]
        
        self.centrosomes = tracked.sort_values(by=["track_id", "T"]).reset_index(drop=True)


if __name__ == "__main__":
    from pathlib import Path
    import tifffile as tiff

    def make_control_image(original, maximas):
        control = np.zeros_like(original, dtype=np.uint16)
        for _, row in maximas.iterrows():
            control[int(row["T"]), int(row["Y"]), int(row["X"])] = int(row["track_id"])
        return control

    DUMP = Path("/home/clement/Documents/projects/nucleation/draft/implementation/dump")

    folder_in = Path("/home/clement/Documents/projects/nucleation/3VPCs")
    filename  = "251119_#4_30_001_016.vsi - C561.tif"
    path_in   = folder_in / filename
    calib     = (1, 0.1083333, 0.1083333)

    img       = tiff.imread(path_in) # T, Y, X
    t_start   = 17
    t_end     = 463
    intensity = img

    starters = np.array([
        [17, 347, 90],
        [17, 337, 132],
        [17, 248, 477],
        [17, 237, 514]
    ])

    operator = FindCentrosomesOperator()
    operator.set_input_image(intensity, calib)
    operator.set_time_range(t_start, t_end)
    operator.set_hint_points(starters)
    operator.run()

    tracked_centrosomes = operator.get_centrosomes()
    tracked_centrosomes.to_csv(DUMP / "centrosome_points_track.csv", index=False)

    control = make_control_image(intensity, tracked_centrosomes)
    tiff.imwrite(DUMP / "control.tif", control)