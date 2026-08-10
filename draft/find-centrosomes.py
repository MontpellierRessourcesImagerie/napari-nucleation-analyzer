import xarray as xr
import pandas as pd
import numpy as np
import trackpy as tp
from skimage.morphology import h_maxima
from concurrent.futures import ThreadPoolExecutor
from scipy.ndimage import (
    grey_opening,
    gaussian_laplace,
    median_filter
)


class FindCentrosomesOperator:
    def __init__(self):
        tp.quiet(True)
        self.input_image = None
        self.hints = {}
        self.centrosomes = pd.DataFrame(columns=[
            "centriole_id", 
            "T", 
            "Y", 
            "X", 
            "centrosome_id"
        ])
        self.prominence = FindCentrosomesOperator.default_prominence()
        self.searching_range = FindCentrosomesOperator.default_searching_range()
        self.memory = FindCentrosomesOperator.default_memory()
        self.max_binding_distance = FindCentrosomesOperator.default_max_binding_distance()

    @staticmethod
    def default_prominence():
        return 0.15

    @staticmethod
    def default_searching_range():
        return 3.25

    def get_searching_range_pxl(self):
        if self.input_image is None:
            raise ValueError("Input image must be set before computing searching range in pixels.")
        return self.searching_range / self.input_image.attrs["scale"]['X']

    @staticmethod
    def default_memory():
        return 10

    def get_memory_frames(self):
        if self.input_image is None:
            raise ValueError("Input image must be set before computing memory in frames.")
        return int(self.memory / self.input_image.attrs["scale"]['T'])

    @staticmethod
    def default_max_binding_distance():
        return 0.6

    def get_max_binding_distance_pxl(self):
        if self.input_image is None:
            raise ValueError("Input image must be set before computing max binding distance in pixels.")
        return self.max_binding_distance / self.input_image.attrs["scale"]['X']

    def set_input_image(self, image_arr, calibration, units):
        if image_arr.ndim != 3:
            raise ValueError("Input image must be a 3D array (T, Y, X).")
        self.input_image = xr.DataArray(
            image_arr, 
            dims=["T", "Y", "X"],
            attrs={
                "scale": calibration,
                "units": units
            }
        )

    def set_hints(self, hints):
        if self.input_image is None:
            raise ValueError("Input image must be set before setting hint points.")
        for centrosome_id, points_info in hints.items():
            if not isinstance(points_info, dict):
                raise ValueError(f"Hint points for centrosome {centrosome_id} must be a dictionary.")
            if "start" not in points_info or "end" not in points_info or "points" not in points_info:
                raise ValueError(f"Hint points for centrosome {centrosome_id} must contain 'start', 'end', and 'points' keys.")
            if not isinstance(points_info["points"], np.ndarray) or points_info["points"].shape[1] != 2:
                raise ValueError(f"'points' for centrosome {centrosome_id} must be a 2D numpy array with shape (N, 2).")
            if points_info["start"] < 0 or points_info["end"] >= self.input_image.sizes["T"]:
                raise ValueError(f"'start' and 'end' for centrosome {centrosome_id} must be within the time range of the input image (0 to {self.input_image.sizes['T'] - 1}).")
        self.hints = hints

    def get_centrosomes(self):
        if self.centrosomes is None:
            raise ValueError("Centrosomes have not been computed yet. Run the operator first.")
        return self.centrosomes

    def _tracking_pre_processing(self):
        if self.input_image is None:
            raise ValueError("Input image must be set before running the operator.")
        raw_image = self.input_image.transpose("T", "Y", "X").values.astype(np.float32)
        log_res = np.zeros_like(raw_image)
        print("Running preprocessing...")

        def _apply_for_frame(t):
            median = median_filter(raw_image[t], size=(3, 3))
            background = grey_opening(median, size=(20, 20))
            no_bg = median - background
            no_bg[no_bg < 0] = 0.0
            log_res[t] = gaussian_laplace(no_bg, sigma=(5, 5))

        with ThreadPoolExecutor() as executor:
            executor.map(_apply_for_frame, range(raw_image.shape[0]))
    
        log_res -= np.min(log_res)
        log_res /= np.max(log_res)
        log_res = log_res * -1.0 + 1.0
        return log_res

    def _find_maximas(self, log_res):
        print("Searching for local maxima...")

        def _process_frame(t):
            coords = h_maxima(log_res[t], h=self.prominence)
            coords = np.argwhere(coords)
            return [{"T": t, "Y": c[0], "X": c[1]} for c in coords]

        with ThreadPoolExecutor() as executor:
            results = list(executor.map(_process_frame, range(log_res.shape[0])))

        df_buffer = [item for sublist in results for item in sublist]
        return pd.DataFrame(df_buffer)
    
    def _track_spots(self, spots_df):
        search_range = self.get_searching_range_pxl()
        memory = self.get_memory_frames()
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
        spots_df["centriole_id"] = linked["particle"].astype(int) + 1
        return spots_df
    
    def _bind_tracks_to_hints(self, tracked):
        bindings = {}
        binding_dist = self.get_max_binding_distance_pxl()
        for centrosome_id, hint_info in self.hints.items():
            bindings[centrosome_id] = [None, None]
            start_frame = hint_info["start"]
            hint_points = hint_info["points"]
            candidates = tracked[tracked["T"] == start_frame]
            for i, hint_point in enumerate(hint_points):
                best_candidate = None
                min_distance = float('inf')
                for _, candidate in candidates.iterrows():
                    distance = np.linalg.norm(candidate[["Y", "X"]].values - hint_point)
                    if distance < min_distance and distance <= binding_dist:
                        min_distance = distance
                        best_candidate = candidate
                if best_candidate is not None:
                    bindings[centrosome_id][i] = best_candidate["centriole_id"]
        return bindings
    
    def _filter_by_hint_points(self, tracked):
        bindings = self._bind_tracks_to_hints(tracked)
        tracked = tracked.copy()
        tracked["centrosome_id"] = 0
        
        for centrosome_id, centriole_ids in bindings.items():
            t1, t2 = centriole_ids
            if t1 is not None:
                tracked.loc[tracked["centriole_id"] == t1, "centrosome_id"] = centrosome_id
            if t2 is not None:
                tracked.loc[tracked["centriole_id"] == t2, "centrosome_id"] = centrosome_id
            if t1 is None or t2 is None:
                raise ValueError(f"Failed to bind centrosomes to hints for centrosome {centrosome_id}.")

        tracked = tracked[tracked["centrosome_id"] != 0]
        return tracked
    
    def _interpolate_missing_time_points(self, tracked):
        """
        All the tracks must be complete in terms of locations.
        This function filss the gap by creating new rows for missing time points.
        The nex coordinates are interpolated linearly between the two closest known points.
        """
        complete_tracks = []
        for centriole_id, group in tracked.groupby("centriole_id"):
            minT = group["T"].min()
            maxT = group["T"].max()
            group = group.set_index("T").reindex(range(minT, maxT + 1))
            group["centriole_id"] = centriole_id
            group[["Y", "X"]] = group[["Y", "X"]].interpolate(method='linear')
            group["centrosome_id"] = group["centrosome_id"].ffill().bfill()
            complete_tracks.append(group.reset_index())
        df = pd.concat(complete_tracks, ignore_index=True)
        for centrosome_id, hint_info in self.hints.items():
            start_frame = hint_info["start"]
            end_frame = hint_info["end"]
            df = df[~((df["centrosome_id"] == centrosome_id) & ((df["T"] < start_frame) | (df["T"] > end_frame)))]
        return df

    def run(self):
        if self.input_image is None:
            raise ValueError("Input image must be set before running the operator.")
        if len(self.hints) == 0:
            raise ValueError("Hint points must be set before running the operator.")
        
        preprocessed = self._tracking_pre_processing()
        maximas = self._find_maximas(preprocessed)

        tracked = self._track_spots(maximas)
        tracked = self._filter_by_hint_points(tracked)
        tracked = self._interpolate_missing_time_points(tracked)
        
        self.centrosomes = tracked.sort_values(by=["centrosome_id", "centriole_id", "T"]).reset_index(drop=True)


def launch_full_process():
    from pathlib import Path
    import tifffile as tiff
    import napari

    def make_control_image(original, maximas):
        control = np.zeros_like(original, dtype=np.uint16)
        for _, row in maximas.iterrows():
            control[int(row["T"]), int(row["Y"]), int(row["X"])] = int(row["centriole_id"])
        return control

    DUMP = Path("/home/clement/Documents/projects/nucleation/draft/implementation/dump")

    folder_in = Path("/home/clement/Documents/projects/nucleation/3VPCs")
    filename = "251119_#4_30_001_016.vsi - C561.tif"
    path_in = folder_in / filename

    calib = {'T': 1, 'Y': 0.1083333, 'X': 0.1083333}
    units = {'T': 's', 'Y': 'um', 'X': 'um'}

    img = tiff.imread(path_in) # (T, Y, X)
    hints = {
        1: {
            'start': 17,
            'end': 463,
            'points': np.array([
                [347, 90],
                [337, 132]
            ]),
            'color': "#1f77b4"
        },
        2: {
            'start': 17,
            'end': 463,
            'points': np.array([
                [248, 477],
                [237, 514]
            ]),
            'color': "#ff7f0e"
        }
    }

    operator = FindCentrosomesOperator()
    operator.set_input_image(img, calib, units)
    operator.set_hints(hints)
    operator.run()

    centrosomes = operator.get_centrosomes()
    centrosomes.to_csv(DUMP / "centrosome_points_track.csv", index=False)

    viewer = napari.Viewer()
    viewer.add_image(img, name="raw")
    viewer.add_points(
        centrosomes[['T', 'Y', 'X']].values,
        name="maximas",
        size=7,
        face_color="transparent",
        border_color="red",
        opacity=0.75,
        visible=True
    )
    viewer.add_tracks(
        centrosomes[["centriole_id", "T", "Y", "X"]].values,
        name="tracks",
        features=centrosomes
    )
    napari.run()


if __name__ == "__main__":
    launch_full_process()