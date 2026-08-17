import xarray as xr
import pandas as pd
import numpy as np
import trackpy as tp
from skimage.morphology import h_maxima
from concurrent.futures import ThreadPoolExecutor
from scipy.ndimage import grey_opening, gaussian_laplace, median_filter
from pprint import pprint


class FindCentrosomesOperator:
    def __init__(self):
        self.input_image = None
        self.hints = {}
        self.prominence = self.default_prominence()
        self.searching_range = self.default_searching_range()
        self.memory = self.default_memory()
        self.max_binding_distance = self.default_max_binding_distance()
        self.centrosomes = None
        tp.quiet(True)

    @staticmethod
    def default_prominence() -> float:
        return 10.0

    @staticmethod
    def default_searching_range() -> float:
        return 3.25

    @staticmethod
    def default_memory() -> int:
        return 10

    @staticmethod
    def default_max_binding_distance() -> float:
        return 0.6

    def set_input_image(self, image_arr: np.ndarray, calibration: dict, units: dict):
        if image_arr.ndim != 3:
            raise ValueError("Input image must be a 3D array (T, Y, X).")
        self.input_image = xr.DataArray(
            image_arr,
            dims=["T", "Y", "X"],
            attrs={"scale": calibration, "units": units},
        )

    def set_hints(self, hints: dict):
        if self.input_image is None:
            raise ValueError("Input image must be set before setting hint points.")
        for centrosome_id, points_info in hints.items():
            if len(set(points_info.keys()).intersection({"start", "end", "points"})) != 3:
                raise ValueError(
                    f"[{centrosome_id}]: 'start', 'end', and 'points' keys are required."
                )
            if points_info["points"].ndim != 2:
                raise ValueError(
                    f"[{centrosome_id}]: points must be a 2D numpy array with shape (N, 2)."
                )
            if (points_info["points"].shape[1] != 2):
                raise ValueError(
                    f"[{centrosome_id}]: points must be a 2D numpy array with shape (N, 2)."
                )
            if (
                points_info["start"] < 0
                or points_info["end"] >= self.input_image.sizes["T"]
            ):
                raise ValueError(
                    f"'start' and 'end' for centrosome {centrosome_id} must be within the time range of the input image (0 to {self.input_image.sizes['T'] - 1})."
                )
        self.hints = hints

    def set_prominence(self, prominence: float):
        if prominence <= 0:
            raise ValueError("Prominence must be a strictly positive value.")
        self.prominence = prominence

    def set_searching_range(self, searching_range: float):
        if searching_range <= 0:
            raise ValueError("Searching range must be a positive value.")
        self.searching_range = searching_range

    def set_memory(self, memory: float | int):
        if memory < 0:
            raise ValueError("Memory must be a non-negative value.")
        self.memory = memory

    def set_max_binding_distance(self, max_binding_distance: float):
        if max_binding_distance <= 0:
            raise ValueError("Max binding distance must be a positive value.")
        self.max_binding_distance = max_binding_distance

    def get_input_image(self) -> xr.DataArray:
        if self.input_image is None:
            raise ValueError("Input image has not been set.")
        return self.input_image

    def get_hints(self) -> dict:
        if not self.hints:
            raise ValueError("Hint points have not been set.")
        return self.hints

    def get_prominence(self) -> float:
        return self.prominence

    def get_searching_range(self) -> float:
        return self.searching_range

    def get_searching_range_pxl(self) -> int:
        if self.input_image is None:
            raise ValueError(
                "Input image must be set before computing searching range in pixels."
            )
        return int(np.ceil(self.searching_range / self.input_image.attrs["scale"]["X"]))

    def get_memory(self) -> float | int:
        return self.memory

    def get_memory_frames(self) -> int:
        if self.input_image is None:
            raise ValueError(
                "Input image must be set before computing memory in frames."
            )
        return int(self.memory / self.input_image.attrs["scale"]["T"])

    def get_max_binding_distance(self) -> float:
        return self.max_binding_distance

    def get_max_binding_distance_pxl(self) -> int:
        if self.input_image is None:
            raise ValueError(
                "Input image must be set before computing max binding distance in pixels."
            )
        return int(
            np.ceil(self.max_binding_distance / self.input_image.attrs["scale"]["X"])
        )

    def get_centrosomes(self) -> pd.DataFrame:
        if self.centrosomes is None:
            raise ValueError(
                "Centrosomes have not been computed yet. Run the operator first."
            )
        return self.centrosomes

    def _tracking_pre_processing(self) -> np.ndarray:
        img = self.get_input_image()
        raw_image = img.transpose("T", "Y", "X").values.astype(np.float32)
        log_res = np.zeros_like(raw_image)

        def _apply_for_frame(t):
            median = median_filter(raw_image[t], size=(3, 3))
            background = grey_opening(median, size=(20, 20))
            no_bg = median - background
            no_bg[no_bg < 0] = 0.0
            log_res[t] = gaussian_laplace(no_bg, sigma=(5.0, 5.0))

        with ThreadPoolExecutor() as executor:
            executor.map(_apply_for_frame, range(raw_image.shape[0]))

        log_res -= np.min(log_res)
        log_res /= np.max(log_res)
        log_res = log_res * -1.0 + 1.0
        return log_res

    def _find_maximas(self, log_res) -> pd.DataFrame:
        std_dev = np.std(log_res)
        prom = std_dev * self.prominence
        print(f"Prominence for centrioles: {prom:.4f}")
        def _process_frame(t):
            coords = h_maxima(log_res[t], h=prom)
            coords = np.argwhere(coords)
            return [{"T": t, "Y": c[0], "X": c[1]} for c in coords]

        with ThreadPoolExecutor() as executor:
            results = list(executor.map(_process_frame, range(log_res.shape[0])))

        df_buffer = [item for sublist in results for item in sublist]
        return pd.DataFrame(df_buffer)

    def _track_spots(self, spots_df) -> pd.DataFrame:
        search_range = self.get_searching_range_pxl()
        memory = self.get_memory_frames()
        predictor = tp.predict.NearestVelocityPredict()
        linked = predictor.link_df(
            spots_df,
            search_range=search_range,
            memory=memory,
            pos_columns=["Y", "X"],
            t_column="T",
            adaptive_stop=0.01,
            adaptive_step=0.75,
        )
        spots_df["centriole_id"] = linked["particle"].astype(int) + 1
        return spots_df

    def _bind_tracks_to_hints(self, tracked) -> dict:
        bindings = {} # centrosome_id -> [centriole_id1, centriole_id2]
        binding_dist = self.get_max_binding_distance_pxl()
        used_centrioles = set()
        tracks_span = {
            centriole_id: group["T"].max() - group["T"].min() 
            for centriole_id, group in tracked.groupby("centriole_id")
        }
        print("Tracks span (in frames):")
        pprint(tracks_span)

        for centrosome_id, hint_info in self.hints.items():
            bindings[centrosome_id] = [None, None]
            start_frame = hint_info["start"]
            end_frame = hint_info["end"]
            hint_points = hint_info["points"]
            candidates = tracked[tracked["T"] == start_frame]
            
            for i, hint_point in enumerate(hint_points):
                # sorted tuples: (distance, centriole_id) for the starting frame
                sorted_points = sorted([
                    (np.linalg.norm(c[["Y", "X"]].values - hint_point), c['centriole_id']) 
                    for _, c in candidates.iterrows()
                ], key=lambda x: x[0])
                # filter: remove already used; remove too far away
                sorted_points = [
                    (dist, cid) 
                    for dist, cid in sorted_points 
                    if cid not in used_centrioles and dist <= binding_dist
                ]
                # candidate centriole id
                best_candidate = None
                # desired length of the track for this hint
                current_hint_duration = end_frame - start_frame
                
                for _, centriole_id in sorted_points:
                    if tracks_span[centriole_id] >= current_hint_duration:
                        best_candidate = centriole_id
                        break

                if best_candidate is not None:
                    used_centrioles.add(best_candidate)

                bindings[centrosome_id][i] = best_candidate
                
        return bindings

    def _find_centrosomes(self, tracked) -> pd.DataFrame:
        bindings = self._bind_tracks_to_hints(tracked)
        tracked = tracked.copy()
        tracked["centrosome_id"] = 0

        for centrosome_id, centriole_ids in bindings.items():
            t1, t2 = centriole_ids
            if t1 is None or t2 is None:
                raise ValueError(
                    f"Failed to bind centrosomes to hints for centrosome {centrosome_id}."
                )
            tracked.loc[tracked["centriole_id"] == t1, "centrosome_id"] = (
                centrosome_id
            )
            tracked.loc[tracked["centriole_id"] == t2, "centrosome_id"] = (
                centrosome_id
            )
            

        tracked = tracked[tracked["centrosome_id"] != 0]
        return tracked

    def _interpolate_missing_time_points(self, tracked) -> pd.DataFrame:
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
            group[["Y", "X"]] = group[["Y", "X"]].interpolate(method="linear")
            group["centrosome_id"] = group["centrosome_id"].ffill().bfill()
            complete_tracks.append(group.reset_index())
        df = pd.concat(complete_tracks, ignore_index=True)
        for centrosome_id, hint_info in self.hints.items():
            start_frame = hint_info["start"]
            end_frame = hint_info["end"]
            df = df[
                ~(
                    (df["centrosome_id"] == centrosome_id)
                    & ((df["T"] < start_frame) | (df["T"] > end_frame))
                )
            ]
        return df

    def run(self):
        if self.input_image is None:
            raise ValueError("Input image must be set before running the operator.")
        if len(self.hints) == 0:
            raise ValueError("Hint points must be set before running the operator.")

        preprocessed = self._tracking_pre_processing()
        maximas = self._find_maximas(preprocessed)

        tracked = self._track_spots(maximas)
        tracked = self._find_centrosomes(tracked)
        tracked = self._interpolate_missing_time_points(tracked)

        self.centrosomes = tracked.sort_values(
            by=["centrosome_id", "centriole_id", "T"]
        ).reset_index(drop=True)

    @staticmethod
    def as_lines(
        centrosomes_df: pd.DataFrame, track_colors: dict
    ) -> tuple[list, list, list]:
        lines = []
        colors = []
        centrosome_ids = []
        for centrosome_id, group in centrosomes_df.groupby("centrosome_id"):
            centriole_ids = group["centriole_id"].unique()
            components = []
            for centriole_id in centriole_ids:
                centriole_points = group[group["centriole_id"] == centriole_id][
                    ["T", "Y", "X"]
                ]
                centriole_points = centriole_points.sort_values(by="T")
                centriole_points = centriole_points.values
                components.append(centriole_points)
            points = np.stack(components, axis=1)
            lines.append([p for p in points])
            colors.append(
                [track_colors.get(centrosome_id, "#ffffff") for _ in range(len(points))]
            )
            centrosome_ids.append(centrosome_id)
        return centrosome_ids, lines, colors


def launch_full_process():
    from pathlib import Path
    import tifffile as tiff
    import napari

    DUMP = Path("/home/clement/Documents/projects/nucleation/draft/implementation/dump")

    folder_in = Path("/home/clement/Documents/projects/nucleation/3VPCs")
    filename = "251119_#4_30_001_016.vsi - C561.tif"
    path_in = folder_in / filename

    calib = {"T": 1, "Y": 0.1083333, "X": 0.1083333}
    units = {"T": "s", "Y": "um", "X": "um"}

    img = tiff.imread(path_in)  # (T, Y, X)
    hints = {
        1: {"start": 17, "end": 463, "points": np.array([[347, 90], [337, 132]])},
        2: {"start": 17, "end": 463, "points": np.array([[248, 477], [237, 514]])},
    }

    operator = FindCentrosomesOperator()
    operator.set_input_image(img, calib, units)
    operator.set_hints(hints)
    operator.set_prominence(10.0)
    operator.run()

    centrosomes = operator.get_centrosomes()
    centrosomes.to_csv(DUMP / "centrosome_points_track.csv", index=False)

    viewer = napari.Viewer()
    viewer.add_image(img, name="raw")
    viewer.add_points(
        centrosomes[["T", "Y", "X"]].values,
        name="maximas",
        size=7,
        face_color="transparent",
        border_color="red",
        opacity=0.75,
        visible=True,
    )
    viewer.add_tracks(
        centrosomes[["centriole_id", "T", "Y", "X"]].values,
        name="tracks",
        features=centrosomes,
    )
    napari.run()


if __name__ == "__main__":
    launch_full_process()
