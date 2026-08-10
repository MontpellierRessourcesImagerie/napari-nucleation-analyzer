import numpy as np
import pandas as pd
import xarray as xr


class BuildArcsOperator:
    """
    Arcs are stored in a dictionary. The keys are the centrioles IDs as they appear in the DataFrame.
    For a given centriole, the arc over time is represented by a unique np.array of shape (duration, num_points, 3).
    'duration' is the number of time points for which the centriole has been tracked.
    'num_points' is the number of points processed to form the arc at each time point.
    The last dimension of size 3 corresponds to the coordinates (T, Y, X) of each point in the arc.

    The original DataFrame is not modified by this operator, but the copy held internally contains extra information after execution.
    - "VX": The X component of the unit vector pointing towards the other centriole in the pair.
    - "VY": The Y component of the unit vector pointing towards the other centriole in the pair.
    - "distance": The distance between the two centrioles of the pair at each time point.

    The get_arcs method returns a copy of the arcs with the last dimension reordered to match the order of dimensions in the input image.
    """

    def __init__(self):
        self.angle_deg = None
        self.radius_um = None
        self.input_image = None
        self.centrosomes = None
        self.arcs = None
        self.pairs = None

    def set_angle(self, angle_deg):
        if angle_deg < 0 or angle_deg > 360:
            raise ValueError("Angle must be between 0 and 360 degrees.")
        self.angle_deg = angle_deg

    def set_radius(self, radius_um):
        if radius_um <= 0:
            raise ValueError("Radius must be a positive value.")
        self.radius_um = radius_um

    def set_input_image(self, image_arr, calibration, units):
        if image_arr.ndim != 3:
            raise ValueError("Input image must be a 3D array (T, Y, X).")
        self.input_image = xr.DataArray(
            image_arr,
            dims=["T", "Y", "X"],
            attrs={"scale": calibration, "units": units},
        )

    def set_centrosomes(self, centrosomes):
        required_columns = {"Y", "X", "T", "centriole_id", "centrosome_id"}
        if not required_columns.issubset(centrosomes.columns):
            raise ValueError(
                f"Centrosomes DataFrame must contain the following columns: {required_columns}"
            )
        self.centrosomes = centrosomes.copy()

    def get_radius_pixels(self):
        if self.radius_um is None or self.input_image is None:
            raise ValueError(
                "Both radius and input image must be set before computing radius in pixels."
            )
        return self.radius_um / self.input_image.attrs["scale"]["X"]

    def get_arcs(self):
        """
        By default, arcs have a shape of (duration, num_points, 3) where the last dimension corresponds to (T, Y, X).
        However, the dimensions of the input image are not necessarily in the same order.
        This function reorders a copy of the arcs according to what is present in input_image.dims.
        """
        if self.arcs is None or self.input_image is None:
            raise ValueError(
                "Arcs have not been built yet or input image is not set. Please run the operator first."
            )
        arcs_copy = {}
        axes = [str(ax) for ax in self.input_image.dims]
        dim_order = [axes.index(d) for d in ["T", "Y", "X"]]
        for centrosome_id, arc_data in self.arcs.items():
            arcs_copy[centrosome_id] = arc_data[..., dim_order]
        return arcs_copy

    def get_centrosomes(self) -> pd.DataFrame:
        if self.centrosomes is None:
            raise ValueError(
                "Centrosomes have not been set or processed yet. Please run the operator first."
            )
        return self.centrosomes.copy()

    def get_pairs(self, flipped=False):
        """
        Returns a dictionary where the keys are centrosome IDs and the values are tuples of centriole IDs.
        If 'flipped' is True, the dictionary associates each centriole ID to its corresponding centrosome ID.
        """
        if self.pairs is None:
            raise ValueError(
                "Pairs have not been built yet. Please run the operator first."
            )
        if flipped:
            flipped_pairs = {}
            for k, v in self.pairs.items():
                flipped_pairs[v[0]] = k
                flipped_pairs[v[1]] = k
            return flipped_pairs
        else:
            return self.pairs.copy()

    def process_n_points(self, radius_pxl, angle_degrees):
        """
        Computes how many points are needed to sample an arc of a given angle and radius.
        """
        arc_length = 2.1 * np.pi * radius_pxl * (angle_degrees / 360.0)
        num_points = max(2, int(np.ceil(arc_length)))  # At least 2 points
        return num_points

    def _check_valid_pairs(self, centrosomes):
        """
        Returns a dict of valid centrosomes with pairs of valid centrioles.
        Validity criteria:
        - For a centrosome ID, there must be exactly two centrioles IDs.
        - The time points for both centrioles must overlap exactly in time.
        - A given centriole ID must not be part of several centrosomes.
        """
        centrosome_groups = centrosomes.groupby("centrosome_id")
        valid_pairs = {}
        used_centriole_ids = set()
        for centrosome_id, group in centrosome_groups:
            centriole_ids = group["centriole_id"].unique()
            if len(centriole_ids) != 2:
                print(
                    f"Centrosome {centrosome_id} does not have exactly two centrioles. Skipping."
                )
                continue
            id1, id2 = centriole_ids
            t1 = set(group[group["centriole_id"] == id1]["T"])
            t2 = set(group[group["centriole_id"] == id2]["T"])
            if t1 != t2:
                print(
                    f"Centrosome {centrosome_id} has centrioles with non-overlapping time points. Skipping."
                )
                continue
            if id1 in used_centriole_ids or id2 in used_centriole_ids:
                print(
                    f"Centrioles {id1} or {id2} have already been used in another centrosome. Skipping."
                )
                continue
            valid_pairs[centrosome_id] = (id1, id2)
            used_centriole_ids.update([id1, id2])
        return valid_pairs

    def _build_vectors(self, pairs, tracked_points):
        """
        This function modifies the tracked_points DataFrame to include a new VX and VY columns.
        These are the two components of a vector pointing towards the other elements of the pair.
        These vectors are normalized to have a unit length.
        If a pair is incomplete, the corresponding VX and VY values are set to NaN.

        Args:
            pairs: list of Track IDs corresponding to the starters.
            tracked_points: DataFrame containing tracked points with columns ["Y", "X", "T", "centriole_id"].
        """
        tracked_points = tracked_points.copy()
        tracked_points["VX"] = np.nan
        tracked_points["VY"] = np.nan
        for centrosome_id, (id1, id2) in pairs.items():
            point1 = tracked_points[tracked_points["centriole_id"] == id1]
            point2 = tracked_points[tracked_points["centriole_id"] == id2]
            if point1.empty or point2.empty:
                print(
                    f"Centrosome {centrosome_id} has no points in tracked data (IDs: {id1}, {id2})"
                )
                continue
            for t in range(point1["T"].min(), point1["T"].max() + 1):
                p1 = point1[point1["T"] == t]
                p2 = point2[point2["T"] == t]
                if p1.empty or p2.empty:
                    continue
                dy = p2.iloc[0]["Y"] - p1.iloc[0]["Y"]
                dx = p2.iloc[0]["X"] - p1.iloc[0]["X"]
                norm = np.sqrt(dy**2 + dx**2)
                if norm > 0:
                    vx = dx / norm
                    vy = dy / norm
                    tracked_points.loc[
                        (tracked_points["centriole_id"] == id1)
                        & (tracked_points["T"] == t),
                        "VX",
                    ] = vx
                    tracked_points.loc[
                        (tracked_points["centriole_id"] == id1)
                        & (tracked_points["T"] == t),
                        "VY",
                    ] = vy
                    tracked_points.loc[
                        (tracked_points["centriole_id"] == id2)
                        & (tracked_points["T"] == t),
                        "VX",
                    ] = -vx
                    tracked_points.loc[
                        (tracked_points["centriole_id"] == id2)
                        & (tracked_points["T"] == t),
                        "VY",
                    ] = -vy
        return tracked_points

    def _build_arcs(self, tracked_points, radius_pxl, angle_degrees):
        num_points = self.process_n_points(radius_pxl, angle_degrees)
        angle_rad = np.radians(angle_degrees)

        local_angles = np.linspace(-angle_rad / 2, angle_rad / 2, num_points)
        ref_points = (
            np.stack((np.sin(local_angles), np.cos(local_angles)), axis=1) * radius_pxl
        )  # (num_points, 2), format (Y, X)

        arcs = {}
        centriole_ids = tracked_points["centriole_id"].unique()

        for centriole_id in centriole_ids:
            track_data = tracked_points[
                tracked_points["centriole_id"] == centriole_id
            ].sort_values("T")
            t_start = track_data["T"].min()
            t_end = track_data["T"].max()

            vectors = track_data[["VY", "VX"]].values
            origins = track_data[["Y", "X"]].values
            theta = np.arctan2(vectors[:, 0], vectors[:, 1])
            cos_t = np.cos(theta)
            sin_t = np.sin(theta)

            rotation_matrices = np.stack(
                [np.stack([cos_t, sin_t], axis=-1), np.stack([-sin_t, cos_t], axis=-1)],
                axis=1,
            )  # (T, 2, 2)

            rotated = np.einsum("tij,pj->tpi", rotation_matrices, ref_points)
            arcs2d = rotated + origins[:, np.newaxis, :]

            time_range = np.arange(t_start, t_end + 1)
            time_stack = np.repeat(time_range[:, np.newaxis], num_points, axis=1)
            time_stack = time_stack[..., np.newaxis]
            arcs3d = np.concatenate([time_stack, arcs2d], axis=-1)  # (T, num_points, 3)
            arcs[centriole_id] = arcs3d

        return arcs

    def _process_distance(self, df, pairs):
        """
        Computes the distance between the two centrioles of each pair at each time point.
        Returns a DataFrame with columns ["T", "centriole_id", "distance"].
        """
        df = df.set_index(["centriole_id", "T"]).sort_index()

        for _, (c1, c2) in pairs.items():
            sub1 = df.loc[c1]
            sub2 = df.loc[c2]

            dist = np.sqrt((sub1["X"] - sub2["X"]) ** 2 + (sub1["Y"] - sub2["Y"]) ** 2)

            df.loc[c1, "distance"] = dist.values
            df.loc[c2, "distance"] = dist.values

        df = df.reset_index()
        return df

    def run(self):
        if self.angle_deg is None or self.radius_um is None or self.input_image is None:
            raise ValueError(
                "Angle, radius, and image must be set before running the operator."
            )

        if self.centrosomes is None:
            raise ValueError("Centrosomes must be set before running the operator.")

        self.pairs = self._check_valid_pairs(self.centrosomes)
        self.centrosomes = self._build_vectors(self.pairs, self.centrosomes)
        self.centrosomes = self._process_distance(self.centrosomes, self.pairs)

        self.arcs = self._build_arcs(
            self.centrosomes, self.get_radius_pixels(), self.angle_deg
        )


if __name__ == "__main__":
    from pathlib import Path
    import tifffile as tiff
    import napari

    # Importing the image
    folder_img = Path("/home/clement/Documents/projects/nucleation/3VPCs")
    filename = "251119_#4_30_001_016.vsi - C561.tif"
    path_in = folder_img / filename

    calib = {"T": 1, "Y": 0.1083333, "X": 0.1083333}
    units = {"T": "s", "Y": "um", "X": "um"}

    img = tiff.imread(path_in)  # (T, Y, X)

    # Centrosomes DataFrame
    folder_df = Path(
        "/home/clement/Documents/projects/nucleation/draft/implementation/dump"
    )
    df_name = "centrosome_points_track.csv"
    df_path = folder_df / df_name

    centrosomes = pd.read_csv(df_path)

    # Running the operator
    op = BuildArcsOperator()
    op.set_angle(90.0)
    op.set_radius(2.0)
    op.set_input_image(img, calib, units)
    op.set_centrosomes(centrosomes)
    op.run()

    # Exporting the updated centrosomes DataFrame with additional columns
    updated_df_path = folder_df / "centrosomes_arcs.csv"
    op.get_centrosomes().to_csv(updated_df_path, index=False)

    # Showing controls in Napari
    viewer = napari.Viewer()

    image_layer = viewer.add_image(img, name="2D+t image")

    arcs = op.get_arcs()
    for centriole_id, arc in arcs.items():
        viewer.add_shapes(
            arc,
            shape_type='path',
            edge_color='red',
            edge_width=1,
            face_color='transparent',
            opacity=0.75,
            name=f"Arc {centriole_id}"
        )

    napari.run()
