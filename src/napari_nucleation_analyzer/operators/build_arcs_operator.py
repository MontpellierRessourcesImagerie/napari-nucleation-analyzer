import numpy as np
import pandas as pd
from scipy.spatial import KDTree

class BuildArcsOperator:

    def __init__(self):
        self.angle_deg   = None
        self.radius_um   = None
        self.calibration = None
        self.hints       = None
        self.centrosomes = None
        self.arcs        = None
        self.time_start  = None
        self.num_points  = 50

    def set_angle(self, angle_deg):
        if angle_deg < 0 or angle_deg > 360:
            raise ValueError("Angle must be between 0 and 360 degrees.")
        self.angle_deg = angle_deg

    def set_radius(self, radius_um):
        if radius_um <= 0:
            raise ValueError("Radius must be a positive value.")
        self.radius_um = radius_um

    def set_calibration(self, calibration):
        if len(calibration) != 3:
            raise ValueError("Calibration must be a tuple of three values (T, Y, X).")
        self.calibration = calibration

    def set_hints(self, hints):
        if not isinstance(hints, np.ndarray):
            raise ValueError("Hints must be a numpy array.")
        if hints.ndim != 2 or hints.shape[1] not in [2, 3]:
            raise ValueError("Hints must be a 2D array with shape (N, 2) or (N, 3) for (T, Y, X).")
        if hints.shape[1] == 2:
            self.hints = hints
        else:
            self.hints = hints[:, 1:]

    def set_centrosomes(self, centrosomes):
        self.centrosomes = centrosomes.copy()

    def set_num_points(self, num_points):
        if num_points <= 3:
            raise ValueError("Number of points must be greater than 3.")
        self.num_points = num_points

    def get_radius_pixels(self):
        if self.radius_um is None or self.calibration is None:
            raise ValueError("Both radius and calibration must be set before computing radius in pixels.")
        return self.radius_um / self.calibration[-1]
    
    def get_arcs(self):
        if self.arcs is None:
            raise ValueError("Arcs have not been built yet. Please run the operator first.")
        return self.arcs
    
    def set_time_start(self, time_start):
        if not isinstance(time_start, int) or time_start < 0:
            raise ValueError("Time start must be a non-negative integer.")
        self.time_start = time_start
    
    def _build_pairs(self, starters, tracked_points, thr=5.0):
        """
        The goal of this function is to use the user-defined starters to rebuild pairs of points
        from the tracked points. It uses a KDTree to find the closest tracked point to each starter point.

        Args:
            starters: list of pairs of points (Y, X) in the first frame (T=0) defined by the user.
            tracked_points: DataFrame containing tracked points with columns ["Y", "X", "T", "track_id"].

        Returns:
            A numpy array containing pairs of Track IDs.
        """
        points_t0 = tracked_points[tracked_points["T"] == self.time_start]
        tree = KDTree(starters)
        pairs = [-1 for _ in range(len(starters))]
        
        for r in points_t0.itertuples(index=False):
            dist, idx = tree.query([r.Y, r.X])
            if dist > thr:
                print(f"  No close starter point found (Distance={dist:.2f})")
                continue
            print(f"  Closest starter point: (Y={starters[idx][0]}, X={starters[idx][1]}), Distance={dist:.2f}")
            pairs[idx] = r.track_id
        
        return np.array(pairs)
    
    def _build_vectors(self, pairs, tracked_points):
        """
        This function modifies the tracked_points DataFrame to include a new VX and VY columns.
        These are the two components of a vector pointing towards the other elements of the pair.
        These vectors are normalized to have a unit length.
        If a pair is incomplete, the corresponding VX and VY values are set to NaN.

        Args:
            pairs: list of Track IDs corresponding to the starters.
            tracked_points: DataFrame containing tracked points with columns ["Y", "X", "T", "track_id"].
        """
        tracked_points["VX"] = np.nan
        tracked_points["VY"] = np.nan
        for i in range(0, len(pairs), 2):
            id1 = pairs[i]
            id2 = pairs[i + 1]
            if id1 == -1 or id2 == -1:
                print(f"Pair {i//2} is incomplete (IDs: {id1}, {id2})")
                continue
            point1 = tracked_points[tracked_points["track_id"] == id1]
            point2 = tracked_points[tracked_points["track_id"] == id2]
            if point1.empty or point2.empty:
                print(f"Pair {i//2} has missing points in tracked data (IDs: {id1}, {id2})")
                continue
            for t in range(max(point1["T"].min(), point2["T"].min()), min(point1["T"].max(), point2["T"].max()) + 1):
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
                    tracked_points.loc[(tracked_points["track_id"] == id1) & (tracked_points["T"] == t), "VX"] = vx
                    tracked_points.loc[(tracked_points["track_id"] == id1) & (tracked_points["T"] == t), "VY"] = vy
                    tracked_points.loc[(tracked_points["track_id"] == id2) & (tracked_points["T"] == t), "VX"] = -vx
                    tracked_points.loc[(tracked_points["track_id"] == id2) & (tracked_points["T"] == t), "VY"] = -vy

    def _build_arcs(self, tracked_points, radius_pxl, angle_degrees, num_points):
        angle_rad = np.radians(angle_degrees)

        local_angles = np.linspace(-angle_rad / 2, angle_rad / 2, num_points)
        ref_points = np.stack((np.sin(local_angles), np.cos(local_angles)), axis=1) * radius_pxl  # (num_points, 2), format (Y, X)

        arcs = {}
        track_ids = tracked_points["track_id"].unique()

        for track_id in track_ids:
            track_data = tracked_points[tracked_points["track_id"] == track_id].sort_values("T")

            vectors = track_data[["VY", "VX"]].values
            origins = track_data[["Y", "X"]].values
            theta   = np.arctan2(vectors[:, 0], vectors[:, 1])
            cos_t   = np.cos(theta)
            sin_t   = np.sin(theta)
            
            rotation_matrices = np.stack([
                np.stack([cos_t,  sin_t], axis=-1),
                np.stack([-sin_t, cos_t], axis=-1)
            ], axis=1)  # (T, 2, 2)

            rotated = np.einsum('tij,pj->tpi', rotation_matrices, ref_points)
            arcs[track_id] = rotated + origins[:, np.newaxis, :]

        return arcs
    
    def run(self):
        if self.angle_deg is None or self.radius_um is None or self.calibration is None:
            raise ValueError("Angle, radius, and calibration must be set before running the operator.")
        
        if self.hints is None or self.centrosomes is None:
            raise ValueError("Both hints and centrosomes must be set before running the operator.")
        
        pairs = self._build_pairs(self.hints, self.centrosomes)
        self._build_vectors(pairs, self.centrosomes)
        self.arcs = self._build_arcs(
            self.centrosomes, 
            self.get_radius_pixels(), 
            self.angle_deg, 
            self.num_points
        )

    @staticmethod
    def as_napari_shapes(arcs_dict, time_start):
        """
        Converts a dictionary indexed by track ID and containing arcs into a list of numpy arrays
        suitable for visualization in Napari.
        Each arc in the dict has the shape (T, num_points, 2).
        Napari expects the data's root to be a Python list.
        Each element inside it must be a NumPy array of shape (num_points, 3)
        """
        napari_shapes = []
        track_ids = []
        for track_id, arc in arcs_dict.items():
            T, num_points, _ = arc.shape
            for t in range(T):
                track_ids.append(track_id)
                points_3d = np.zeros((num_points, 3))
                points_3d[:, 0] = t + time_start
                points_3d[:, 1:] = arc[t]
                napari_shapes.append(points_3d)
        return napari_shapes, pd.DataFrame({'track_id': track_ids})


if __name__ == "__main__":
    import pandas as pd
    from pathlib import Path

    folder = Path("/home/clement/Documents/projects/nucleation/draft/implementation/dump")
    df_name = "centrosome_points_track.csv"
    df_path = folder / df_name

    centrosomes = pd.read_csv(df_path)

    hints = np.array([
        [17, 347, 90],
        [17, 337, 132],
        [17, 248, 477],
        [17, 237, 514]
    ])

    op = BuildArcsOperator()
    op.set_angle(90.0)
    op.set_radius(2.0)
    op.set_calibration((1, 0.1083333, 0.1083333))
    op.set_hints(hints)
    op.set_centrosomes(centrosomes)
    op.set_time_start(17)
    op.run()

    arcs = op.get_arcs()
    napari_shapes, track_ids = BuildArcsOperator.as_napari_shapes(arcs, time_start=17)
