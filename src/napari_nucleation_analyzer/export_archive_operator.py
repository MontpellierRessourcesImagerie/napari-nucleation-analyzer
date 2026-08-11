import napari
import numpy as np
from pathlib import Path
import shutil
import json
import tifffile as tiff

class ExportArchiveOperator:

    def __init__(self, widget, viewer: napari.Viewer):
        self.viewer = viewer
        self.widget = widget
        self.root_path = None

    def _gather_settings(self, folder: Path):
        settings = {
            'image-name': self.widget.image_combo.currentText(),
            'prominence': self.widget.prominence_input.value(),
            'searching-range': self.widget.searching_range_input.value(),
            'memory': self.widget.memory_input.value(),
            'max-binding-distance': self.widget.max_binding_distance_input.value(),

            'arc-radius': self.widget.radius_input.value(),
            'arc-angle': self.widget.angle_input.value(),

            'spots-prominence': self.widget.spot_prominence_input.value()
        }
        with open(folder / "settings.json", "w") as f:
            json.dump(settings, f)

    def _gather_centrosomes_tracks(self, folder: Path):
        tracks = {}
        o = self.widget.tracks_manager_widget
        for centrosome_id, row in o._rows.items():
            item = {
                'start': row.frame_start,
                'end'  : row.frame_end,
                'color': row.color,
                'index': row.index
            }
            tracks[centrosome_id] = item
        with open(folder / "centrosomes_tracks.json", "w") as f:
            json.dump(tracks, f)

    def _gather_centrosomes_lines(self, folder: Path):
        o = self.widget.tracks_manager_widget
        lines_folder = folder / "centrosomes_lines"
        lines_folder.mkdir(exist_ok=True)

        for centrosome_id, row in o._rows.items():
            layer_name = o.as_track_layer_name(centrosome_id)
            layer = self.viewer.layers[layer_name] if layer_name in self.viewer.layers else None

            if layer is None:
                continue

            line = layer.data
            np.save(lines_folder / f"{centrosome_id}.npy", line)

    def _gather_tracked_centrosomes(self, folder: Path):
        image_layer_name = self.widget.image_combo.currentText()
        tracked_centrosomes_layer_name = self.widget.centrioles_tracks_prefix + image_layer_name
        tracked_centrosomes_layer = (
            self.viewer.layers[tracked_centrosomes_layer_name] 
            if tracked_centrosomes_layer_name in self.viewer.layers 
            else None
        )
        if tracked_centrosomes_layer is None:
            raise ValueError(f"Tracked centrosomes layer '{tracked_centrosomes_layer_name}' not found in the viewer.")
        df = tracked_centrosomes_layer.features
        df.to_csv(folder / "tracked_centrosomes.csv", index=False)


    def set_root_path(self, path: str | Path):
        self.root_path = Path(path)

    def _as_buffer_folder(self):
        if self.root_path is None:
            raise ValueError("Root path is not set. Please set the root path before running the operator.")
        buffer_name = "_" + self.root_path.name
        buffer_path = self.root_path.parent / buffer_name
        if buffer_path.exists():
            shutil.rmtree(buffer_path)
        buffer_path.mkdir(parents=True, exist_ok=True)
        return buffer_path

    def _gather_arcs(self, folder: Path):
        arcs_folder = folder / "arcs"
        arcs_folder.mkdir(exist_ok=True)
        arcs = self.widget._gather_arcs()
        for centriole_id, arc in arcs.items():
            np.save(arcs_folder / f"{centriole_id}.npy", arc)

    def _gather_kymographs(self, folder: Path):
        kymos_folder = folder / "kymographs"
        kymos_folder.mkdir(exist_ok=True)
        kymographs = self.widget._gather_kymographs()
        for centriole_id, kymo in kymographs.items():
            tiff.imwrite(kymos_folder / f"{centriole_id}.tif", kymo)

    def _gather_spots(self, folder: Path):
        spots_folder = folder / "spots"
        spots_folder.mkdir(exist_ok=True)
        spots = self.widget._gather_spots()
        for centriole_id, spot in spots.items():
            np.save(spots_folder / f"{centriole_id}.npy", spot)

    def run(self):
        if self.root_path is None:
            raise ValueError("Root path is not set. Please set the root path before running the operator.")
        buffer_folder = self._as_buffer_folder()
        print(f">>> {buffer_folder}")
        self._gather_settings(buffer_folder)
        self._gather_centrosomes_tracks(buffer_folder)
        self._gather_tracked_centrosomes(buffer_folder)
        self._gather_centrosomes_lines(buffer_folder)
        self._gather_arcs(buffer_folder)
        self._gather_kymographs(buffer_folder)
        self._gather_spots(buffer_folder)
        shutil.make_archive(str(self.root_path), "zip", str(buffer_folder))
        shutil.rmtree(buffer_folder)
    