import napari
import numpy as np
import pandas as pd
from pathlib import Path
import shutil
import json
import tifffile as tiff


def _stem_to_id(path: Path):
    try:
        return int(path.stem)
    except ValueError:
        return path.stem


class ImportArchiveOperator:

    def __init__(self, viewer: napari.Viewer, widget=None):
        self.viewer = viewer
        self.widget = widget
        self.root_path = None
        self.settings = {}
        self.centrosomes_tracks = None
        self.tracked_centrosomes = None
        self.centrosomes_lines = None
        self.arcs = None
        self.kymographs = None
        self.spots = None
        self.kymos_order = []
        self.triplets = []

    def _get_widget(self):
        if self.widget is not None:
            return self.widget
        from napari_nucleation_analyzer.widget import CentrosomesWidget
        for widget in self.viewer.window.dock_widgets.values():
            if type(widget) == CentrosomesWidget:
                return widget
        raise ValueError("CentrosomesWidget not found in the viewer.")

    def set_root_path(self, path: str | Path):
        self.root_path = Path(path)

    def _as_buffer_folder(self):
        if self.root_path is None:
            raise ValueError("Root path is not set. Please set the root path before running the operator.")
        buffer_name = "_" + self.root_path.stem
        buffer_path = self.root_path.parent / buffer_name
        if buffer_path.exists():
            shutil.rmtree(buffer_path)
        buffer_path.mkdir(parents=True, exist_ok=True)
        shutil.unpack_archive(str(self.root_path), str(buffer_path), "zip")
        return buffer_path

    def _restore_settings(self, folder: Path):
        with open(folder / "settings.json", "r") as f:
            return json.load(f)

    def _check_image_present(self):
        layer_name = self.settings["image-name"]
        if layer_name not in self.viewer.layers:
            raise ValueError(f"Image layer '{layer_name}' not found in the viewer.")

    def _restore_centrosomes_tracks(self, folder: Path):
        with open(folder / "centrosomes_tracks.json", "r") as f:
            raw = json.load(f)
        return {int(k): v for k, v in raw.items()}

    def _restore_tracked_centrosomes(self, folder: Path):
        return pd.read_csv(folder / "tracked_centrosomes.csv")

    def _restore_centrosomes_lines(self, folder: Path):
        lines_folder = folder / "centrosomes_lines"
        lines = {}
        for path in lines_folder.glob("*.npy"):
            lines[_stem_to_id(path)] = np.load(path)
        return lines

    def _restore_arcs(self, folder: Path):
        arcs_folder = folder / "arcs"
        arcs = {}
        for path in arcs_folder.glob("*.npy"):
            arcs[_stem_to_id(path)] = np.load(path)
        return arcs

    def _restore_kymographs(self, folder: Path):
        kymos_folder = folder / "kymographs"
        kymographs = {}
        for path in kymos_folder.glob("*.tif"):
            kymographs[_stem_to_id(path)] = tiff.imread(path)
        return kymographs

    def _restore_spots(self, folder: Path):
        spots_folder = folder / "spots"
        spots = {}
        for path in spots_folder.glob("*.npy"):
            spots[_stem_to_id(path)] = np.load(path)
        return spots

    def _inject_centrosomes_lines(self, widget):
        if self.centrosomes_tracks is None or self.centrosomes_lines is None:
            raise ValueError("Centrosomes tracks or lines not set. Please set them before injecting.")
        for c_id, properties in self.centrosomes_tracks.items():
            centrosome_id = int(c_id)
            start = properties["start"]
            end = properties["end"]
            line = self.centrosomes_lines.get(centrosome_id, None)
            if line is None:
                raise ValueError(f"No line found for centriole ID {centrosome_id}.")
            widget.tracks_manager_widget.insert_track(centrosome_id, line, start, end)

    def _inject_settings(self, widget):
        if self.settings is None:
            raise ValueError("Settings not set. Please set settings before injecting.")
        if 'prominence' in self.settings:
            widget.prominence_input.setValue(self.settings['prominence'])
        if 'searching-range' in self.settings:
            widget.searching_range_input.setValue(self.settings['searching-range'])
        if 'memory' in self.settings:
            widget.memory_input.setValue(self.settings['memory'])
        if 'max-binding-distance' in self.settings:
            widget.max_binding_distance_input.setValue(self.settings['max-binding-distance'])
        if 'arc-radius' in self.settings:
            widget.radius_input.setValue(self.settings['arc-radius'])
        if 'arc-angle' in self.settings:
            widget.angle_input.setValue(self.settings['arc-angle'])
        if 'spots-prominence' in self.settings:
            widget.spot_prominence_input.setValue(self.settings['spots-prominence'])

    # def _inject_tracked_centrosomes(self, widget):
    #     result = self.tracked_centrosomes
    #     centrosomes_tracks_layer_name = widget.centrioles_tracks_prefix + self.settings["image-name"]
    #     image_layer = self.viewer.layers[self.settings["image-name"]]

    #     self.viewer.add_tracks(
    #         result[['centriole_id', 'T', 'Y', 'X']],
    #         name=centrosomes_tracks_layer_name,
    #         scale=image_layer.scale,
    #         features=result,
    #         graph=None,
    #         tail_length=8,
    #         hide_completed_tracks=True,
    #         tail_width=3,
    #         units=image_layer.units
    #     )

    def _inject_tracked_centrosomes(self, widget):
        result = self.tracked_centrosomes
        centrosomes_tracks_layer_name = widget.centrioles_tracks_prefix + self.settings["image-name"]
        image_layer = self.viewer.layers[self.settings["image-name"]]

        if result is None or result.empty:
            raise ValueError("Tracked centrosomes data is empty. Please ensure that tracked centrosomes are set before injecting.")

        self.triplets.append((
            result[['centriole_id', 'T', 'Y', 'X']],
            {
                'name': centrosomes_tracks_layer_name,
                'scale': image_layer.scale,
                'features': result,
                'graph': None,
                'tail_length': 8,
                'hide_completed_tracks': True,
                'tail_width': 3,
                'units': image_layer.units
            },
            'tracks'
        ))

    def _rebuild_pairs(self):
        if self.tracked_centrosomes is None:
            raise ValueError("Tracked centrosomes not set. Please set tracked centrosomes before rebuilding pairs.")
        unique_centrioles = self.tracked_centrosomes['centriole_id'].unique()
        pairs = {}
        for centriole_id in unique_centrioles:
            centriole_data = self.tracked_centrosomes[self.tracked_centrosomes['centriole_id'] == centriole_id]
            unique_centrosomes = centriole_data['centrosome_id'].unique()
            centrosomes = unique_centrosomes.tolist()
            if len(centrosomes) != 1:
                raise ValueError(f"Centriole {centriole_id} has multiple centrosomes: {centrosomes}")
            pairs[centriole_id] = centrosomes[0]
        return pairs

    # def _inject_arcs(self, widget):
    #     if self.arcs is None:
    #         raise ValueError("Arcs not set. Please set arcs before injecting.")
    #     image_layer = self.viewer.layers[self.settings["image-name"]]

    #     for centriole_id, arc in self.arcs.items():
    #         centrosome_id = int(self.pairs[centriole_id])
    #         color = widget.tracks_manager_widget._rows[centrosome_id].color
    #         self.viewer.add_shapes(
    #             arc,
    #             shape_type='path',
    #             edge_color=color,
    #             edge_width=2,
    #             face_color='transparent',
    #             opacity=0.75,
    #             name=f"{widget.arcs_shapes_prefix}{centriole_id}",
    #             scale=image_layer.scale,
    #             units=image_layer.units
    #         )

    def _inject_arcs(self, widget):
        if self.arcs is None:
            raise ValueError("Arcs not set. Please set arcs before injecting.")
        image_layer = self.viewer.layers[self.settings["image-name"]]

        for centriole_id, arc in self.arcs.items():
            centrosome_id = int(self.pairs[centriole_id])
            color = widget.tracks_manager_widget._rows[centrosome_id].color
            self.triplets.append((
                arc,
                {
                    'shape_type': 'path',
                    'edge_color': color,
                    'edge_width': 2,
                    'face_color': 'transparent',
                    'opacity': 0.75,
                    'name': f"{widget.arcs_shapes_prefix}{centriole_id}",
                    'scale': image_layer.scale,
                    'units': image_layer.units
                },
                'shapes'
            ))

    # def _inject_kymographs(self, widget):
    #     if self.kymographs is None:
    #         raise ValueError("Kymographs not set. Please set kymographs before injecting.")

    #     centriole_ids = [(k, v) for k, v in self.pairs.items()]
    #     centriole_ids.sort(key=lambda x: x[1]) # sorted by centrosome_id
    #     centrosome_ids = [int(v) for _, v in centriole_ids]
    #     centriole_ids = [int(k) for k, _ in centriole_ids]

    #     kymographs = self.kymographs
    #     padding = 10

    #     # Hide all the layers except the original image
    #     for layer in self.viewer.layers:
    #         layer.visible = False

    #     # create a list of polygons
    #     polygons = []
    #     for i, centriole_id in enumerate(centriole_ids):
    #         kymograph = kymographs[centriole_id]
    #         T, Y = kymograph.shape
    #         polygon = np.array([
    #             [0, i * (Y + padding)],
    #             [T - 1, i * (Y + padding)],
    #             [T - 1, i * (Y + padding) + Y - 1],
    #             [0, i * (Y + padding) + Y - 1],
    #         ])
    #         polygons.append(polygon)

    #     # create features
    #     features = {
    #         'centriole_id': centriole_ids,
    #         'centrosome_id': centrosome_ids
    #     }

    #     text = {
    #         'string': 'C{centrosome_id} -> c{centriole_id}',
    #         'anchor': 'upper_left',
    #         'translation': [-5, 0],
    #         'size': 16,
    #         'color': 'white',
    #     }

    #     colors = [widget.tracks_manager_widget._rows[int(self.pairs[centriole_id])].color for centriole_id in centriole_ids]

    #     self.viewer.add_shapes(
    #         polygons,
    #         features=features,
    #         shape_type='polygon',
    #         edge_width=3,
    #         opacity=1.0,
    #         edge_color=colors,
    #         face_color='transparent',
    #         text=text,
    #         name='kymo_outlines'
    #     )

    #     # add the images
    #     for i, centriole_id in enumerate(centriole_ids):
    #         kymograph = kymographs[centriole_id]
    #         self.viewer.add_image(
    #             kymograph,
    #             name=f"{widget.kymo_prefix}{centriole_id}",
    #             scale=(1, 1),
    #             translate=(0, i * (kymograph.shape[1] + padding)),
    #             colormap='turbo'
    #         )

    def _inject_kymographs(self, widget):
        if self.kymographs is None:
            raise ValueError("Kymographs not set. Please set kymographs before injecting.")

        centriole_ids = [(k, v) for k, v in self.pairs.items()]
        centriole_ids.sort(key=lambda x: x[1]) # sorted by centrosome_id
        centrosome_ids = [int(v) for _, v in centriole_ids]
        centriole_ids = [int(k) for k, _ in centriole_ids]
        self.kymos_order = centriole_ids  # Store the order of kymographs for later use

        kymographs = self.kymographs
        padding = 10

        # create a list of polygons
        polygons = []
        for i, centriole_id in enumerate(centriole_ids):
            kymograph = kymographs[centriole_id]
            T, Y = kymograph.shape
            polygon = np.array([
                [0, i * (Y + padding)],
                [T - 1, i * (Y + padding)],
                [T - 1, i * (Y + padding) + Y - 1],
                [0, i * (Y + padding) + Y - 1],
            ])
            polygons.append(polygon)

        # create features
        features = {
            'centriole_id': centriole_ids,
            'centrosome_id': centrosome_ids
        }

        text = {
            'string': 'C{centrosome_id} -> c{centriole_id}',
            'anchor': 'upper_left',
            'translation': [-5, 0],
            'size': 16,
            'color': 'white',
        }

        colors = [widget.tracks_manager_widget._rows[int(self.pairs[centriole_id])].color for centriole_id in centriole_ids]
        self.triplets.append((
            polygons,
            {
                'features': features,
                'shape_type': 'polygon',
                'edge_width': 3,
                'opacity': 1.0,
                'edge_color': colors,
                'face_color': 'transparent',
                'text': text,
                'name': 'kymo_outlines'
            },
            'shapes'
        ))

        # add the images
        for i, centriole_id in enumerate(centriole_ids):
            kymograph = kymographs[centriole_id]
            self.triplets.append((
                kymograph,
                {
                    'name': f"{widget.kymo_prefix}{centriole_id}",
                    'scale': (1, 1),
                    'translate': (0, i * (kymograph.shape[1] + padding)),
                    'colormap': 'turbo'
                },
                'image'
            ))

    def _inject_spots(self, widget):
        if self.spots is None:
            raise ValueError("Spots not set. Please set spots before injecting.")

        coordinates = self.spots
        if len(self.kymos_order) != len(coordinates):
            raise ValueError("Mismatch between the number of kymographs and spots. Please ensure that both are set correctly.")

        padding = 10
        kymographs = self.kymographs

        for i, centriole_id in enumerate(self.kymos_order):
            coords = coordinates[centriole_id]
            spots_layer_name = f"{widget.spots_layer_prefix}{centriole_id}"
            kymograph = kymographs[centriole_id]

            self.triplets.append((
                coords,
                {
                    'name': spots_layer_name,
                    'scale': (1, 1),
                    'face_color': 'transparent',
                    'size': 5,
                    'translate': (0, i * (kymograph.shape[1] + padding))
                },
                'points'
            ))

    def _inject_values(self):
        widget = self._get_widget()
        self._inject_settings(widget)
        self._inject_tracked_centrosomes(widget)
        self._inject_centrosomes_lines(widget)
        self._inject_arcs(widget)
        self._inject_kymographs(widget)
        self._inject_spots(widget)

    def run(self):
        if self.root_path is None:
            raise ValueError("Root path is not set. Please set the root path before running the operator.")
        buffer_folder = self._as_buffer_folder()
        self.settings = self._restore_settings(buffer_folder)
        self._check_image_present()
        self.centrosomes_tracks = self._restore_centrosomes_tracks(buffer_folder)
        self.tracked_centrosomes = self._restore_tracked_centrosomes(buffer_folder)
        self.centrosomes_lines = self._restore_centrosomes_lines(buffer_folder)
        self.arcs = self._restore_arcs(buffer_folder)
        self.kymographs = self._restore_kymographs(buffer_folder)
        self.spots = self._restore_spots(buffer_folder)
        self.pairs = self._rebuild_pairs()
        shutil.rmtree(buffer_folder)
        self._inject_values()