from qtpy.QtWidgets import (
    QWidget,
    QVBoxLayout,
    QHBoxLayout,
    QGroupBox,
    QSpinBox,
    QDoubleSpinBox,
    QPushButton,
    QComboBox,
    QLabel,
    QFileDialog
)

from qtpy.QtGui import QFont
import napari
from pathlib import Path

from napari.qt.threading import create_worker
from napari.utils.notifications import show_info, show_warning

from napari.layers.labels.labels import Labels
from napari.layers.points.points import Points
from napari.layers.shapes.shapes import Shapes
from napari.layers.image.image import Image
from napari.layers.tracks.tracks import Tracks

import numpy as np
import pandas as pd

from .operators import (
    FindCentrosomesOperator,
    BuildArcsOperator,
    MakeKymographOperator,
    FindSpotsOperator,
    SummaryOperator
)
from .tracks_manager_widget import TracksManagerWidget
from .export_archive_operator import ExportArchiveOperator

class CentrosomesWidget(QWidget):

    centrioles_tracks_prefix = "_Tracked-centrioles-"
    arcs_shapes_prefix = "_Arc "
    kymo_prefix = "_Kymo "
    spots_layer_prefix = "_Spots "

    def __init__(self, viewer: "napari.viewer.Viewer"): # type: ignore
        super().__init__()
        self.viewer = viewer
        self.image_comboboxes = []
        self.label_comboboxes = []
        self.shape_comboboxes = []
        self.points_comboboxes = []
        self.track_comboboxes = []
        self.custom_font = QFont()
        self.custom_font.setFamily("Arial Unicode MS, Segoe UI Emoji, Apple Color Emoji, Noto Color Emoji")
        self.init_ui()
        self.current_operator = None
        self.init_callbacks()

    # -------- UI: ----------------------------------

    def init_ui(self):
        layout = QVBoxLayout()
        self.set_tracks_panel(layout)
        self.find_centrosomes_panel(layout)
        self.build_arcs_panel(layout)
        self.kymographs_panel(layout)
        self.setLayout(layout)

    def init_callbacks(self):
        self.viewer.layers.events.inserted.connect(self.update_comboboxes)
        self.viewer.layers.events.removed.connect(self.update_comboboxes)
        self.viewer.layers.events.renamed.connect(self.update_comboboxes)
        self.update_comboboxes()

    def _update_comboboxes(self, comboboxes, layer_type):
        for combobox in comboboxes:
            current_selection = combobox.currentText()
            combobox.clear()
            for layer in self.viewer.layers:
                if layer.name.startswith("_"):
                    continue
                if isinstance(layer, layer_type):
                    combobox.addItem(layer.name)
            index = combobox.findText(current_selection)
            if index >= 0:
                combobox.setCurrentIndex(index)

    def update_comboboxes(self):
        self._update_comboboxes(self.image_comboboxes, Image)
        self._update_comboboxes(self.label_comboboxes, Labels)
        self._update_comboboxes(self.points_comboboxes, Points)
        self._update_comboboxes(self.track_comboboxes, Tracks)
        self._update_comboboxes(self.shape_comboboxes, Shapes)

    def set_tracks_panel(self, parent_layout):
        self.tracks_group = QGroupBox("Centrosome tracks")
        layout = QVBoxLayout()

        self.tracks_manager_widget = TracksManagerWidget(self.viewer, self)
        layout.addWidget(self.tracks_manager_widget)

        self.tracks_group.setLayout(layout)
        parent_layout.addWidget(self.tracks_group)

    def find_centrosomes_panel(self, parent_layout):
        self.find_centrosomes_group = QGroupBox("Find centrosomes")
        layout = QVBoxLayout()

        # Image
        h_layout = QHBoxLayout()
        h_layout.addWidget(QLabel("Image"))
        self.image_combo = QComboBox()
        self.image_comboboxes.append(self.image_combo)
        h_layout.addWidget(self.image_combo)
        layout.addLayout(h_layout)

        # Prominence threshold
        h_layout = QHBoxLayout()
        h_layout.addWidget(QLabel("Prominence"))
        self.prominence_input = QDoubleSpinBox()
        self.prominence_input.setRange(0.0, 999.0)
        self.prominence_input.setSingleStep(0.5)
        self.prominence_input.setValue(FindCentrosomesOperator.default_prominence())
        h_layout.addWidget(self.prominence_input)
        layout.addLayout(h_layout)

        # Searching range
        h_layout = QHBoxLayout()
        h_layout.addWidget(QLabel("Searching range"))
        self.searching_range_input = QDoubleSpinBox()
        self.searching_range_input.setRange(0.0, 100.0)
        self.searching_range_input.setValue(FindCentrosomesOperator.default_searching_range())
        h_layout.addWidget(self.searching_range_input)
        layout.addLayout(h_layout)

        # Memory
        h_layout = QHBoxLayout()
        h_layout.addWidget(QLabel("Memory"))
        self.memory_input = QSpinBox()
        self.memory_input.setRange(0, 100)
        self.memory_input.setValue(FindCentrosomesOperator.default_memory())
        h_layout.addWidget(self.memory_input)
        layout.addLayout(h_layout)

        # Max binding distance
        h_layout = QHBoxLayout()
        h_layout.addWidget(QLabel("Max binding distance"))
        self.max_binding_distance_input = QDoubleSpinBox()
        self.max_binding_distance_input.setRange(0.0, 100.0)
        self.max_binding_distance_input.setValue(FindCentrosomesOperator.default_max_binding_distance())
        h_layout.addWidget(self.max_binding_distance_input)
        layout.addLayout(h_layout)

        # Vertical spacing
        layout.addSpacing(10)

        # Find centrosomes
        self.find_centrosomes_button = QPushButton("Find centrosomes")
        self.find_centrosomes_button.setFont(self.custom_font)
        self.find_centrosomes_button.clicked.connect(self.launch_find_centrosomes)
        layout.addWidget(self.find_centrosomes_button)

        self.find_centrosomes_group.setLayout(layout)
        parent_layout.addWidget(self.find_centrosomes_group)

    def get_image_calibration(self):
        image_layer_name = self.image_combo.currentText()
        if image_layer_name not in self.viewer.layers:
            raise ValueError("Selected image layer not found.")
        image_layer = self.viewer.layers[image_layer_name]
        s, u = image_layer.scale, image_layer.units
        a = image_layer.axis_labels
        return a, s, u

    def build_arcs_panel(self, parent_layout):
        self.build_arcs_group = QGroupBox("Build arcs")
        layout = QVBoxLayout()

        # Radius (µm)
        h_layout = QHBoxLayout()
        h_layout.addWidget(QLabel("Radius (µm)"))
        self.radius_input = QDoubleSpinBox()
        self.radius_input.setRange(0.0, 100.0)
        self.radius_input.setValue(2.0)
        self.radius_input.setSingleStep(0.1)
        h_layout.addWidget(self.radius_input)
        layout.addLayout(h_layout)

        # Angle (°)
        h_layout = QHBoxLayout()
        h_layout.addWidget(QLabel("Angle (°)"))
        self.angle_input = QDoubleSpinBox()
        self.angle_input.setRange(0.0, 360.0)
        self.angle_input.setValue(90.0)
        self.angle_input.setSingleStep(5.0)
        h_layout.addWidget(self.angle_input)
        layout.addLayout(h_layout)

        # Build arcs
        self.build_arcs_button = QPushButton("Build arcs")
        self.build_arcs_button.setFont(self.custom_font)
        self.build_arcs_button.clicked.connect(self.launch_build_arcs)
        layout.addWidget(self.build_arcs_button)

        self.build_arcs_group.setLayout(layout)
        parent_layout.addWidget(self.build_arcs_group)

    def kymographs_panel(self, parent_layout):
        self.kymographs_group = QGroupBox("Kymographs")
        layout = QVBoxLayout()

        # Build kymographs
        self.build_kymographs_button = QPushButton("Build kymographs")
        self.build_kymographs_button.setFont(self.custom_font)
        self.build_kymographs_button.clicked.connect(self.launch_build_kymographs)
        layout.addWidget(self.build_kymographs_button)

        # Locate spots
        h_layout = QHBoxLayout()

        self.spot_prominence_input = QDoubleSpinBox()
        self.spot_prominence_input.setRange(0.0, 10.0)
        self.spot_prominence_input.setSingleStep(0.05)
        self.spot_prominence_input.setValue(FindSpotsOperator.default_std_coeff())
        h_layout.addWidget(QLabel("Spot prominence"))
        h_layout.addWidget(self.spot_prominence_input)

        self.locate_spots_button = QPushButton("Locate spots")
        self.locate_spots_button.setFont(self.custom_font)
        self.locate_spots_button.clicked.connect(self.launch_locate_spots)
        h_layout.addWidget(self.locate_spots_button)

        layout.addLayout(h_layout)

        # Export

        h_layout = QHBoxLayout()

        self.export_summary_button = QPushButton("Export summary")
        self.export_summary_button.setFont(self.custom_font)
        self.export_summary_button.clicked.connect(self.launch_export_summary)
        h_layout.addWidget(self.export_summary_button)

        self.export_archive_button = QPushButton("Export archive")
        self.export_archive_button.setFont(self.custom_font)
        self.export_archive_button.clicked.connect(self.launch_export_archive)
        h_layout.addWidget(self.export_archive_button)

        layout.addLayout(h_layout)

        self.kymographs_group.setLayout(layout)
        parent_layout.addWidget(self.kymographs_group)

    # -------- Callbacks: ----------------------------------

    def set_enabled(self, enabled: bool):
        self.tracks_group.setEnabled(enabled)
        self.find_centrosomes_group.setEnabled(enabled)
        self.build_arcs_group.setEnabled(enabled)
        self.kymographs_group.setEnabled(enabled)

    # -------- Utils: --------------------------------------

    def launch_find_centrosomes(self):
        image_layer_name = self.image_combo.currentText()

        if image_layer_name not in self.viewer.layers:
            show_warning("Selected image layer not found.")
            return
        
        layer = self.viewer.layers[image_layer_name]
        hints = (
            self.tracks_manager_widget.as_hints()
            if 'hints' not in layer.metadata
            else layer.metadata['hints']
        )
        
        if hints is None:
            show_warning("No hints found.")
            return

        layer.metadata['hints'] = hints
        
        self.set_enabled(False)

        op = FindCentrosomesOperator()

        image = layer.data
        calib = {a: s for a, s in zip(layer.axis_labels, layer.scale)}
        units = {a: u for a, u in zip(layer.axis_labels, layer.units)}

        if set(calib.keys()) != {'T', 'Y', 'X'}:
            show_warning("Image must have axes T, Y, X.")
            self.set_enabled(True)
            return

        prominence = self.prominence_input.value()
        searching_range = self.searching_range_input.value()
        memory = self.memory_input.value()
        max_binding_distance = self.max_binding_distance_input.value()

        op.set_input_image(image, calib, units)
        op.set_hints(hints)
        op.set_prominence(prominence)
        op.set_searching_range(searching_range)
        op.set_memory(memory)
        op.set_max_binding_distance(max_binding_distance)
        
        self.current_operator = op

        def runner():
            try:
                return self.current_operator.run() if self.current_operator is not None else None
            except Exception as e:
                show_warning(f"Error during centrosome finding: {str(e)}")
                self.set_enabled(True)
                self.current_operator = None

        worker = create_worker(
            runner,
            _progress={
                "desc": "Building centrosomes..."
            },
        )
        worker.finished.connect(self.finished_find_centrosomes)
        worker.start()

    def finished_find_centrosomes(self, *args):
        if self.current_operator is None:
            return

        self.set_enabled(True)

        if not isinstance(self.current_operator, FindCentrosomesOperator):
            raise ValueError("Current operator is not a FindCentrosomesOperator.")
        
        image_layer_name = self.image_combo.currentText()
        centrosomes_tracks_layer_name = self.centrioles_tracks_prefix + image_layer_name
        image_layer = self.viewer.layers[image_layer_name]

        result = self.current_operator.get_centrosomes()
        required_cols = ["centriole_id", "T", "Y", "X", "centrosome_id"]
        if not all(col in result.columns for col in required_cols):
            raise ValueError(f"Resulting centrosomes dataframe must contain columns: {required_cols}")
        other_cols = [c for c in result.columns if c not in required_cols]
        result = result[required_cols + other_cols]

        # Adding/updating the tracks layer
        if centrosomes_tracks_layer_name in self.viewer.layers:
            layer = self.viewer.layers[centrosomes_tracks_layer_name]
            layer.data = result[['centriole_id', 'T', 'Y', 'X']]
            layer.features = result
        else:
            self.viewer.add_tracks(
                result[['centriole_id', 'T', 'Y', 'X']],
                name=centrosomes_tracks_layer_name,
                scale=image_layer.scale,
                features=result,
                graph=None,
                tail_length=8,
                hide_completed_tracks=True,
                tail_width=3,
                units=image_layer.units
            )

        track_colors = {k: v.color for k, v in self.tracks_manager_widget._rows.items()}
        centrosome_ids, lines, colors = FindCentrosomesOperator.as_lines(result, track_colors)

        for c_id, line, color in zip(centrosome_ids, lines, colors):
            layer_name = self.tracks_manager_widget.as_track_layer_name(c_id)
            if layer_name in self.viewer.layers:
                layer = self.viewer.layers[layer_name]
                self.viewer.layers.remove(layer)
            self.viewer.add_shapes(
                line,
                shape_type="line",
                edge_color=color,
                name=layer_name,
                opacity=0.75,
                ndim=3,
                scale=image_layer.scale,
                units=image_layer.units,
                edge_width=1
            )

        self.current_operator = None

    def launch_build_arcs(self):
        tracked_centrosomes_layer_name = self.centrioles_tracks_prefix + self.image_combo.currentText()
        if tracked_centrosomes_layer_name not in self.viewer.layers:
            show_warning("No tracked centrosomes layer found. Please run 'Find centrosomes' first.")
            return
        tracked_centrosomes_layer = self.viewer.layers[tracked_centrosomes_layer_name]
        tracked_centrosomes = tracked_centrosomes_layer.features

        image_layer_name = self.image_combo.currentText()
        layer = self.viewer.layers[image_layer_name]
        image = layer.data
        calib = {a: s for a, s in zip(layer.axis_labels, layer.scale)}
        units = {a: u for a, u in zip(layer.axis_labels, layer.units)}

        angle = self.angle_input.value()
        radius = self.radius_input.value()

        self.set_enabled(False)
        op = BuildArcsOperator()
        op.set_input_image(image, calib, units)
        op.set_angle(angle)
        op.set_radius(radius)
        op.set_centrosomes(tracked_centrosomes)

        self.current_operator = op

        def runner():
            try:
                return self.current_operator.run() if self.current_operator is not None else None
            except Exception as e:
                show_warning(f"Error during arc building: {str(e)}")
                self.set_enabled(True)
                self.current_operator = None

        worker = create_worker(
            runner,
            _progress={
                "desc": "Building arcs..."
            },
        )
        worker.finished.connect(self.finished_build_arcs)
        worker.start()

    def _distances_as_features(self, centrosome_id, df, scale):
        df = df[df['centrosome_id'] == centrosome_id]
        if df.empty:
            raise ValueError(f"No centrioles found for centrosome_id {centrosome_id}.")
        df = df.sort_values(by=['centriole_id', 'T'])
        df = df.iloc[len(df) // 2:]
        df = df['distance'] * scale['X'][0]
        df = df.reset_index(drop=True)
        return df

    def _add_centrosome_distance(self, centrosome_id, df, scale):
        distances = self._distances_as_features(centrosome_id, df, scale)
        centrosome_layer_name = self.tracks_manager_widget.as_track_layer_name(centrosome_id)
        centrosome_layer = self.viewer.layers[centrosome_layer_name]
        centrosome_layer.features['distance'] = np.round(distances, 2)
        unit = scale['X'][2]
        text = {
            'string': '{distance}'+format(unit, '~'),
            'anchor': 'center',
            'size': 10,
            'color': 'white',
        }
        centrosome_layer.text = text

    def finished_build_arcs(self, *args):
        if self.current_operator is None:
            return
        
        self.set_enabled(True)
        if not isinstance(self.current_operator, BuildArcsOperator):
            raise ValueError("Current operator is not a BuildArcsOperator.")
        
        arcs = self.current_operator.get_arcs()
        pairs = self.current_operator.get_pairs(flipped=True)
        df = self.current_operator.get_centrosomes()

        image_layer_name = self.image_combo.currentText()
        image_layer = self.viewer.layers[image_layer_name]
        image_layer.metadata['centrioles-pairs'] = pairs
        dims = self._get_image_dims()

        centrioles_layer_name = self.centrioles_tracks_prefix + image_layer_name
        centrioles_layer = self.viewer.layers[centrioles_layer_name]
        centrioles_layer.features = df

        for centriole_id, arc in arcs.items():
            centrosome_id = int(pairs[centriole_id])
            color = self.tracks_manager_widget._rows[centrosome_id].color
            layer_name = f"{self.arcs_shapes_prefix}{centriole_id}"
            if layer_name in self.viewer.layers:
                layer = self.viewer.layers[layer_name]
                self.viewer.layers.remove(layer)
            self._add_centrosome_distance(centrosome_id, df, dims)
            self.viewer.add_shapes(
                arc,
                shape_type='path',
                edge_color=color,
                edge_width=1,
                face_color='transparent',
                opacity=0.75,
                name=f"{self.arcs_shapes_prefix}{centriole_id}",
                scale=image_layer.scale,
                units=image_layer.units
            )

        self.current_operator = None

    def _gather_arcs(self):
        arcs = {}
        for layer in self.viewer.layers:
            if layer.name.startswith(self.arcs_shapes_prefix):
                centriole_id = int(layer.name.replace(self.arcs_shapes_prefix, ""))
                arcs[centriole_id] = np.array(layer.data)
        return arcs

    def launch_build_kymographs(self):
        image_layer_name = self.image_combo.currentText()
        if image_layer_name not in self.viewer.layers:
            show_warning("Selected image layer not found.")
            return
        
        centrosomes_layer_name = self.centrioles_tracks_prefix + image_layer_name
        if centrosomes_layer_name not in self.viewer.layers:
            show_warning("No tracked centrosomes layer found. Please run 'Find centrosomes' first.")
            return

        image_layer = self.viewer.layers[image_layer_name]
        image = image_layer.data
        scale = {a: s for a, s in zip(image_layer.axis_labels, image_layer.scale)}
        units = {a: u for a, u in zip(image_layer.axis_labels, image_layer.units)}

        arcs = self._gather_arcs()
        centrosomes = self.viewer.layers[centrosomes_layer_name].features
        self.set_enabled(False)

        op = MakeKymographOperator()
        op.set_input_image(image, scale, units)
        op.set_arcs(arcs)
        op.set_centrosomes(centrosomes)

        self.current_operator = op

        def runner():
            try:
                return self.current_operator.run() if self.current_operator is not None else None
            except Exception as e:
                show_warning(f"Error during kymograph building: {str(e)}")
                self.set_enabled(True)
                self.current_operator = None

        worker = create_worker(
            runner,
            _progress={
                "desc": "Building kymographs..."
            },
        )
        worker.finished.connect(self.finished_build_kymographs)
        worker.start()

    def _get_image_dims(self):
        image_layer_name = self.image_combo.currentText()
        if image_layer_name not in self.viewer.layers:
            raise ValueError("Selected image layer not found.")
        image_layer = self.viewer.layers[image_layer_name]
        dims = {a: (c, s, u) for a, c, s, u in zip(image_layer.axis_labels, image_layer.scale, image_layer.data.shape, image_layer.units)}
        return dims

    def finished_build_kymographs(self, *args):
        if self.current_operator is None:
            return
        
        self.set_enabled(True)

        if not isinstance(self.current_operator, MakeKymographOperator):
            raise ValueError("Current operator is not a MakeKymographOperator.")

        image_layer_name = self.image_combo.currentText()
        image_layer = self.viewer.layers[image_layer_name]
        pairs = image_layer.metadata.get('centrioles-pairs', None)

        if pairs is None:
            raise ValueError("Centrioles pairs not found in image layer metadata. Please run 'Build arcs' first.")

        centriole_ids = [(k, v) for k, v in pairs.items()]
        centriole_ids.sort(key=lambda x: x[1]) # sorted by centrosome_id
        centrosome_ids = [int(v) for _, v in centriole_ids]
        centriole_ids = [int(k) for k, _ in centriole_ids]

        kymographs = self.current_operator.get_kymographs()
        padding = 10
        dims = self._get_image_dims()
        scale = (dims['Y'][0], dims['X'][0])
        units = (dims['Y'][2], dims['X'][2])
        start_x = dims['X'][1] + padding

        # create a list of polygons
        polygons = []
        for i, centriole_id in enumerate(centriole_ids):
            kymograph = kymographs[centriole_id]
            T, Y = kymograph.shape
            polygon = np.array([
                [0, i * (Y + padding) + start_x],
                [T - 1, i * (Y + padding) + start_x],
                [T - 1, i * (Y + padding) + Y - 1 + start_x],
                [0, i * (Y + padding) + Y - 1 + start_x],
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
            'size': 10,
            'color': 'white',
        }

        colors = [self.tracks_manager_widget._rows[int(pairs[centriole_id])].color for centriole_id in centriole_ids]

        if "kymo_outlines" in self.viewer.layers:
            layer = self.viewer.layers["kymo_outlines"]
            self.viewer.layers.remove(layer)

        self.viewer.add_shapes(
            polygons,
            opacity=1.0,
            scale=scale,
            units=units,
            features=features,
            shape_type='polygon',
            edge_width=3,
            edge_color=colors,
            face_color='transparent',
            text=text,
            name='kymo_outlines'
        )

        # add the images
        for i, centriole_id in enumerate(centriole_ids):
            kymograph = kymographs[centriole_id]
            name = f"{self.kymo_prefix}{centriole_id}"
            if name in self.viewer.layers:
                layer = self.viewer.layers[name]
                self.viewer.layers.remove(layer)
            self.viewer.add_image(
                kymograph,
                name=name,
                scale=scale,
                units=units,
                translate=(0, (i * (kymograph.shape[1] + padding) + start_x) * scale[1]),
                colormap='turbo'
            )

    def _gather_kymographs(self):
        kymographs = {}
        for layer in self.viewer.layers:
            if layer.name.startswith(self.kymo_prefix):
                centriole_id = int(layer.name.replace(self.kymo_prefix, ""))
                kymographs[centriole_id] = layer.data
        return kymographs

    def launch_locate_spots(self):
        all_kymos_layer_names = [layer.name for layer in self.viewer.layers if layer.name.startswith(CentrosomesWidget.kymo_prefix)]
        if not all_kymos_layer_names:
            show_warning("No kymograph layers found. Please run 'Build kymographs' first.")
            return
        
        kymographs = {int(layer.name.replace(CentrosomesWidget.kymo_prefix, '')): self.viewer.layers[layer.name].data for layer in self.viewer.layers if layer.name.startswith(CentrosomesWidget.kymo_prefix)}
        self.set_enabled(False)

        spots_promi = self.spot_prominence_input.value()

        op = FindSpotsOperator()
        op.set_std_coeff(spots_promi)
        op.set_kymographs(kymographs)

        self.current_operator = op

        def runner():
            try:
                return self.current_operator.run() if self.current_operator is not None else None
            except Exception as e:
                show_warning(f"Error during spot location: {str(e)}")
                self.set_enabled(True)
                self.current_operator = None

        worker = create_worker(
            runner,
            _progress={
                "desc": "Locating spots in kymographs..."
            },
        )
        worker.finished.connect(self.finished_locate_spots)
        worker.start()

    def finished_locate_spots(self, *args):
        if self.current_operator is None:
            return
        
        self.set_enabled(True)

        if not isinstance(self.current_operator, FindSpotsOperator):
            raise ValueError("Current operator is not a FindSpotsOperator.")

        coordinates = self.current_operator.get_coordinates()

        for centriole_id, coords in coordinates.items():
            if coords.size == 0:
                continue
            
            kymo_layer_name = f"{CentrosomesWidget.kymo_prefix}{centriole_id}"
            if kymo_layer_name not in self.viewer.layers:
                show_warning(f"Kymograph layer for track {centriole_id} not found.")
                continue
            
            kymo_layer = self.viewer.layers[kymo_layer_name]
            spots_layer_name = f"{CentrosomesWidget.spots_layer_prefix}{centriole_id}"
            if spots_layer_name in self.viewer.layers:
                layer = self.viewer.layers[spots_layer_name]
                self.viewer.layers.remove(layer)
            self.viewer.add_points(
                coords,
                name=f"{CentrosomesWidget.spots_layer_prefix}{centriole_id}",
                scale=kymo_layer.scale,
                units=kymo_layer.units,
                face_color='transparent',
                size=3,
                translate=kymo_layer.translate
            )

    def _gather_spots(self):
        spots = {}
        for layer in self.viewer.layers:
            if layer.name.startswith(self.spots_layer_prefix):
                centriole_id = int(layer.name.replace(self.spots_layer_prefix, ""))
                spots[centriole_id] = layer.data
        return spots

    def launch_export_summary(self):
        all_kymos_layer_names = [layer.name for layer in self.viewer.layers if layer.name.startswith(CentrosomesWidget.kymo_prefix)]
        if not all_kymos_layer_names:
            show_warning("No kymograph layers found. Please run 'Build kymographs' first.")
            return
        
        all_spots_layer_names = [layer.name for layer in self.viewer.layers if layer.name.startswith(CentrosomesWidget.spots_layer_prefix)]
        if not all_spots_layer_names:
            show_warning("No spots layers found. Please run 'Locate spots' first.")
            return

        tracked_centrosomes_layer_name = self.centrioles_tracks_prefix + self.image_combo.currentText()
        if tracked_centrosomes_layer_name not in self.viewer.layers:
            show_warning("No tracked centrosomes layer found. Please run 'Find centrosomes' first.")
            return
        
        kymographs = {int(layer.name.replace(CentrosomesWidget.kymo_prefix, '')): self.viewer.layers[layer.name].data for layer in self.viewer.layers if layer.name.startswith(CentrosomesWidget.kymo_prefix)}
        spots = {int(layer.name.replace(CentrosomesWidget.spots_layer_prefix, '')): self.viewer.layers[layer.name].data for layer in self.viewer.layers if layer.name.startswith(CentrosomesWidget.spots_layer_prefix)}
        centrosomes_df = self.viewer.layers[tracked_centrosomes_layer_name].features

        self.set_enabled(False)

        op = SummaryOperator()
        op.set_centrosomes(centrosomes_df)
        op.set_kymographs(kymographs)
        op.set_spots(spots)

        self.current_operator = op

        def runner():
            try:
                return self.current_operator.run() if self.current_operator is not None else None
            except Exception as e:
                show_warning(f"Error during summary export: {str(e)}")
                self.set_enabled(True)
                self.current_operator = None

        worker = create_worker(
            runner,
            _progress={
                "desc": "Exporting summary..."
            },
        )
        worker.finished.connect(self.finished_export_summary)
        worker.start()

    def finished_export_summary(self, *args):
        if self.current_operator is None:
            return
        
        self.set_enabled(True)

        if not isinstance(self.current_operator, SummaryOperator):
            raise ValueError("Current operator is not a SummaryOperator.")

        csv_path, _ = QFileDialog.getSaveFileName(
            self, 
            "Save Summary CSV", 
            self.image_combo.currentText() + ".csv", 
            "CSV Files (*.csv);;All Files (*)"
        )
        if not csv_path:
            show_info("Export cancelled.")
            return
        summary = self.current_operator.get_summary()
        summary.to_csv(csv_path, index=False)

    def launch_export_archive(self):
        archive_path = QFileDialog.getExistingDirectory(
            self,
            "Select Export Directory"
        )

        if not archive_path:
            show_info("Export cancelled.")
            return

        parent_folder = Path(archive_path)
        img_name_as_folder = self.image_combo.currentText().replace(".", "_")
        archive_path = parent_folder / img_name_as_folder
        
        op = ExportArchiveOperator(self, self.viewer)
        op.set_root_path(archive_path)
        
        self.current_operator = op

        def runner():
            try:
                return self.current_operator.run() if self.current_operator is not None else None
            except Exception as e:
                show_warning(f"Error during archive export: {str(e)}")
                self.set_enabled(True)
                self.current_operator = None

        worker = create_worker(
            runner,
            _progress={
                "desc": "Exporting archive..."
            },
        )
        worker.finished.connect(self.finished_export_archive)
        worker.start()

    def finished_export_archive(self, *args):
        if self.current_operator is None:
            return
        
        self.set_enabled(True)

        if not isinstance(self.current_operator, ExportArchiveOperator):
            raise ValueError("Current operator is not an ExportArchiveOperator.")

        show_info(f"Archive exported.")


def run():
    import tifffile as tiff
    from pathlib import Path
    from .import_archive_operator import ImportArchiveOperator

    viewer = napari.Viewer()
    widget = CentrosomesWidget(viewer=viewer)
    viewer.window.add_dock_widget(widget)

    # Loading the time series image and add it to the viewer
    folder_in = Path("/home/clement/Documents/projects/nucleation/3VPCs")
    filename  = "251119_#4_30_001_016.vsi - C561.tif"
    path_in   = folder_in / filename
    image     = tiff.imread(path_in)
    calib = {'T': 1, 'Y': 0.1083333, 'X': 0.1083333}
    units = {'T': 's', 'Y': 'um', 'X': 'um'}

    viewer.add_image(
        image, 
        name=filename, 
        scale=(calib['T'], calib['Y'], calib['X']),
        units=(units['T'], units['Y'], units['X']),
        axis_labels=('T', 'Y', 'X')
    )


    # archive_path = "/home/clement/Desktop/251119_#4_30_001_016_vsi - C561_tif.zip"
    # op = ImportArchiveOperator(viewer, widget)
    # op.set_root_path(archive_path)
    # op.run()

    id1 = widget.tracks_manager_widget.add_track(
        np.array([
            [48, 347, 90],
            [48, 337, 132]
        ])
    )
    widget.tracks_manager_widget.update_starting_frame(id1, 17)
    widget.tracks_manager_widget.update_ending_frame(id1, 463)

    id2 = widget.tracks_manager_widget.add_track(
        np.array([
            [78, 248, 477],
            [78, 237, 514]
        ])
    )
    widget.tracks_manager_widget.update_starting_frame(id2, 17)
    widget.tracks_manager_widget.update_ending_frame(id2, 463)

    # hints = np.array([
    #     [17, 347, 90],
    #     [17, 337, 132],
    #     [17, 248, 477],
    #     [17, 237, 514]
    # ])

    # hints_layer = viewer.add_points(
    #     hints, 
    #     name="hints",
    #     scale=(1, calib, calib),
    # )

    # widget.time_range_start_input.setValue(17)
    # widget.time_range_end_input.setValue(463)

    # dump = Path("/home/clement/Documents/projects/nucleation/draft/implementation/dump")
    # centrosomes_path = dump / "centrosome_points_track.csv"
    # centrosomes = pd.read_csv(centrosomes_path)

    # viewer.add_tracks(
    #     centrosomes[['track_id', 'T', 'Y', 'X']],
    #     name=CentrosomesWidget.centrosomes_tracks_prefix + im_layer_name,
    #     scale=(1, calib, calib),
    #     features=centrosomes,
    #     graph=None,
    #     tail_length=4,
    #     hide_completed_tracks=True,
    #     tail_width=3,
    # )
    # viewer.add_points(
    #     centrosomes[['T', 'Y', 'X']].values,
    #     name=CentrosomesWidget.centrosomes_points_prefix + im_layer_name,
    #     scale=(1, calib, calib),
    #     properties=centrosomes.drop(columns=['T', 'Y', 'X']).to_dict(orient='list'),
    #     border_color='track_id',
    #     face_color='transparent',
    #     size=15,
    #     border_color_cycle=widget.colors[:len(centrosomes['track_id'].unique())]
    # )

    # arcs = {}
    # arcs_folder = dump / "arcs"
    # arcs_content = [f for f in arcs_folder.iterdir() if f.is_file() and f.suffix == ".npy"]
    # for arc_file in arcs_content:
    #     track_id = int(arc_file.name.replace(".npy", ""))
    #     arc = np.load(arc_file)
    #     arcs[track_id] = arc
    
    # as_napari_shapes, features = BuildArcsOperator.as_napari_shapes(arcs, time_start=17)
    # viewer.add_shapes(
    #     as_napari_shapes,
    #     name=CentrosomesWidget.arcs_shapes_prefix + im_layer_name,
    #     shape_type='path',
    #     edge_color='track_id',
    #     features=features,
    #     edge_width=2,
    #     edge_color_cycle=widget.colors[:len(arcs)],
    #     scale=(1, calib, calib),
    #     metadata={'arcs': arcs}
    # )

    napari.run()


def widget_only():
    viewer = napari.Viewer()
    widget = CentrosomesWidget(viewer=viewer)
    viewer.window.add_dock_widget(widget)

    napari.run()


if __name__ == "__main__":
    run()
