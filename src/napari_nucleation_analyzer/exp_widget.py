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
from qtpy.QtCore import Qt
import napari

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


class NucleationWidget(QWidget):

    centrosomes_tracks_prefix = "centrosomes_tracks_"
    centrosomes_points_prefix = "centrosomes_points_"
    arcs_shapes_prefix = "arcs_"
    kymo_prefix = "kymograph_"
    spots_layer_prefix = "spots_kymo_"

    colors = [
        'red', 'lime', 'cyan', 'yellow', 'magenta', 'orange', 'purple', 'brown', 'pink',
        'green', 'teal', 'navy', 'maroon', 'olive', 'gray', 'black', 'white', 'silver', 
        'violet', 'turquoise', 'salmon', 'blue', 'gold', 'indigo'
    ]

    def __init__(self, viewer: "napari.viewer.Viewer"):
        super().__init__()
        self.viewer = viewer
        self.image_comboboxes = []
        self.label_comboboxes = []
        self.shape_comboboxes = []
        self.points_comboboxes = []
        self.track_comboboxes = []
        self.custom_font = QFont()
        self.tracks_manager = None
        self.custom_font.setFamily("Arial Unicode MS, Segoe UI Emoji, Apple Color Emoji, Noto Color Emoji")
        self.init_ui()
        self.current_operator = None
        self.init_callbacks()

    # -------- UI: ----------------------------------

    def init_ui(self):
        layout = QVBoxLayout()
        self.manage_tracks_panel(layout)
        # self.find_centrosomes_panel(layout)
        # self.build_arcs_panel(layout)
        # self.kymographs_panel(layout)
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

    def manage_tracks_panel(self, parent_layout):
        self.tracks_manager = TracksManagerWidget(self.viewer)
        parent_layout.addWidget(self.tracks_manager)

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

        # Hints
        h_layout = QHBoxLayout()
        h_layout.addWidget(QLabel("Hints"))
        self.hints_combo = QComboBox()
        self.points_comboboxes.append(self.hints_combo)
        h_layout.addWidget(self.hints_combo)
        layout.addLayout(h_layout)

        # Time range
        h_layout = QHBoxLayout()
        h_layout.addWidget(QLabel("Start"))
        self.time_range_start_input = QSpinBox()
        self.time_range_start_input.setRange(0, 1000000)
        h_layout.addWidget(self.time_range_start_input)
        h_layout.addWidget(QLabel("End"))
        self.time_range_end_input = QSpinBox()
        self.time_range_end_input.setRange(0, 1000000)
        h_layout.addWidget(self.time_range_end_input)
        layout.addLayout(h_layout)

        # Find centrosomes
        self.find_centrosomes_button = QPushButton("Find centrosomes")
        self.find_centrosomes_button.setFont(self.custom_font)
        self.find_centrosomes_button.clicked.connect(self.launch_find_centrosomes)
        layout.addWidget(self.find_centrosomes_button)

        self.find_centrosomes_group.setLayout(layout)
        parent_layout.addWidget(self.find_centrosomes_group)

    def build_arcs_panel(self, parent_layout):
        self.build_arcs_group = QGroupBox("Build arcs")
        layout = QVBoxLayout()

        # Radius (µm)
        h_layout = QHBoxLayout()
        h_layout.addWidget(QLabel("Radius (µm)"))
        self.radius_input = QDoubleSpinBox()
        self.radius_input.setRange(0.0, 100.0)
        self.radius_input.setValue(2.0)
        h_layout.addWidget(self.radius_input)
        layout.addLayout(h_layout)

        # Angle (°)
        h_layout = QHBoxLayout()
        h_layout.addWidget(QLabel("Angle (°)"))
        self.angle_input = QDoubleSpinBox()
        self.angle_input.setRange(0.0, 360.0)
        self.angle_input.setValue(90.0)
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
        self.locate_spots_button = QPushButton("Locate spots")
        self.locate_spots_button.setFont(self.custom_font)
        self.locate_spots_button.clicked.connect(self.launch_locate_spots)
        layout.addWidget(self.locate_spots_button)

        # Export summary
        self.export_summary_button = QPushButton("Export summary")
        self.export_summary_button.setFont(self.custom_font)
        self.export_summary_button.clicked.connect(self.launch_export_summary)
        layout.addWidget(self.export_summary_button)

        self.kymographs_group.setLayout(layout)
        parent_layout.addWidget(self.kymographs_group)

    # -------- Callbacks: ----------------------------------

    def set_enabled(self, enabled: bool):
        self.find_centrosomes_group.setEnabled(enabled)
        self.build_arcs_group.setEnabled(enabled)
        self.kymographs_group.setEnabled(enabled)

    # -------- Utils: --------------------------------------

    def launch_find_centrosomes(self):
        image_layer_name = self.image_combo.currentText()
        hints_layer_name = self.hints_combo.currentText()
        time_start = self.time_range_start_input.value()
        time_end = self.time_range_end_input.value()

        if image_layer_name not in self.viewer.layers:
            show_warning("Selected image layer not found.")
            return
        
        if hints_layer_name not in self.viewer.layers:
            show_warning("Selected hints layer not found.")
            return
        
        if time_start < 0 or time_end >= self.viewer.layers[image_layer_name].data.shape[0] or time_start > time_end:
            show_warning("Invalid time range.")
            return
        
        self.set_enabled(False)
        op = FindCentrosomesOperator()
        op.set_input_image(
            self.viewer.layers[image_layer_name].data, 
            self.viewer.layers[image_layer_name].scale
        )
        op.set_hint_points(self.viewer.layers[hints_layer_name].data)
        op.set_time_range(time_start, time_end)
        self.current_operator = op

        worker = create_worker(
            self.current_operator.run,
            _progress={
                "desc": "Detecting and tracking centrosomes..."
            },
        )
        worker.finished.connect(self.finished_find_centrosomes)
        worker.start()

    def finished_find_centrosomes(self, *args):
        if self.current_operator is None:
            raise ValueError("No operator is currently running.")
        
        self.set_enabled(True)
        image_layer_name = self.image_combo.currentText()
        centrosomes_tracks_layer_name = self.centrosomes_tracks_prefix + image_layer_name
        centrosomes_points_layer_name = self.centrosomes_points_prefix + image_layer_name
        result = self.current_operator.get_centrosomes()

        if centrosomes_tracks_layer_name in self.viewer.layers:
            layer = self.viewer.layers[centrosomes_tracks_layer_name]
            layer.data = result[['track_id', 'T', 'Y', 'X']]
            layer.features = result
        else:
            self.viewer.add_tracks(
                result[['track_id', 'T', 'Y', 'X']],
                name=centrosomes_tracks_layer_name,
                scale=self.viewer.layers[image_layer_name].scale,
                features=result,
                graph=None,
                tail_length=4,
                hide_completed_tracks=True,
                tail_width=3,
            )

        if centrosomes_points_layer_name in self.viewer.layers:
            layer = self.viewer.layers[centrosomes_points_layer_name]
            layer.data = result[['T', 'Y', 'X']].values
            layer.features = result
            layer.border_color = 'track_id'
            layer.face_color = 'transparent'
            layer.size = 15
            layer.border_color_cycle = self.colors[:len(result['track_id'].unique())]
        else:
            self.viewer.add_points(
                result[['T', 'Y', 'X']].values,
                name=centrosomes_points_layer_name,
                scale=self.viewer.layers[image_layer_name].scale,
                properties=result.drop(columns=['T', 'Y', 'X']).to_dict(orient='list'),
                border_color='track_id',
                face_color='transparent',
                size=15,
                border_color_cycle=self.colors[:len(result['track_id'].unique())]
            )
        self.current_operator = None

    def launch_build_arcs(self):
        tracked_centrosomes_layer_name = self.centrosomes_tracks_prefix + self.image_combo.currentText()
        if tracked_centrosomes_layer_name not in self.viewer.layers:
            show_warning("No tracked centrosomes layer found. Please run 'Find centrosomes' first.")
            return
        tracked_centrosomes_layer = self.viewer.layers[tracked_centrosomes_layer_name]
        tracked_centrosomes = tracked_centrosomes_layer.features

        hints_layer_name = self.hints_combo.currentText()
        if hints_layer_name not in self.viewer.layers:
            show_warning("Selected hints layer not found.")
            return
        hints_layer = self.viewer.layers[hints_layer_name]

        angle = self.angle_input.value()
        radius = self.radius_input.value()
        calib = tracked_centrosomes_layer.scale
        hints = hints_layer.data
        start = self.time_range_start_input.value()

        self.set_enabled(False)
        op = BuildArcsOperator()
        op.set_angle(angle)
        op.set_radius(radius)
        op.set_calibration(calib)
        op.set_hints(hints)
        op.set_centrosomes(tracked_centrosomes)
        op.set_time_start(start)
        self.current_operator = op

        worker = create_worker(
            self.current_operator.run,
            _progress={
                "desc": "Building arcs..."
            },
        )
        worker.finished.connect(self.finished_build_arcs)
        worker.start()

    def finished_build_arcs(self, *args):
        if self.current_operator is None:
            raise ValueError("No operator is currently running.")
        
        self.set_enabled(True)
        arcs = self.current_operator.get_arcs()
        time_start = self.current_operator.time_start
        as_napari_shapes, features = BuildArcsOperator.as_napari_shapes(arcs, time_start=time_start)
        arcs_shapes_layer_name = self.arcs_shapes_prefix + self.image_combo.currentText()
        
        if arcs_shapes_layer_name in self.viewer.layers:
            layer = self.viewer.layers[arcs_shapes_layer_name]
            layer.data = as_napari_shapes
            layer.features = features
            layer.edge_color = 'track_id'
            layer.edge_color_cycle = self.colors[:len(arcs)]
            layer.metadata['arcs'] = arcs
        else:
            self.viewer.add_shapes(
                as_napari_shapes,
                name=arcs_shapes_layer_name,
                shape_type='path',
                edge_color='track_id',
                features=features,
                edge_width=2,
                edge_color_cycle=self.colors[:len(arcs)],
                scale=self.viewer.layers[self.image_combo.currentText()].scale,
                metadata={'arcs': arcs}
            )
        self.current_operator = None

    def launch_build_kymographs(self):
        image_layer_name = self.image_combo.currentText()
        if image_layer_name not in self.viewer.layers:
            show_warning("Selected image layer not found.")
            return
        
        centrosomes_layer_name = self.centrosomes_tracks_prefix + image_layer_name
        if centrosomes_layer_name not in self.viewer.layers:
            show_warning("No tracked centrosomes layer found. Please run 'Find centrosomes' first.")
            return
        
        arcs_layer_name = self.arcs_shapes_prefix + image_layer_name
        if arcs_layer_name not in self.viewer.layers:
            show_warning("No arcs layer found. Please run 'Build arcs' first.")
            return
        
        t_start = self.time_range_start_input.value()
        t_end = self.time_range_end_input.value()
        image = self.viewer.layers[image_layer_name].data
        arcs = self.viewer.layers[arcs_layer_name].metadata['arcs']
        centrosomes = self.viewer.layers[centrosomes_layer_name].features

        op = MakeKymographOperator()
        self.set_enabled(False)
        self.current_operator = op

        op.set_input_image(image[t_start:t_end+1])
        op.set_arcs(arcs)
        op.set_centrosomes(centrosomes)

        worker = create_worker(
            self.current_operator.run,
            _progress={
                "desc": "Building kymographs..."
            },
        )
        worker.finished.connect(self.finished_build_kymographs)
        worker.start()

    def finished_build_kymographs(self, *args):
        if self.current_operator is None:
            raise ValueError("No operator is currently running.")
        
        self.set_enabled(True)
        kymographs = self.current_operator.get_kymographs()
        padding    = 5
        track_ids  = []

        # add the images
        for track_id, kymograph in kymographs.items():
            self.viewer.add_image(
                kymograph,
                name=f"{CentrosomesWidget.kymo_prefix}{track_id}",
                scale=(1, 1),
                translate=(0, track_id * (kymograph.shape[1] + padding)),
                colormap='turbo'
            )
            track_ids.append(track_id)

        # create a list of polygons
        polygons = []
        for track_id, kymograph in kymographs.items():
            T, Y = kymograph.shape
            polygon = np.array([
                [0, track_id * (Y + padding)],
                [T - 1, track_id * (Y + padding)],
                [T - 1, track_id * (Y + padding) + Y - 1],
                [0, track_id * (Y + padding) + Y - 1],
            ])
            polygons.append(polygon)

        # create features
        features = {
            'track_id': track_ids
        }

        text = {
            'string': 'ID: {track_id}',
            'anchor': 'upper_left',
            'translation': [-5, 0],
            'size': 16,
            'color': 'white',
        }

        self.viewer.add_shapes(
            polygons,
            features=features,
            shape_type='polygon',
            edge_width=3,
            edge_color='track_id',
            edge_color_cycle=self.colors[:len(track_ids)],
            face_color='transparent',
            text=text,
            name='kymo_outlines'
        )

    def launch_locate_spots(self):
        all_kymos_layer_names = [layer.name for layer in self.viewer.layers if layer.name.startswith(CentrosomesWidget.kymo_prefix)]
        if not all_kymos_layer_names:
            show_warning("No kymograph layers found. Please run 'Build kymographs' first.")
            return
        
        kymographs = {int(layer.name.replace(CentrosomesWidget.kymo_prefix, '')): self.viewer.layers[layer.name].data for layer in self.viewer.layers if layer.name.startswith(CentrosomesWidget.kymo_prefix)}
        self.set_enabled(False)

        op = FindSpotsOperator()
        self.current_operator = op

        op.set_kymographs(kymographs)

        worker = create_worker(
            self.current_operator.run,
            _progress={
                "desc": "Locating spots in kymographs..."
            },
        )
        worker.finished.connect(self.finished_locate_spots)
        worker.start()

    def finished_locate_spots(self, *args):
        if self.current_operator is None:
            raise ValueError("No operator is currently running.")
        
        self.set_enabled(True)
        coordinates = self.current_operator.get_coordinates()

        for track_id, coords in coordinates.items():
            if coords.size == 0:
                continue
            
            kymo_layer_name = f"{CentrosomesWidget.kymo_prefix}{track_id}"
            if kymo_layer_name not in self.viewer.layers:
                show_warning(f"Kymograph layer for track {track_id} not found.")
                continue
            
            kymo_layer = self.viewer.layers[kymo_layer_name]
            self.viewer.add_points(
                coords,
                name=f"{CentrosomesWidget.spots_layer_prefix}{track_id}",
                scale=kymo_layer.scale,
                face_color='transparent',
                size=5,
                translate=(0, track_id * (kymo_layer.data.shape[1] + 5))
            )

    def launch_export_summary(self):
        all_kymos_layer_names = [layer.name for layer in self.viewer.layers if layer.name.startswith(CentrosomesWidget.kymo_prefix)]
        if not all_kymos_layer_names:
            show_warning("No kymograph layers found. Please run 'Build kymographs' first.")
            return
        
        all_spots_layer_names = [layer.name for layer in self.viewer.layers if layer.name.startswith(CentrosomesWidget.spots_layer_prefix)]
        if not all_spots_layer_names:
            show_warning("No spots layers found. Please run 'Locate spots' first.")
            return
        
        kymographs = {int(layer.name.replace(CentrosomesWidget.kymo_prefix, '')): self.viewer.layers[layer.name].data for layer in self.viewer.layers if layer.name.startswith(CentrosomesWidget.kymo_prefix)}
        spots = {int(layer.name.replace(CentrosomesWidget.spots_layer_prefix, '')): self.viewer.layers[layer.name].data for layer in self.viewer.layers if layer.name.startswith(CentrosomesWidget.spots_layer_prefix)}

        self.set_enabled(False)
        op = SummaryOperator()
        self.current_operator = op

        op.set_kymographs(kymographs)
        op.set_spots(spots)

        worker = create_worker(
            self.current_operator.run,
            _progress={
                "desc": "Exporting summary..."
            },
        )
        worker.finished.connect(self.finished_export_summary)
        worker.start()

    def finished_export_summary(self, *args):
        self.set_enabled(True)
        csv_path, _ = QFileDialog.getSaveFileName(self, "Save Summary CSV", "", "CSV Files (*.csv);;All Files (*)")
        if not csv_path:
            show_info("Export cancelled.")
            return
        summary = self.current_operator.get_summary()
        summary.to_csv(csv_path, index=False)


def run():
    import tifffile as tiff
    from pathlib import Path

    viewer = napari.Viewer()
    widget = CentrosomesWidget(viewer=viewer)
    viewer.window.add_dock_widget(widget)

    folder_in = Path("/home/clement/Documents/projects/nucleation/3VPCs")
    filename  = "251119_#4_30_001_016.vsi - C561.tif"
    path_in   = folder_in / filename
    calib     = 0.1083333 # µm/pixel
    image     = tiff.imread(path_in)

    image_layer = viewer.add_image(
        image, 
        name=filename, 
        scale=(1, calib, calib)
    )
    im_layer_name = image_layer.name

    hints = np.array([
        [17, 347, 90],
        [17, 337, 132],
        [17, 248, 477],
        [17, 237, 514]
    ])

    hints_layer = viewer.add_points(
        hints, 
        name="hints",
        scale=(1, calib, calib),
    )

    widget.time_range_start_input.setValue(17)
    widget.time_range_end_input.setValue(463)

    dump = Path("/home/clement/Documents/projects/nucleation/draft/implementation/dump")
    centrosomes_path = dump / "centrosome_points_track.csv"
    centrosomes = pd.read_csv(centrosomes_path)

    viewer.add_tracks(
        centrosomes[['track_id', 'T', 'Y', 'X']],
        name=CentrosomesWidget.centrosomes_tracks_prefix + im_layer_name,
        scale=(1, calib, calib),
        features=centrosomes,
        graph=None,
        tail_length=4,
        hide_completed_tracks=True,
        tail_width=3,
    )
    viewer.add_points(
        centrosomes[['T', 'Y', 'X']].values,
        name=CentrosomesWidget.centrosomes_points_prefix + im_layer_name,
        scale=(1, calib, calib),
        properties=centrosomes.drop(columns=['T', 'Y', 'X']).to_dict(orient='list'),
        border_color='track_id',
        face_color='transparent',
        size=15,
        border_color_cycle=widget.colors[:len(centrosomes['track_id'].unique())]
    )

    arcs = {}
    arcs_folder = dump / "arcs"
    arcs_content = [f for f in arcs_folder.iterdir() if f.is_file() and f.suffix == ".npy"]
    for arc_file in arcs_content:
        track_id = int(arc_file.name.replace(".npy", ""))
        arc = np.load(arc_file)
        arcs[track_id] = arc
    
    as_napari_shapes, features = BuildArcsOperator.as_napari_shapes(arcs, time_start=17)
    viewer.add_shapes(
        as_napari_shapes,
        name=CentrosomesWidget.arcs_shapes_prefix + im_layer_name,
        shape_type='path',
        edge_color='track_id',
        features=features,
        edge_width=2,
        edge_color_cycle=widget.colors[:len(arcs)],
        scale=(1, calib, calib),
        metadata={'arcs': arcs}
    )

    napari.run()


def widget_only():
    import tifffile as tiff

    viewer = napari.Viewer()
    widget = NucleationWidget(viewer=viewer)
    viewer.window.add_dock_widget(widget)

    path = "/home/clement/Documents/projects/nucleation/3VPCs/251119_#4_30_001_016.vsi - C561.tif"
    img_2dt = tiff.imread(path)

    viewer.add_image(
        img_2dt,
        name="2D+t image"
    )
    
    napari.run()


if __name__ == "__main__":
    widget_only()
