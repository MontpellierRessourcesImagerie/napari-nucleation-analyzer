from __future__ import annotations

from dataclasses import dataclass
from napari.utils.notifications import show_info, show_warning
from numpy.strings import index
from qtpy.QtCore import Signal
from typing import Tuple
from qtpy.QtWidgets import (
    QWidget,
    QVBoxLayout,
    QHBoxLayout,
    QLabel,
    QPushButton,
    QStyle,
    QSizePolicy,
)
import numpy as np


@dataclass
class _TrackRow:
    """Internal bookkeeping for a single track row's widgets."""

    color        : str
    index        : int
    container    : QWidget
    label        : QLabel
    start_button : QPushButton
    end_button   : QPushButton
    delete_button: QPushButton
    frame_start  : int | None = None
    frame_end    : int | None = None


class TracksManagerWidget(QWidget):
    """Widget managing a dynamic table of track rows plus an "Add track" button."""

    trackStartRequested  = Signal(int)
    trackEndRequested    = Signal(int)
    trackDeleteRequested = Signal(int)

    INDEX_PADDING       = 2
    ACTION_BUTTON_WIDTH = 110
    DELETE_BUTTON_WIDTH = 28
    PREFIX              = "_Centrosome "
    PALETTE             = [
        "#e01010",
        "#ff7f0e",
        "#15d5eb",
        "#1bdb1b",
        "#0073c5",
        "#f544c0",
        "#e6e618",
        "#ad59fc",
        "#6bff6b",
        "#eeba0f",
    ]

    def __init__(self, viewer, parent) -> None:
        super().__init__()
        self.viewer = viewer
        self.parent = parent
        self._next_index: int = 1
        self._rows: dict[int, _TrackRow] = {}
        self.create_ui()

    def create_ui(self):
        self._main_layout = QVBoxLayout(self)

        self._rows_container = QWidget(self)
        self._rows_layout = QVBoxLayout(self._rows_container)
        self._rows_layout.setSpacing(5)
        self._main_layout.addWidget(self._rows_container)

        self._add_button = QPushButton("Add centrosome", self)
        self._add_button.clicked.connect(self._on_add_track_clicked)

        self.clear_all_button = QPushButton("Clear all", self)
        self.clear_all_button.clicked.connect(self.clear_all)

        h_layout = QHBoxLayout()

        h_layout.addWidget(self._add_button)
        h_layout.addWidget(self.clear_all_button)

        self._main_layout.addLayout(h_layout)

    def add_track(self, line=None) -> int:
        """Add a new track row. Returns the index assigned to the new track."""
        index = self._next_index
        self._next_index += 1

        row = self._build_row(index)
        self._rows[index] = row
        self._rows_layout.addWidget(row.container)

        self._add_track_layer(index, line)

        return index

    def insert_track(self, index: int, line: list, f_start: int, f_end: int):
        """Insert a track row with a specific index, line, start/end frames, and color."""
        if index in self._rows:
            raise ValueError(f"Track with index {index} already exists.")

        self._next_index = max(self._next_index, index + 1)
        row = self._build_row(index)

        self._rows[index] = row
        self._rows_layout.addWidget(row.container)

        self._add_track_layer(index, line)

        self.set_starting_frame(index, f_start)
        self.set_ending_frame(index, f_end)

    def clear_all(self):
        """Remove all track rows and their corresponding layers."""
        for index in list(self._rows.keys()):
            self.remove_track(index)
        self._next_index = 1

    def remove_track(self, index: int) -> None:
        """Remove the row corresponding to `index`, if it exists."""
        row = self._rows.pop(index, None)
        if row is None:
            return

        self._rows_layout.removeWidget(row.container)
        row.container.deleteLater()
        self._remove_track_layer(index)

    def track_count(self) -> int:
        """Number of currently displayed tracks."""
        return len(self._rows)

    def update_starting_frame(self, index: int, frame: int):
        """Set the starting frame for a track row, updating its label."""
        if not self._has_a_line(index):
            return

        self._freeze_line_on_layer(index, frame)
        self.set_starting_frame(index, frame)

    def set_starting_frame(self, index: int, frame: int):

        row = self._rows.get(index)
        if row is None:
            return 

        end_frame = row.frame_end if row.frame_end is not None else frame
        row.frame_start = frame
        row.frame_end = max(end_frame, frame)
        
        start, end = (row.frame_start, row.frame_end)

        row.start_button.setText(f"From: {start}")
        row.end_button.setText(f"To: {end}")
        row.end_button.setEnabled(True)

    def update_ending_frame(self, index: int, frame: int):
        """Set the ending frame for a track row, updating its label."""
        if not self._has_a_line(index):
            return

        row = self._rows.get(index)
        if row is None:
            return

        start_frame = row.frame_start if row.frame_start is not None else frame
        self._freeze_line_on_layer(index, start_frame)
        self.set_ending_frame(index, frame)

    def set_ending_frame(self, index: int, frame: int):
        row = self._rows.get(index)
        if row is None:
            return
        
        start_frame = row.frame_start if row.frame_start is not None else frame
        row.frame_end = frame
        row.frame_start = min(start_frame, frame)
        start, end = (row.frame_start, row.frame_end)

        row.start_button.setText(f"From: {start}")
        row.end_button.setText(f"To: {end}")

    def as_track_layer_name(self, index: int) -> str:
        """Return the name of the track layer corresponding to `index`."""
        index = int(index)
        return f"{self.PREFIX}{index:0{self.INDEX_PADDING}d}"

    def as_hints(self) -> dict:
        """Convert tracks to hint format.
        
        Returns a dictionary where keys are track indices and values are dicts with:
        - 'start': starting frame
        - 'end': ending frame
        - 'points': numpy array of points (spatial coordinates)
        - 'color': color string
        """
        hints = {}
        for index, row in self._rows.items():
            # Get points from the layer if it exists
            layer_name = self.as_track_layer_name(index)
            points = np.array([])
            if layer_name in self.viewer.layers:
                layer = self.viewer.layers[layer_name]
                if len(layer.data) == 1:
                    line = layer.data[0]
                    if line.shape[1] >= 3:
                        points = line[:, 1:]  # Remove time dimension
                    else:
                        points = line
                else:
                    show_warning(f"Layer {layer_name} has {len(layer.data)} shapes, expected 1.")

            if row.frame_start is None or row.frame_end is None:
                show_warning(f"Track {index} has undefined start or end frame.")
                continue

            if row.frame_end - row.frame_start <= 0:
                show_warning(f"Track {index} has non-positive duration: start={row.frame_start}, end={row.frame_end}.")
                continue
            
            hints[index] = {
                'start': row.frame_start,
                'end': row.frame_end,
                'points': points,
                'color': row.color,
            }
        
        return hints

    def _next_color(self) -> str:
        """Return a color for the next track row, cycling through a palette."""
        return self.PALETTE[(self._next_index - 1) % len(self.PALETTE)]

    @staticmethod
    def _lock_width(button: QPushButton, width: int) -> None:
        """Force a button to keep a constant width, whatever its text/icon."""
        button.setFixedWidth(width)
        policy = button.sizePolicy()
        policy.setHorizontalPolicy(QSizePolicy.Fixed)
        button.setSizePolicy(policy)

    def _build_row(self, index: int, color=None) -> _TrackRow:
        container = QWidget(self._rows_container)
        layout = QHBoxLayout(container)
        layout.setContentsMargins(0, 0, 0, 0)

        color_dot = QLabel(container)
        color_dot.setFixedSize(12, 12)
        color_dot.setStyleSheet(f"background-color: {color or self._next_color()}; border-radius: 6px;")

        label = QLabel(self._format_label(index), container)

        start_button = QPushButton("Start", container)
        end_button = QPushButton("End", container)
        end_button.setEnabled(False)
        self._lock_width(start_button, self.ACTION_BUTTON_WIDTH)
        self._lock_width(end_button, self.ACTION_BUTTON_WIDTH)

        delete_button = QPushButton(container)
        delete_button.setIcon(self.style().standardIcon(QStyle.SP_DialogCancelButton))
        delete_button.setToolTip("Delete")
        self._lock_width(delete_button, self.DELETE_BUTTON_WIDTH)

        layout.addWidget(color_dot)
        layout.addWidget(label)
        layout.addWidget(start_button)
        layout.addWidget(end_button)
        layout.addWidget(delete_button)

        start_button.clicked.connect(lambda _checked, i=index: self._on_start_clicked(i))
        end_button.clicked.connect(lambda _checked, i=index: self._on_end_clicked(i))
        delete_button.clicked.connect(lambda _checked, i=index: self._on_delete_clicked(i))

        return _TrackRow(
            color=color or self._next_color(),
            index=index,
            container=container,
            label=label,
            start_button=start_button,
            end_button=end_button,
            delete_button=delete_button,
        )

    def _format_label(self, index: int) -> str:
        return f"Centrosome {index:0{self.INDEX_PADDING}d}:"

    def _add_track_layer(self, track_id, line=None):
        layer_name = f"{self.PREFIX}{track_id:0{self.INDEX_PADDING}d}"
        if layer_name in self.viewer.layers:
            return
        a, s, u = self.parent.get_image_calibration()
        return self.viewer.add_shapes(
            [] if line is None else line,
            name=layer_name, 
            shape_type="line", 
            edge_color=self._next_color(), 
            edge_width=1, 
            face_color='transparent', 
            opacity=0.75, 
            visible=True,
            ndim=3,
            axis_labels=a,
            scale=s,
            units=u
        )

    def _remove_track_layer(self, track_id):
        layer_name = f"{self.PREFIX}{track_id:0{self.INDEX_PADDING}d}"
        if layer_name in self.viewer.layers:
            self.viewer.layers.remove(self.viewer.layers[layer_name])

    def _on_add_track_clicked(self) -> None:
        layer = self.viewer.layers.selection.active
        if layer is None:
            show_warning("An active image layer is required.")
            return
        self.add_track()

    def _freeze_line_on_layer(self, index, start_t):
        layer_name = f"{self.PREFIX}{index:0{self.INDEX_PADDING}d}"
        if layer_name not in self.viewer.layers:
            return
        layer = self.viewer.layers[layer_name]
        line = layer.data[0]

        if line.shape[1] == 2:
            line = np.insert(line, 0, start_t, axis=1)
        else:
            line[0, 0] = start_t
            line[1, 0] = start_t
            
        layer.data = [line]

    def _has_a_line(self, index) -> bool:
        layer_name = f"{self.PREFIX}{index:0{self.INDEX_PADDING}d}"
        if layer_name not in self.viewer.layers:
            return False
        layer = self.viewer.layers[layer_name]
        d = layer.data
        if len(d) != 1:
            show_warning(f"Layer {layer_name} has {len(d)} shapes, expected 1.")
            return False
        line = layer.data[0]
        return line.shape[1] >= 3

    def _on_start_clicked(self, index: int) -> None:
        f = int(self.viewer.dims.point[0])
        self.update_starting_frame(index, f)
        self.trackStartRequested.emit(index)

    def _on_end_clicked(self, index: int) -> None:
        f = int(self.viewer.dims.point[0])
        self.update_ending_frame(index, f)
        self.trackEndRequested.emit(index)

    def _on_delete_clicked(self, index: int) -> None:
        # Callback(s) run first (via the signal), THEN the row is removed,
        # as requested.
        self.trackDeleteRequested.emit(index)
        self.remove_track(index)