from __future__ import annotations

from dataclasses import dataclass

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

    color: str
    index: int
    container: QWidget
    label: QLabel
    start_button: QPushButton
    end_button: QPushButton
    delete_button: QPushButton
    frame_start: int | None = None
    frame_end: int | None = None


class TracksManagerWidget(QWidget):
    """Widget managing a dynamic table of track rows plus an "Add track" button."""

    trackStartRequested = Signal(int)
    trackEndRequested = Signal(int)
    trackDeleteRequested = Signal(int)

    INDEX_PADDING = 2
    ACTION_BUTTON_WIDTH = 110
    DELETE_BUTTON_WIDTH = 28
    PREFIX = "_Centrosome "

    PALETTE = [
        "#1f77b4",  # muted blue
        "#ff7f0e",  # safety orange
        "#2ca02c",  # cooked asparagus green
        "#d62728",  # brick red
        "#9467bd",  # muted purple
        "#8c564b",  # chestnut brown
        "#e377c2",  # raspberry yogurt pink
        "#7f7f7f",  # middle gray
        "#bcbd22",  # curry yellow-green
        "#17becf",  # blue-teal
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
        self._rows_layout.setSpacing(2)
        self._main_layout.addWidget(self._rows_container)

        self._add_button = QPushButton("Add centrosome", self)
        self._add_button.clicked.connect(self._on_add_track_clicked)
        self._main_layout.addWidget(self._add_button)

        self._main_layout.addStretch(1)

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def add_track(self, line=None) -> int:
        """Add a new track row. Returns the index assigned to the new track."""
        index = self._next_index
        self._next_index += 1

        row = self._build_row(index)
        self._rows[index] = row
        self._rows_layout.addWidget(row.container)

        self._add_track_layer(index, line)

        return index

    def remove_track(self, index: int) -> None:
        """Remove the row corresponding to `index`, if it exists."""
        row = self._rows.pop(index, None)
        if row is None:
            return  # nothing to do, already removed or never existed

        self._rows_layout.removeWidget(row.container)
        row.container.deleteLater()
        self._remove_track_layer(index)

    def track_count(self) -> int:
        """Number of currently displayed tracks."""
        return len(self._rows)

    def update_starting_frame(self, index: int, frame: int) -> Tuple[int, int] | None:
        """Set the starting frame for a track row, updating its label."""
        row = self._rows.get(index)
        if row is None:
            return None
        end_frame = row.frame_end if row.frame_end is not None else frame
        row.frame_start = frame
        row.frame_end = max(end_frame, frame)
        self._freeze_line_on_layer(index, frame)
        return (row.frame_start, row.frame_end)

    def update_ending_frame(self, index: int, frame: int) -> Tuple[int, int] | None:
        """Set the ending frame for a track row, updating its label."""
        row = self._rows.get(index)
        if row is None:
            return None
        start_frame = row.frame_start if row.frame_start is not None else frame
        row.frame_end = frame
        row.frame_start = min(start_frame, frame)
        self._freeze_line_on_layer(index, start_frame)
        return (row.frame_start, row.frame_end)

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
                if len(layer.data) > 0:
                    line = layer.data[0]
                    # Extract spatial coordinates (skip the time dimension)
                    if line.shape[1] >= 3:
                        points = line[:, 1:]  # Remove time dimension
                    else:
                        points = line
            
            hints[index] = {
                'start': row.frame_start,
                'end': row.frame_end,
                'points': points,
                'color': row.color,
            }
        
        return hints

    # ------------------------------------------------------------------
    # Row construction
    # ------------------------------------------------------------------

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

    def _build_row(self, index: int) -> _TrackRow:
        container = QWidget(self._rows_container)
        layout = QHBoxLayout(container)
        layout.setContentsMargins(0, 0, 0, 0)

        color_dot = QLabel(container)
        color_dot.setFixedSize(12, 12)
        color_dot.setStyleSheet(f"background-color: {self._next_color()}; border-radius: 6px;")

        label = QLabel(self._format_label(index), container)

        start_button = QPushButton("Start", container)
        end_button = QPushButton("End", container)
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

        start_button.clicked.connect(lambda _checked, i=index, b1=start_button, b2=end_button: self._on_start_clicked(i, b1, b2))
        end_button.clicked.connect(lambda _checked, i=index, b1=start_button, b2=end_button: self._on_end_clicked(i, b1, b2))
        delete_button.clicked.connect(lambda _checked, i=index: self._on_delete_clicked(i))

        return _TrackRow(
            color=self._next_color(),
            index=index,
            container=container,
            label=label,
            start_button=start_button,
            end_button=end_button,
            delete_button=delete_button,
        )

    def _format_label(self, index: int) -> str:
        return f"Centrosome {index:0{self.INDEX_PADDING}d}:"

    # ------------------------------------------------------------------
    # Slots
    # ------------------------------------------------------------------

    def _add_track_layer(self, track_id, line=None):
        layer_name = f"{self.PREFIX}{track_id:0{self.INDEX_PADDING}d}"
        if layer_name in self.viewer.layers:
            return
        a, s, u = self.parent.get_image_calibration()
        self.viewer.add_shapes(
            [] if line is None else [line],
            name=layer_name, 
            shape_type="line", 
            edge_color=self._next_color(), 
            edge_width=2, 
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

    def _on_start_clicked(self, index: int, start_button: QPushButton, end_button: QPushButton) -> None:
        f = int(self.viewer.dims.point[0])
        res = self.update_starting_frame(index, f)
        if res is None:
            return
        start, end = res
        start_button.setText(f"From: {start}")
        end_button.setText(f"To: {end}")
        self.trackStartRequested.emit(index)

    def _on_end_clicked(self, index: int, start_button: QPushButton, end_button: QPushButton) -> None:
        f = int(self.viewer.dims.point[0])
        res = self.update_ending_frame(index, f)
        if res is None:
            return
        start, end = res
        start_button.setText(f"From: {start}")
        end_button.setText(f"To: {end}")
        self.trackEndRequested.emit(index)

    def _on_delete_clicked(self, index: int) -> None:
        # Callback(s) run first (via the signal), THEN the row is removed,
        # as requested.
        self.trackDeleteRequested.emit(index)
        self.remove_track(index)