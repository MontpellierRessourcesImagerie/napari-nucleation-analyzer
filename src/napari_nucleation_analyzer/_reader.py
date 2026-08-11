from .import_archive_operator import ImportArchiveOperator
from pathlib import Path
import napari

def _reader_function(archive_path):
    viewer = napari.current_viewer()
    op = ImportArchiveOperator(viewer)
    op.set_root_path(archive_path)
    op.run()

def napari_get_reader(path):
    if isinstance(path, list):
        path = path[0]

    path = Path(path)
    if not path.exists():
        return None

    extension = path.suffix.lower()
    if extension == ".zip":
        return _reader_function
    return None