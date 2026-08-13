import numpy as np
import pandas as pd
import pytest
import urllib.request
from pathlib import Path
import tifffile as tiff

TESTING_IMAGE_URL = "https://sdrive.cnrs.fr/s/oK2erWnk6xM9dzb/download"

from napari_nucleation_analyzer.operators.find_centrosomes_operator import FindCentrosomesOperator


# ----------------------------------------------------------------------
# Fixtures
# ----------------------------------------------------------------------

@pytest.fixture
def operator():
    return FindCentrosomesOperator()


@pytest.fixture
def calibration():
    return {"T": 1.0, "Y": 0.5, "X": 0.5}


@pytest.fixture
def units():
    return {"T": "s", "Y": "um", "X": "um"}


@pytest.fixture
def calibrated_operator(operator, calibration, units):
    img = np.zeros((10, 50, 50), dtype=np.float32)
    operator.set_input_image(img, calibration, units)
    return operator


# ----------------------------------------------------------------------
# Valeurs par défaut
# ----------------------------------------------------------------------

def test_default_values(operator):
    assert operator.prominence == FindCentrosomesOperator.default_prominence()
    assert operator.searching_range == FindCentrosomesOperator.default_searching_range()
    assert operator.memory == FindCentrosomesOperator.default_memory()
    assert operator.max_binding_distance == FindCentrosomesOperator.default_max_binding_distance()
    assert operator.input_image is None
    assert operator.hints == {}
    assert operator.centrosomes is None


# ----------------------------------------------------------------------
# set_input_image
# ----------------------------------------------------------------------

class TestSetInputImage:

    def test_valid_image(self, operator, calibration, units):
        img = np.zeros((5, 20, 20), dtype=np.float32)
        operator.set_input_image(img, calibration, units)
        assert operator.input_image is not None
        assert operator.input_image.dims == ("T", "Y", "X")
        assert operator.input_image.attrs["scale"] == calibration
        assert operator.input_image.attrs["units"] == units

    @pytest.mark.parametrize("shape", [(20, 20), (2, 5, 20, 20)])
    def test_wrong_ndim_raises(self, operator, calibration, units, shape):
        img = np.zeros(shape, dtype=np.float32)
        with pytest.raises(ValueError, match="3D array"):
            operator.set_input_image(img, calibration, units)


# ----------------------------------------------------------------------
# set_hints
# ----------------------------------------------------------------------

class TestSetHints:

    def test_requires_input_image_first(self, operator):
        hints = {1: {"start": 0, "end": 1, "points": np.array([[0, 0]])}}
        with pytest.raises(ValueError, match="Input image must be set"):
            operator.set_hints(hints)

    def test_valid_hints(self, calibrated_operator):
        hints = {
            1: {"start": 0, "end": 5, "points": np.array([[10, 10], [12, 12]])},
        }
        calibrated_operator.set_hints(hints)
        assert calibrated_operator.hints == hints

    @pytest.mark.parametrize("missing_key", ["start", "end", "points"])
    def test_missing_key_raises(self, calibrated_operator, missing_key):
        hint = {"start": 0, "end": 5, "points": np.array([[10, 10]])}
        del hint[missing_key]
        with pytest.raises(ValueError, match="keys are required"):
            calibrated_operator.set_hints({1: hint})

    def test_points_wrong_shape_raises(self, calibrated_operator):
        hint = {"start": 0, "end": 5, "points": np.array([1, 2, 3])}
        with pytest.raises(ValueError, match="shape \\(N, 2\\)"):
            calibrated_operator.set_hints({1: hint})

    def test_start_negative_raises(self, calibrated_operator):
        hint = {"start": -1, "end": 5, "points": np.array([[10, 10]])}
        with pytest.raises(ValueError, match="time range"):
            calibrated_operator.set_hints({1: hint})

    def test_end_out_of_range_raises(self, calibrated_operator):
        # calibrated_operator image has T=10, so valid end is 0..9
        hint = {"start": 0, "end": 10, "points": np.array([[10, 10]])}
        with pytest.raises(ValueError, match="time range"):
            calibrated_operator.set_hints({1: hint})


# ----------------------------------------------------------------------
# Setters numériques (prominence, searching_range, memory, max_binding_distance)
# ----------------------------------------------------------------------

class TestNumericSetters:

    @pytest.mark.parametrize("value", [0, 1, -0.1, 1.1])
    def test_prominence_invalid(self, operator, value):
        with pytest.raises(ValueError):
            operator.set_prominence(value)

    def test_prominence_valid(self, operator):
        operator.set_prominence(0.3)
        assert operator.prominence == 0.3

    @pytest.mark.parametrize("value", [0, -1.0])
    def test_searching_range_invalid(self, operator, value):
        with pytest.raises(ValueError):
            operator.set_searching_range(value)

    def test_searching_range_valid(self, operator):
        operator.set_searching_range(2.0)
        assert operator.searching_range == 2.0

    def test_memory_negative_invalid(self, operator):
        with pytest.raises(ValueError):
            operator.set_memory(-1)

    def test_memory_zero_is_valid(self, operator):
        # borne exacte : < 0 lève une erreur, mais 0 est accepté
        operator.set_memory(0)
        assert operator.memory == 0

    def test_memory_positive_valid(self, operator):
        operator.set_memory(5)
        assert operator.memory == 5

    @pytest.mark.parametrize("value", [0, -1.0])
    def test_max_binding_distance_invalid(self, operator, value):
        with pytest.raises(ValueError):
            operator.set_max_binding_distance(value)

    def test_max_binding_distance_valid(self, operator):
        operator.set_max_binding_distance(1.0)
        assert operator.max_binding_distance == 1.0


# ----------------------------------------------------------------------
# Getters de conversion pixel
# ----------------------------------------------------------------------

class TestPixelConversionGetters:

    def test_get_input_image_raises_if_unset(self, operator):
        with pytest.raises(ValueError):
            operator.get_input_image()

    def test_get_hints_raises_if_empty(self, operator):
        with pytest.raises(ValueError):
            operator.get_hints()

    def test_get_searching_range_pxl(self, calibrated_operator):
        calibrated_operator.set_searching_range(2.0)  # calibration X=0.5
        assert calibrated_operator.get_searching_range_pxl() == 4

    def test_get_memory_frames(self, calibrated_operator):
        calibrated_operator.set_memory(5.0)  # calibration T=1.0
        assert calibrated_operator.get_memory_frames() == 5

    def test_get_max_binding_distance_pxl(self, calibrated_operator):
        calibrated_operator.set_max_binding_distance(1.0)  # calibration X=0.5
        assert calibrated_operator.get_max_binding_distance_pxl() == 2

    def test_pixel_getters_raise_without_image(self, operator):
        with pytest.raises(ValueError):
            operator.get_searching_range_pxl()
        with pytest.raises(ValueError):
            operator.get_memory_frames()
        with pytest.raises(ValueError):
            operator.get_max_binding_distance_pxl()


# ----------------------------------------------------------------------
# _bind_tracks_to_hints
# ----------------------------------------------------------------------

class TestBindTracksToHints:

    def _make_operator_with_unit_scale(self):
        # calibration X=1.0 => distance en pixel == distance physique,
        # plus simple à raisonner dans les tests
        op = FindCentrosomesOperator()
        img = np.zeros((10, 100, 100), dtype=np.float32)
        op.set_input_image(img, {"T": 1.0, "Y": 1.0, "X": 1.0}, {"T": "s", "Y": "px", "X": "px"})
        op.set_max_binding_distance(5.0)
        return op

    def test_binds_closest_candidate(self):
        op = self._make_operator_with_unit_scale()
        op.hints = {
            1: {"start": 0, "end": 5, "points": np.array([[10.0, 10.0]])},
        }
        tracked = pd.DataFrame({
            "T": [0, 0, 0],
            "Y": [10.5, 50.0, 11.0],
            "X": [10.5, 50.0, 11.0],
            "centriole_id": [1, 2, 3],
        })
        bindings = op._bind_tracks_to_hints(tracked)
        # le centriole 1 (distance ~0.7) doit être choisi plutôt que le 3 (distance ~1.4)
        assert bindings[1][0] == 1

    def test_no_candidate_within_distance_gives_none(self):
        op = self._make_operator_with_unit_scale()
        op.hints = {
            1: {"start": 0, "end": 5, "points": np.array([[10.0, 10.0]])},
        }
        tracked = pd.DataFrame({
            "T": [0],
            "Y": [90.0],
            "X": [90.0],
            "centriole_id": [1],
        })
        bindings = op._bind_tracks_to_hints(tracked)
        assert bindings[1][0] is None


# ----------------------------------------------------------------------
# _filter_by_hint_points
# ----------------------------------------------------------------------

class TestFilterByHintPoints:

    def test_raises_if_binding_incomplete(self, calibrated_operator):
        calibrated_operator.set_max_binding_distance(1000.0)
        calibrated_operator.hints = {
            1: {"start": 0, "end": 5, "points": np.array([[1.0, 1.0], [2.0, 2.0]])},
        }
        # une seule ligne au T de start => un seul point pourra être bindé, l'autre non
        tracked = pd.DataFrame({
            "T": [0],
            "Y": [1.0],
            "X": [1.0],
            "centriole_id": [1],
        })
        with pytest.raises(ValueError, match="Failed to bind"):
            calibrated_operator._filter_by_hint_points(tracked)

    def test_filters_correctly_on_success(self, calibrated_operator):
        calibrated_operator.set_max_binding_distance(1000.0)
        calibrated_operator.hints = {
            1: {"start": 0, "end": 5, "points": np.array([[1.0, 1.0], [50.0, 50.0]])},
        }
        tracked = pd.DataFrame({
            "T": [0, 0, 0],
            "Y": [1.0, 50.0, 99.0],
            "X": [1.0, 50.0, 99.0],
            "centriole_id": [1, 2, 3],
        })
        result = calibrated_operator._filter_by_hint_points(tracked)
        assert set(result["centriole_id"]) == {1, 2}
        assert set(result["centrosome_id"]) == {1}


# ----------------------------------------------------------------------
# _interpolate_missing_time_points
# ----------------------------------------------------------------------

class TestInterpolateMissingTimePoints:

    def test_fills_gap_linearly(self, calibrated_operator):
        calibrated_operator.hints = {
            1: {"start": 0, "end": 4, "points": np.array([[0, 0]])},
        }
        tracked = pd.DataFrame({
            "T": [0, 2, 4],
            "Y": [0.0, 2.0, 4.0],
            "X": [0.0, 20.0, 40.0],
            "centriole_id": [1, 1, 1],
            "centrosome_id": [1, 1, 1],
        })
        result = calibrated_operator._interpolate_missing_time_points(tracked)
        result = result.sort_values("T").reset_index(drop=True)

        assert list(result["T"]) == [0, 1, 2, 3, 4]
        assert result.loc[result["T"] == 1, "Y"].iloc[0] == pytest.approx(1.0)
        assert result.loc[result["T"] == 1, "X"].iloc[0] == pytest.approx(10.0)
        assert result.loc[result["T"] == 3, "Y"].iloc[0] == pytest.approx(3.0)

    def test_trims_outside_hint_range(self, calibrated_operator):
        calibrated_operator.hints = {
            1: {"start": 1, "end": 3, "points": np.array([[0, 0]])},
        }
        tracked = pd.DataFrame({
            "T": [0, 1, 2, 3, 4],
            "Y": [0.0, 1.0, 2.0, 3.0, 4.0],
            "X": [0.0, 1.0, 2.0, 3.0, 4.0],
            "centriole_id": [1, 1, 1, 1, 1],
            "centrosome_id": [1, 1, 1, 1, 1],
        })
        result = calibrated_operator._interpolate_missing_time_points(tracked)
        assert set(result["T"]) == {1, 2, 3}


# ----------------------------------------------------------------------
# as_lines (staticmethod)
# ----------------------------------------------------------------------

class TestAsLines:

    def test_structure_and_default_color_fallback(self):
        df = pd.DataFrame({
            "T": [0, 1, 0, 1, 0, 1],
            "Y": [0, 1, 10, 11, 0, 1],
            "X": [0, 1, 10, 11, 5, 6],
            "centriole_id": [1, 1, 2, 2, 3, 3],
            "centrosome_id": [10, 10, 10, 10, 20, 20],
        })
        centrosome_ids, lines, colors = FindCentrosomesOperator.as_lines(df, track_colors={10: "#ff0000"})

        assert set(centrosome_ids) == {10, 20}

        idx_10 = centrosome_ids.index(10)
        idx_20 = centrosome_ids.index(20)

        # centrosome 10 a 2 centrioles -> 2 composantes par point temporel
        assert lines[idx_10][0].shape == (2, 3)
        assert all(c == "#ff0000" for c in colors[idx_10])

        # centrosome 20 n'a pas de couleur définie -> fallback blanc
        assert all(c == "#ffffff" for c in colors[idx_20])


# ----------------------------------------------------------------------
# run() -- test d'intégration léger
# ----------------------------------------------------------------------

class TestRunIntegration:

    def test_raises_without_image(self, operator):
        with pytest.raises(ValueError, match="Input image must be set"):
            operator.run()

    def test_raises_without_hints(self, calibrated_operator):
        with pytest.raises(ValueError, match="Hint points must be set"):
            calibrated_operator.run()

    @staticmethod
    def get_testing_image() -> np.ndarray:
        cache_dir = Path(__file__).parent / "data"
        cache_dir.mkdir(exist_ok=True)
        image_path = cache_dir / "testing.tif"

        if not image_path.exists():
            urllib.request.urlretrieve(TESTING_IMAGE_URL, image_path)

        return tiff.imread(image_path)

    def test_run_smoke_test(self):
        """
        Test d'intégration léger : deux blobs qui bougent sur quelques frames.
        Vérifie que le pipeline complet tourne sans erreur et produit un
        résultat cohérent -- pas un test de précision sur les positions
        exactes, les paramètres par défaut n'étant pas garantis optimaux
        pour cette image synthétique.
        """
        img = self.get_testing_image()

        op = FindCentrosomesOperator()
        op.set_input_image(
            img,
            calibration={"T": 1.0, "Y": 0.103, "X": 0.103},
            units={"T": "s", "Y": "µm", "X": "µm"},
        )
        op.set_max_binding_distance(5.0)
        op.set_searching_range(5.0)
        op.set_hints({
            1: {"start": 0, "end": img.shape[0] - 1, "points": np.array([[56.0, 33.0], [43.0, 72.0]])},
        })

        op.run()
        result = op.get_centrosomes()

        assert not result.empty
        assert set(result.columns) == {"centriole_id", "T", "Y", "X", "centrosome_id"}
        assert set(result["centrosome_id"]) == {1}