import numpy as np
import pandas as pd

class SummaryOperator:

    def __init__(self):
        self.kymographs = None
        self.summary    = None
        self.spots      = None

    def set_kymographs(self, kymographs):
        self.kymographs = kymographs

    def set_spots(self, spots):
        self.spots = spots

    def _build_summary(self, all_kymographs, all_spots):
        """
        Counts the number of spots per frame.
        Time is the vertical axis.
        """
        buffer = {}
        for track_id, kymo in all_kymographs.items():
            spots = all_spots[track_id]
            canvas = np.zeros_like(kymo, dtype=np.uint8)
            canvas[spots[:, 0], spots[:, 1]] = 1
            counts = np.sum(canvas, axis=1)
            buffer[f"Track: {track_id}"] = counts
        return pd.DataFrame(buffer)
    
    def get_summary(self):
        if self.summary is None:
            raise ValueError("Summary not computed. Please run the operator first.")
        return self.summary

    def run(self):
        if self.kymographs is None:
            raise ValueError("Kymographs not set. Please set kymographs before running the operator.")
        if self.spots is None:
        
            raise ValueError("Spots not set. Please set spots before running the operator.")
        
        set_spots = set(self.spots.keys())
        set_kymos = set(self.kymographs.keys())
        
        if set_spots != set_kymos:
            raise ValueError(f"Mismatch between spots and kymographs keys. Spots keys: {set_spots}, Kymographs keys: {set_kymos}")
        
        self.summary = self._build_summary(self.kymographs, self.spots)
        