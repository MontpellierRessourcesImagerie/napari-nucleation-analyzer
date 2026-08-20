import numpy as np
import pandas as pd

class SummaryOperator:

    def __init__(self):
        self.kymographs  = None
        self.summary     = None
        self.spots       = None
        self.pairs = None

    def set_kymographs(self, kymographs):
        self.kymographs = kymographs

    def set_spots(self, spots):
        self.spots = spots

    def set_pairs(self, pairs):
        required_columns = {"centrosome_id", "T", "Y", "X", "pair_id"}
        if not required_columns.issubset(pairs.columns):
            raise ValueError(f"Pairs dataframe must contain columns: {required_columns}")
        self.pairs = pairs

    def get_pairs(self) -> pd.DataFrame:
        if self.pairs is None:
            raise ValueError("Pairs dataframe has not been set.")
        return self.pairs

    def _build_summary(self, all_kymographs, all_spots):
        """
        Counts the number of spots per frame.
        Time is the vertical axis.
        """
        blocks = []
        pairs = self.get_pairs()
        matchs = self._build_matchs(pairs)

        for centrosome_id, kymo in all_kymographs.items():
            t_start = pairs[pairs['centrosome_id'] == centrosome_id]['T'].min()
            spots = all_spots[centrosome_id]
            canvas = np.zeros_like(kymo, dtype=np.uint8)
            canvas[spots[:, 0], spots[:, 1]] = 1
            counts = np.sum(canvas, axis=1)
            identifier = f"P{int(matchs[centrosome_id])} → c{int(centrosome_id)}"

            block = pd.DataFrame({
                f"Count [{identifier}]": counts,
                f"T [{identifier}]": np.arange(len(counts)) + t_start,
            })
            blocks.append(block)

        return pd.concat(blocks, axis=1)

    def _build_matchs(self, pairs_df):
        unique_centrosomes = pairs_df['centrosome_id'].unique()
        matchs = {}
        for centrosome_id in unique_centrosomes:
            centrosome_data = pairs_df[pairs_df['centrosome_id'] == centrosome_id]
            unique_pairs = centrosome_data['pair_id'].unique()
            pairs = unique_pairs.tolist()
            if len(pairs) != 1:
                raise ValueError(f"Centrosome {centrosome_id} has multiple pairs: {pairs}")
            matchs[centrosome_id] = pairs[0]
        return matchs
    
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


if __name__ == "__main__":
    from pathlib import Path
    import numpy as np
    import tifffile as tiff

    # Loading the produced kymographs and spots
    kymos_path = Path("/home/clement/Documents/projects/nucleation/draft/implementation/dump/kymos")
    spots_path = kymos_path

    kymographs = {int(f.stem.split('_')[1]): tiff.imread(f) for f in kymos_path.glob("kymo_*.tif")}
    spots = {int(f.stem.split('_')[1]): np.load(f) for f in spots_path.glob("spots_*.npy")}
    print("Loaded kymographs:", list(kymographs.keys()))
    print("Loaded spots:", list(spots.keys()))

    # Loading the pairs dataframe
    pairs_path = Path("/home/clement/Documents/projects/nucleation/draft/implementation/dump/pairs_arcs.csv")
    pairs = pd.read_csv(pairs_path)

    # Running the summary operator
    summary_operator = SummaryOperator()
    summary_operator.set_kymographs(kymographs)
    summary_operator.set_spots(spots)
    summary_operator.set_pairs(pairs)
    summary_operator.run()

    summary_df = summary_operator.get_summary()
    print(summary_df)