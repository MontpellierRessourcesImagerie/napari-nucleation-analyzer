import numpy as np
import pandas as pd

class SummaryOperator:

    def __init__(self):
        self.kymographs  = None
        self.summary     = None
        self.spots       = None
        self.centrosomes = None

    def set_kymographs(self, kymographs):
        self.kymographs = kymographs

    def set_spots(self, spots):
        self.spots = spots

    def set_centrosomes(self, centrosomes):
        required_columns = {"centriole_id", "T", "Y", "X", "centrosome_id"}
        if not required_columns.issubset(centrosomes.columns):
            raise ValueError(f"Centrosomes dataframe must contain columns: {required_columns}")
        self.centrosomes = centrosomes

    def get_centrosomes(self) -> pd.DataFrame:
        if self.centrosomes is None:
            raise ValueError("Centrosomes dataframe has not been set.")
        return self.centrosomes

    def _build_summary(self, all_kymographs, all_spots):
        """
        Counts the number of spots per frame.
        Time is the vertical axis.
        """
        blocks = []
        centrosomes = self.get_centrosomes()
        pairs = self._build_pairs(centrosomes)

        for centriole_id, kymo in all_kymographs.items():
            t_start = centrosomes[centrosomes['centriole_id'] == centriole_id]['T'].min()
            spots = all_spots[centriole_id]
            canvas = np.zeros_like(kymo, dtype=np.uint8)
            canvas[spots[:, 0], spots[:, 1]] = 1
            counts = np.sum(canvas, axis=1)
            identifier = f"C{int(pairs[centriole_id])} → c{int(centriole_id)}"

            block = pd.DataFrame({
                f"Count [{identifier}]": counts,
                f"T [{identifier}]": np.arange(len(counts)) + t_start,
            })
            blocks.append(block)

        return pd.concat(blocks, axis=1)

    def _build_pairs(self, centrosomes_df):
        unique_centrioles = centrosomes_df['centriole_id'].unique()
        pairs = {}
        for centriole_id in unique_centrioles:
            centriole_data = centrosomes_df[centrosomes_df['centriole_id'] == centriole_id]
            unique_centrosomes = centriole_data['centrosome_id'].unique()
            centrosomes = unique_centrosomes.tolist()
            if len(centrosomes) != 1:
                raise ValueError(f"Centriole {centriole_id} has multiple centrosomes: {centrosomes}")
            pairs[centriole_id] = centrosomes[0]
        return pairs
    
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

    # Loading the centrosomes dataframe
    centrosomes_path = Path("/home/clement/Documents/projects/nucleation/draft/implementation/dump/centrosomes_arcs.csv")
    centrosomes = pd.read_csv(centrosomes_path)

    # Running the summary operator
    summary_operator = SummaryOperator()
    summary_operator.set_kymographs(kymographs)
    summary_operator.set_spots(spots)
    summary_operator.set_centrosomes(centrosomes)
    summary_operator.run()

    summary_df = summary_operator.get_summary()
    print(summary_df)