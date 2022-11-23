import os
import numpy as np
import pandas as pd


# --------------------------------------------------------------------------------
# Parameters
#
# --------------------------------------------------------------------------------


EXPT_DIRS = [
    "experiments/2021-11-14_strain-validation-flongle-lfb",
    "experiments/2021-11-20_strain-validation-16plex",
    "experiments/2022-03-15_large-lab-validation-r1",
]
OUTPUT_PATH = "scripts/tables/fig3.expt_aggregate.csv"


# --------------------------------------------------------------------------------
# Experiment Class
#
# --------------------------------------------------------------------------------


class Experiment:

    mapping_levels = ["pf_mapped", "hs_mapped", "unmapped"]
    thresh = 100
    focus_columns = [
        "barcode",
        "sample_id",
        "parasites_per_ul",
        "parasites_per_ul_log10",
        "per_pf_mapped",
        "per_hs_mapped",
        "n_reads",
        "n_reads_log10",
        "imb_factor",
        "imb_factor_log10",
        "n_reads_ontarget_total",
        "n_reads_ontarget_total_log10",
        "n_reads_ontarget_mean",
        "n_reads_ontarget_mean_log10",
        "n_targets",
        f"n_targets_gr{thresh}",
    ]

    def __init__(self, expt_dir):
        """
        Co-ordinate metadata, mapping and balance data for a single
        NOMADS nanopore experiment

        """
        # Define directories
        self.expt_dir = expt_dir
        self.nomadic_dir = f"{expt_dir}/nomadic/guppy/hac/single_end/"

        # Define paths
        self.metadata_path = f"{expt_dir}/metadata/sample_info.csv"
        self.mapping_path = (
            f"{self.nomadic_dir}/qc-bams/table.mapping.primary_state.csv"
        )
        self.extraction_path = (
            f"{self.nomadic_dir}/target-extraction/table.target_coverage.overview.csv"
        )

    def _load_metadata(self):
        """Load metadata"""

        self.metadata = pd.read_csv(self.metadata_path)

    def _load_mapping(self):
        """Load mapping data"""

        # LOAD
        self.mapping_df = pd.read_csv(self.mapping_path)

        # SUMMARY STATISTICS
        # No. reads
        self.mapping_df.insert(
            4, "n_reads", self.mapping_df[self.mapping_levels].sum(1)
        )
        self.mapping_df.insert(
            5, "n_reads_log10", np.log10(self.mapping_df["n_reads"] + 1)
        )

        # Percentages
        for i, mapping_level in enumerate(self.mapping_levels):
            self.mapping_df.insert(
                6 + i,
                f"per_{mapping_level}",
                100 * self.mapping_df[mapping_level] / self.mapping_df["n_reads"],
            )
        self.mapping_df.insert(
            self.mapping_df.shape[1],
            "parasites_per_ul_log10",
            np.log10(self.mapping_df["parasites_per_ul"]),
        )

    def _load_balance(self):
        """
        Load balance data

        TODO:
        - Compute median coverage?

        """

        # LOAD
        self.extraction_df = pd.read_csv(self.extraction_path)

        # RESHAPE
        self.balance_df = pd.pivot_table(
            data=self.extraction_df.query("overlap == 'any'"),
            values="reads_mapped",
            index="sample_id",
            columns="gene_name",
        )

        # SUMMARY STATISTICS
        # Ontarget
        n_ontarget_total = self.balance_df.sum(1)
        n_ontarget_mean = self.balance_df.mean(1)

        # Imbalance factor
        min_ontarget = self.balance_df.min(1)
        max_ontarget = self.balance_df.max(1)
        imb_factor = max_ontarget / (min_ontarget + 1)

        # Above threshold
        n_targets = self.balance_df.shape[1]
        targets_abv_thresh = self.balance_df.apply(
            lambda x: (x > self.thresh).sum(), axis=1
        )

        # Store
        self.balance_df.insert(
            self.balance_df.shape[1], "n_reads_ontarget_total", n_ontarget_total
        )
        self.balance_df.insert(
            self.balance_df.shape[1],
            "n_reads_ontarget_total_log10",
            np.log10(n_ontarget_total),
        )
        self.balance_df.insert(
            self.balance_df.shape[1], "n_reads_ontarget_mean", n_ontarget_mean
        )
        self.balance_df.insert(
            self.balance_df.shape[1],
            "n_reads_ontarget_mean_log10",
            np.log10(n_ontarget_mean),
        )
        self.balance_df.insert(self.balance_df.shape[1], "imb_factor", imb_factor)
        self.balance_df.insert(
            self.balance_df.shape[1], "imb_factor_log10", np.log10(imb_factor)
        )
        self.balance_df.insert(self.balance_df.shape[1], "n_targets", n_targets)
        self.balance_df.insert(
            self.balance_df.shape[1], f"n_targets_gr{self.thresh}", targets_abv_thresh
        )

    def _merge_dfs(self):
        """
        Merge metadata, mapping data, and balance data

        """
        n_barcodes = self.metadata.shape[0]
        assert self.balance_df.shape[0] == n_barcodes
        assert self.mapping_df.shape[0] == n_barcodes

        self.merged_df = pd.merge(
            left=self.mapping_df, right=self.balance_df, on="sample_id", how="inner"
        )

        assert self.merged_df.shape[0] == n_barcodes

    def load_data(self):
        """Load all data and merge"""

        self._load_metadata()
        self._load_mapping()
        self._load_balance()
        self._merge_dfs()

    def get_focus_df(self):
        """Return dataframe with only column of interest"""

        return self.merged_df[self.focus_columns]


# --------------------------------------------------------------------------------
# Main Script
#
# --------------------------------------------------------------------------------


def main():
    """
    Load data for `EXPT_DIRS` and write as a single
    tall data frame

    """

    # Load summary data for each experiment
    dfs = []
    for j, expt_dir in enumerate(EXPT_DIRS):

        # Instantiate experiment class and load data
        print(f"Loading experiment: {expt_dir}")
        assert os.path.exists(expt_dir), "Experiment directory does not exist."
        expt = Experiment(expt_dir)
        expt.load_data()

        # Store key data for experiment
        df = expt.get_focus_df()
        print(f"  Discovered {df.shape[0]} samples for this experiment.")
        df.insert(0, "expt_name", os.path.basename(expt_dir))
        df.insert(0, "expt_ix", j)
        dfs.append(df)
        print("Done.\n")

    # Aggregate and write
    print("Aggregating all experiments...")
    summary_df = pd.concat(dfs)
    print(f"  Total of {summary_df.shape[0]} samples found.")
    print(f"  Storing {summary_df.shape[1]} features for each sample.")
    summary_df.to_csv(OUTPUT_PATH, index=False)
    print(f"  Wrote to: {OUTPUT_PATH}")
    print("Done.\n")


if __name__ == "__main__":
    main()
