import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec


class MutationPanelPlot:

    # Plotting aesthetics
    edge_buffer = 0.75
    ec = "black"
    lw = 0.5

    # Point sizes
    pass_size = 60
    fail_size = 20

    # Scaling
    size_scaling = 0.25
    gap_scaling = 0.2

    def __init__(self, input_df):
        """
        Create a plot of a panel of mutations, optionally focussing on
        a single gene

        params:
            input_df: Pandas DataFrame
                Input data frame where each row is a mutation X sample_id.
                Requires the following columns...
                    gene_name -> Name of genes
                    mutation -> Mutation, e.g. K76T

        returns:
            None

        """

        self.df = self._prepare_dataframe(input_df)
        self.n_mutations = self._count_unique(self.df, "mutation")
        self.n_samples = self._count_unique(self.df, "sample_id")

    @staticmethod
    def _prepare_dataframe(df):
        """
        Prepare the provided dataframe for plotting

        """

        pdf = df.copy()

        required_columns = ["sample_id", "gene_name", "mutation"]
        for r in required_columns:
            assert r in pdf.columns, f"Dataframe must have {r} as column."

        if not "position" in pdf.columns:
            pdf["position"] = [int(m[1:-1]) for m in pdf["mutation"]]

        return pdf

    @staticmethod
    def _count_unique(df, col):
        """Count number of unique elements in a specific dataframe column"""

        return len(df[col].unique())

    def _get_gene_data(self, focus_gene, sort_mutations=True):
        """Get data for a specific gene"""

        gene_df = self.df.query("gene_name == @focus_gene")

        if sort_mutations:
            gene_df = gene_df.sort_values(["position", "mutation"])

        return gene_df

    def plot_gene(
        self, focus_gene, focus_stat="gt", cmap="Reds", ax=None, output_path=None
    ):
        """Plot a specific gene"""

        gene_df = self._get_gene_data(focus_gene)
        n_gene_mutations = self._count_unique(gene_df, "mutation")

        if ax is None:
            fig, ax = plt.subplots(
                figsize=(
                    self.size_scaling * n_gene_mutations,
                    self.size_scaling * self.n_samples,
                )
            )

        # Plot
        sc = ax.scatter(
            x=gene_df["mutation"],
            y=gene_df["sample_id"],
            c=gene_df[focus_stat],
            s=self.pass_size,
            cmap=cmap,
            ec=self.ec,
            lw=self.lw,
            clip_on=False,
            vmin=0,
            vmax=2,
        )

        # Axes and limits
        ax.set_ylim(-self.edge_buffer, self.n_samples + self.edge_buffer - 1)
        ax.set_xlim(-self.edge_buffer, n_gene_mutations + self.edge_buffer - 1)
        ax.invert_yaxis()

        # Ticks
        ax.xaxis.set_ticks_position("top")
        ax.tick_params(left=False, top=False)
        ax.tick_params(axis="x", rotation=90)
        ax.set_facecolor("None")

        # Labels
        ax.set_title(focus_gene, loc="left", rotation=90 if n_gene_mutations < 3 else 0)

        # Box
        for side in ["top", "right", "bottom", "left"]:
            ax.spines[side].set_visible(False)

        # Optionally save
        if output_path is not None:
            fig.savefig(output_path, bbox_inches="tight", pad_inches=0.5)

        return sc

    def plot_multiple_genes(
        self,
        focus_genes=["CRT1", "DHFR", "DHPS", "K13", "MDR1"],
        cmaps=["Reds"] * 5,
        output_path=None,
    ):
        """Plot multiple genes in a sigle plot"""

        # Compute statistics required to set up plot geometry
        n_gene_mutations = [
            self._count_unique(self._get_gene_data(focus_gene), "mutation")
            for focus_gene in focus_genes
        ]

        # Compute edges
        edges = []
        ix = 0
        for m in n_gene_mutations:
            edges.append((ix, ix + m))
            ix += m

        # Adjust spacing
        xsize = self.n_mutations * self.size_scaling
        ysize = (self.n_samples + 1) * self.size_scaling

        # Plot
        fig, ax = plt.subplots(figsize=(xsize, ysize), sharey=True)
        fig.subplots_adjust(wspace=self.gap_scaling)
        gs = GridSpec(nrows=self.n_samples + 1, ncols=self.n_mutations)

        axes = []
        for s, e in edges:
            axes.append(plt.subplot(gs[:-1, s:e]))

        for ax, focus_gene, cmap in zip(axes, focus_genes, cmaps):
            sc = self.plot_gene(focus_gene, cmap="Reds", ax=ax)
            ax.label_outer()
            ax.set_axisbelow(True)

        # Optionally save
        if output_path is not None:
            fig.savefig(output_path, bbox_inches="tight", pad_inches=0.5)

        return None
