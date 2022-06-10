import configparser
import pandas as pd
from .exceptions import MetadataError


# ================================================================
# Parse parameters into an BasicArgument class
#
# ================================================================


class BasicArguments:
    """
    Class to hold arguments

    """

    def __init__(self, expt_dir, config, barcode):
        self.expt_dir = expt_dir
        self.config = config
        self.barcode = barcode


def build_parameter_dict(expt_dir, config, barcode):
    """
    Build a parameter dictionary for NOMADIC, based on the commad line inputs

    """
    # Parse core arguments
    args = BasicArguments(expt_dir, config, barcode)

    # Add derived arguments, order matters here
    args = add_config(args)
    args = add_directories(args)
    args = add_metadata(args)

    # Format barcode parsing
    if args.barcode:
        args.focus_barcode = f"barcode{args.barcode:02d}"

    return args.__dict__


# ================================================================
# Add specific groups of arguments to arguments object
#
# ================================================================


def add_config(args):
    """
    Parse the configuration file

    I think this function, and perhaps even the configuration file,
    can be largely deprectiated.

    Could be kept, for example, if the goal is to 'runall';
    then all arguments are housed in an `.ini` file.

    But for specific scripts, probably want to allow more flexible argument
    passing

    """

    # Parse config
    config = configparser.ConfigParser()
    config.read(args.config)

    # [Experiment]
    args.metadata = config.get("Experiment", "metadata")
    args.basecalling = config.get("Experiment", "basecalling")

    # [Genes]
    args.target_ids = [g.strip() for g in config.get("Genes", "target_ids").split(",")]
    args.target_names = [
        g.strip() for g in config.get("Genes", "target_names").split(",")
    ]
    args.name_dt = {t: n for t, n in zip(args.target_ids, args.target_names)}

    # [Mutations]
    if config.has_section("Mutations"):
        mutations = config.get("Mutations", "csv")
        args.mutations = pd.read_csv(mutations)
        args.mutation_dt = {
            target: gdf["mutation"].tolist()
            for target, gdf in args.mutations.groupby("target")
        }
    else:
        args.mutation_dt = {}

    return args


def add_directories(args):
    """Add experiment directories to arguments"""

    # Define key directories
    args.fastq_dir = f"{args.expt_dir}/{args.basecalling}"
    if args.basecalling == "minknow":
        args.fastq_dir += "/fastq_pass"
    args.nomadic_dir = f"{args.expt_dir}/nomadic/{args.basecalling}"
    args.barcodes_dir = f"{args.nomadic_dir}/barcodes"

    return args


def add_metadata(args, include_unclassified=False):
    """Get metadata"""

    # Load metadata
    metadata_path = f"{args.expt_dir}/{args.metadata}"
    args.metadata = pd.read_csv(metadata_path)
    args.barcodes = args.metadata.barcode.tolist()
    if include_unclassified:
        args.barcodes += ["unclassified"]

    # Sanity checks
    required_columns = ["sample_id", "barcode"]
    for rc in required_columns:
        if not rc in args.metadata.columns:
            raise MetadataError(f"Metadata file {metadata_path} must have a {rc} column.")
        if not len(args.metadata[rc]) == len(args.metadata[rc].unique()):
            raise MetadataError(f"All entries in {rc} column must be unique.")

    return args
