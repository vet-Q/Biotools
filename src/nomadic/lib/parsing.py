import os
import argparse
import configparser
import pandas as pd


# ================================================================
# Available parsers
#
# ================================================================


def basic_input_parser(cli, descrip):
    """ Get basic set of argumets for NOMADIC """
    
    # Build basic parser
    parser = build_basic_parser(descrip)
    
    # Parse core arguments
    args = parser.parse_args(cli[1:])
    
    # Add derived arguments, order matters here
    args = add_config(args)
    args = add_directories(args)
    args = add_metadata(args)
    
    return args.__dict__


def barcode_input_parser(cli, descrip):
    """ Get arguments for a specific barcode """
    
    # Build basic parser
    parser = build_basic_parser(descrip)
    
    # Add specific arguments
    parser.add_argument("-b", 
                        "--barcode",
                        type=int,
                        help="Barcode on which to focus")
    
    # Parse core arguments
    args = parser.parse_args(cli[1:])
    
    # Add derived arguments, order matters here
    args = add_config(args)
    args = add_directories(args)
    args = add_metadata(args)
    
    # Format barcode parsing
    if args.barcode:
        args.focus_barcode = f"barcode{args.barcode:02d}"
    
    return args.__dict__


# ================================================================
# Components
#
# ================================================================


def build_basic_parser(descrip):
    """ Build the parser object """


    parser = argparse.ArgumentParser(descrip)
    parser.add_argument("-e", 
                        "--expt_dir",
                        help="Path to directory containing experiment data.")
    parser.add_argument("-c", 
                        "--config",
                        default="configs/default.ini",
                        help="Path to configuration file (.ini).") 

    return parser


def add_config(args):
    """ Parse the configuration file """

    # Parse config
    config = configparser.ConfigParser()
    config.read(args.config)
    
    # [Experiment]
    args.metadata = config.get("Experiment", "metadata")
    args.basecalling = config.get("Experiment", "basecalling")

    # [Genes]
    args.target_ids = [g.strip() for g in config.get("Genes", "target_ids").split(",")]
    args.target_names = [g.strip() for g in config.get("Genes", "target_names").split(",")]
    args.name_dt = {t: n for t, n in zip(args.target_ids, args.target_names)}

    # [Mutations]
    if config.has_section("Mutations"):
        mutations = config.get("Mutations", "csv")
        args.mutations = pd.read_csv(mutations)
        args.mutation_dt = {target: gdf["mutation"].tolist() 
                                 for target, gdf in args.mutations.groupby("target")}
    else:
        args.mutation_dt = {}

    # [Files]
    args.gff_path = config.get('Files', 'gff')
    args.fasta_path = config.get('Files', 'fasta')
    
    return args


def add_directories(args):
    """ Add experiment directories to arguments"""
                                                      
    # Define key directories
    args.fastq_dir = f"{args.expt_dir}/{args.basecalling}"
    if args.basecalling == "minknow":
            args.fastq_dir += "/fastq_pass"
    args.nomadic_dir = f"{args.expt_dir}/nomadic/{args.basecalling}"
    args.barcodes_dir = f"{args.nomadic_dir}/barcodes"
    
    return args


def add_metadata(args, include_unclassified=False):
    """ Get metadata """

    # Load metadata
    args.metadata = pd.read_csv(f"{args.expt_dir}/{args.metadata}")
    args.barcodes = args.metadata.barcode.tolist()
    if include_unclassified:
        args.barcodes += ["unclassified"]
        
    return args


# ================================================================
# Old -- too much going on in one function
#
# ================================================================


def parse_inputs(cli, descrip, include_unclassified=False, verbose=True):
    """
    Parse command line inputs
    
    params
        cli : list
            Command line arguments of user, i.e. from sys.argv[1:]
        descrip : str
            Text description of current script.
        include_unclassified : bool
            Add `unclassified` to list of barcodes?
    
    returns
        params : dt
            Dictionary holding information about experiment,
            barcodes, target genes, &c.
            
    """
    
    # CLI
    parser = argparse.ArgumentParser(descrip)
    parser.add_argument("-e", "--expt_dir",
                        help="Path to directory containing experiment data.")
    parser.add_argument("-c", "--config",
                        default="configs/default.ini",
                        help="Path to configuration file (.ini).")
                        
    args = parser.parse_args(cli[1:])
    
    # Configuration
    config = configparser.ConfigParser()
    config.read(args.config)
    params = {}
    params["config_file"] = args.config
    
    # [Experiment]
    metadata = config.get("Experiment", "metadata")
    params["basecalling"] = config.get("Experiment", "basecalling")

    # [Genes]
    params["target_ids"] = [g.strip() for g in config.get("Genes", "target_ids").split(",")]
    params["target_names"] = [g.strip() for g in config.get("Genes", "target_names").split(",")]
    params["name_dt"] = {t: n for t, n in zip(params["target_ids"], params["target_names"])}
    
    # [Mutations]
    if config.has_section("Mutations"):
        mutations = config.get("Mutations", "csv")
        params["mutations"] = pd.read_csv(mutations)
        params["mutation_dt"] = mutation_dt = { target: gdf["mutation"].tolist() 
                                               for target, gdf in params["mutations"].groupby("target") }
    else:
        params["mutation_dt"] = {}

    # [Files]
    params["gff_path"] = config.get('Files', 'gff')
    params["fasta_path"] = config.get('Files', 'fasta')
    
    # Define key directories
    params["expt_dir"] = args.expt_dir
    params["fastq_dir"] = os.path.join(args.expt_dir, params["basecalling"])
    if params["basecalling"] == "minknow":
            params["fastq_dir"] += "/fastq_pass"
    params["nomadic_dir"] = os.path.join(args.expt_dir, "nomadic", params["basecalling"])
    params["barcodes_dir"] = os.path.join(params["nomadic_dir"], "barcodes")

    # Load metadata
    params["metadata"] = pd.read_csv(os.path.join(params["expt_dir"], metadata))
    params["barcodes"] = params["metadata"]["barcode"].tolist()
    if include_unclassified:
        params["barcodes"] += ["unclassified"]
    
    if verbose:
        print("Parsing input...")
        print("  Experiment directory: %s" % params["expt_dir"])
        print("  Configuration file: %s" % args.config)
        print("  Basecalling method: %s" % params["basecalling"])
        print("  Metadata file: %s" % metadata)
        print("  No. barcodes: %d" % len(params["barcodes"]))
        print("  No. target genes: %d" % len(params["target_ids"]))
        print("  NOMADIC directory: %s" % params["nomadic_dir"])
        print("  Barcodes directory: %s" % params["barcodes_dir"])
        if config.has_section("Mutations"):
            print("  Mutation file: %s" % mutations)
        print("Done.")
        print("")
        
    return params

