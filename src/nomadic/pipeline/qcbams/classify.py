import pandas as pd


# ================================================================
# Classify reads based on their alignment state
# - Should be a cleaner way to do below
#
# ================================================================

def convert_to_padded_binary(val):
    return f"{val:012b}"[::-1]

def check_bit(bit_str, pos):
    return bit_str[pos] == '1'

def get_read_mapping_state(flags):
    """
    Classify a read into a given `mapping_state` given the `flags`
    of all aligned segments associated with that read

    params:
        flags: list of ints
            List of FLAG fields associated with aligned segments
            belonging to a single read.

    returns:
        mapping_state: str
            A single keyword summary of the mapping state of the
            read given all of its aligned segments.

    """

    bflags = [convert_to_padded_binary(f) for f in flags]

    # NB: Chimera arbitrarily dominant, given chimera and supp.
    rules = {
        "unmapped": any([check_bit(bf, 2) for bf in bflags]),
        "chim_mapped": any([check_bit(bf, 11) for bf in bflags]),
        "supp_mapped": any([check_bit(bf, 8) for bf in bflags]),
        "uniq_mapped": True
    }

    for mapping_state, answer in rules.items():
        if answer:
            return mapping_state


def get_read_mapping_state_depreciated(flags):
    """
    Classify a read into a given `mapping_state` given the `flags`
    of all aligned segments associated with that read

    params:
        flags: list of ints
            List of FLAG fields associated with aligned segments
            belonging to a single read.

    returns:
        mapping_state: str
            A single keyword summary of the mapping state of the
            read given all of its aligned segments.

    """
    # Probably the best thing to do is to actually convert back to bits;
    # then look for specific flags for each below
    # Or can easily create an unclassified state



    # NB: Chimera arbitrarily dominant, given chimera and supp.
    rules = {
        "uniq_mapped": (0 in flags or 16 in flags) and len(flags) == 1,
        "chim_mapped": (2048 in flags or 2064 in flags) and len(flags) > 1,
        "supp_mapped": (256 in flags or 272 in flags) and len(flags) > 1,
        "unmapped": 4 in flags,
    }

    for mapping_state, answer in rules.items():
        if answer:
            return mapping_state

    raise ValueError(f"Flags {flags} do not conform to any mapping state.")


def get_read_state_summary(species, mapping_state, detailed=False):
    """
    Get a single string summary of the read mapping state, taking into
    account `species` and `mapping_state` information

    params
        species: str
            Species to which the read is mapped.
        mapping_state: str
            Mapping state for the read, i.e. unmapped, uniquely
            mapped, chimera mapped, supplementary mapped, or both.
        detailed: bool
            Return more detailed summary

    return
        _: str
            Summary of the species and mapping state.

    """

    # Handle unmapped
    if mapping_state == "unmapped":
        return mapping_state

    return f"{species}_mapped" if not detailed else f"{species}_{mapping_state}"


def reduce_to_read_dataframe(alignments_df):
    """
    Reduce an alignments dataframe to a read-level dataframe,
    annotating mapping state

    Pretty painful function

    """

    # Annotate every alignment with a mapping state
    mapping_state_dt = {
        read_id: get_read_mapping_state(adf["flag"].tolist())
        for read_id, adf in alignments_df.groupby("read_id")
    }
    alignments_df.insert(
        8,
        "mapping_state",
        [mapping_state_dt[read_id] for read_id in alignments_df["read_id"]],
    )

    # Reduce to primary alignment for every read
    read_df = alignments_df.query("flag in [0, 16, 4]")

    # Annotate with summary of species + mapping state
    read_df.insert(
        9,
        "primary_state",
        [
            get_read_state_summary(s, m)
            for s, m in zip(read_df["species"], read_df["mapping_state"])
        ],
    )
    read_df.insert(
        10,
        "secondary_state",
        [
            get_read_state_summary(s, m, detailed=True)
            for s, m in zip(read_df["species"], read_df["mapping_state"])
        ],
    )

    return read_df


def convert_column_to_ordered_category(df, column, category_order):
    """
    Convert a `column` of a dataframe `df` to an ordered Categorical
    column, given by the order in `category_order`

    """

    # This will occur inplace
    df[column] = pd.Categorical(
        values=df[column], categories=category_order, ordered=True
    )
