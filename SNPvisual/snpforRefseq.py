import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import matplotlib.colors as mcl
import seaborn as sns

import re

def SNPvalues(sheet_name,seq):
    df = pd.read_excel(r"C:\Users\user\Desktop\LAO234_variant table_refseq.xlsx",header=2,sheet_name=sheet_name)

    columns = list(df.columns.values)
    ref = 'Reference\n(NC_044959.1)'
    p = re.compile('Reference\n')
    for i in columns:
        m = p.match(i)
        if m != None:
            ref = i
        else: continue

    # Classify the data to 3 category (REF, ALT, InDel)
    for i in df.index:
        for j in seq:
            if df.loc[i, ref] == df.loc[i, j]:
                df.loc[i, j] = 0
            elif df.loc[i, ref] != df.loc[i,j] and df.loc[i,"Type"] in ["SNV","MNV"]:
                df.loc[i, j] = 7
            else: df.loc[i, j] = 12


    ax0_sum = df[seq].astype(bool).sum(axis=0)
    ax1_sum = df[seq].astype(bool).sum(axis=1)
    ax2_sum = ax1_sum + df['Reference Position']
    conc = pd.concat([df['Reference Position'], ax1_sum], axis=1, keys = ['Reference Position', sheet_name])

    return df['Reference Position'], conc, ax0_sum


NCpos, NC_ax1sum, NC_ax0sum  = SNPvalues("maptoNC",['LAO2 general', 'LAO3 geeneral','LAO4 general'])
MKpos, MK_ax1sum, MK_ax0sum  = SNPvalues("maptoMK",['LAO2 general', 'LAO3 geeneral','LAO4 general'])
OPpos, OP_ax1sum, OP_ax0sum  = SNPvalues("maptoOP",['LAO2 general', 'LAO3 geeneral','LAO4 general'])
MTpos, MT_ax1sum, MT_ax0sum  = SNPvalues("maptoMT",['LAO2 general', 'LAO3 geeneral','LAO4 general'])


# Set the style globally using sns.set_style()
sns.set_style("darkgrid")

# Define a list of color palettes for each subplot
palettes = ['red', 'blue', 'orange', 'grey']

fig, ax = plt.subplots(ncols=2, nrows=2, figsize=(30, 5))

for i in range(2):
    for j in range(2):
        sns.histplot(
            NC_ax1sum if i == 0 and j == 0
            else MK_ax1sum if i == 0 and j == 1
            else OP_ax1sum if i == 1 and j == 0
            else MT_ax1sum,
            x="Reference Position",
            bins=100,
            element="poly",
            kde=True,
            ax=ax[i, j],
            color=palettes[i * 2 + j],  # Use a different palette for each subplot
        )
        ax[i, j].set_title(
            "NC044959" if i == 0 and j == 0
            else "MK_543947.1" if i == 0 and j == 1
            else "OP_467597.1" if i == 1 and j == 0
            else "MT_872723.1",
            fontsize=14,
        )
        ax[i, j].set_xlabel("")  # x축 label 제거
        ax[i, j].set_ylim(0, 30)

plt.tight_layout()
plt.show()