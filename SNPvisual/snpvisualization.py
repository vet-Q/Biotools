import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import matplotlib.colors as mcl
import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap

df = pd.read_excel(r"C:\Users\user\Desktop\230705_ASFV APQA.alignment.variation2.xlsx",header=1)

seq = ['1st', '2nd', '3rd', '4th', '5th', '6th', '7th', '8th', '9th', '10th', '11th', '12th', '13th', '14th', '15th','16th', '17th', '18th', '19th', '20th', '21st']

for i in df.index:
    for j in seq:
        if df.loc[i, "NC_044959"] == df.loc[i, j]:
            df.loc[i, j] = 0
        elif df.loc[i,j] == '-':
            df.loc[i, j] = 7
        else: df.loc[i, j] = 12

h = 24
s = 0.99
v = 1

colors = [
    mcl.hsv_to_rgb((h / 360, 1, v)),
    mcl.hsv_to_rgb((h / 360, 0.5, v)),
    mcl.hsv_to_rgb((h / 360, 0, v)),
]
cmap = LinearSegmentedColormap.from_list('my_cmap', colors, gamma=3)

df_sorted = df.drop(['Type', 'Variation', 'CDS gene', 'Mutation','NC_044959'],axis='columns')
df_sorted.index = df_sorted['Position']
df_sorted.rename(columns={'Position':'NC_044959'},inplace=True)
df_sorted['NC_044959'] = 0
df_sorted[seq] = df_sorted[seq].astype(float)

# Build figure and axes
fig, axs = plt.subplots(2, 2, sharex="col", sharey="row", figsize=(40,20),
    gridspec_kw=dict(height_ratios=[1,1.5],width_ratios=[1.5, 1]))
axs[0, 1].set_visible(False)
axs[0, 0].set_box_aspect(1/25)
axs[1, 1].set_box_aspect(1/1)


# subplot간의 공간을 없앰
plt.subplots_adjust(wspace=0, hspace=0)


green = sns.light_palette("seagreen", reverse=True, as_cmap=True)
green.set_over('tomato')
sns.set(font_scale=0.7)
sns.color_palette("crest",3)
ax = sns.heatmap(df_sorted, square=False, linewidths=0.6, annot=False, cmap=cmap, cbar=False,
                 annot_kws={'fontsize': 12, 'fontstyle': 'italic', 'color':'b', 'alpha': 0.6,
                       'rotation': 'vertical', 'verticalalignment': 'center', 'backgroundcolor': 'w'},ax=axs[1,0])
ax.tick_params(axis='y', length=0, labelrotation=0, pad=10)
legend_handles = [Patch(color=colors[0], label='REF.'),
                  Patch(color=colors[1], label='ALT.'),  # alt
                  Patch(color=colors[2], label='InDel')]  # Indel

plt.legend(handles=legend_handles, ncol=3, bbox_to_anchor=[0.5, -0.1], loc='lower center', fontsize=8, handlelength=.8)

# Make means by axis
ax0_sum = df_sorted[seq].astype(bool).sum(axis=0)
ax1_sum = df_sorted[seq].astype(bool).sum(axis=1)

# # Rotate the tick labels and set their alignment.
# axs[1, 1].barh(y=ax1_sum.index, width=ax1_sum.values)
sns.kdeplot(y=ax1_sum.values, ax=axs[1,1])
axs[0, 0].bar(x=ax0_sum.index, height=ax0_sum.values)

from matplotlib.transforms import Bbox

(x0m, y0m), (x1m, y1m) = axs[1, 0].get_position().get_points()  # main heatmap
(x0h, y0h), (x1h, y1h) = axs[0, 0].get_position().get_points()  # horizontal histogram
axs[0, 0].set_position(Bbox([[x0m, y0h], [x1m, y1h]]))
(x0v, y0v), (x1v, y1v) = axs[1, 1].get_position().get_points()  # vertical histogram
axs[1, 1].set_position(Bbox([[x0v, y0m], [x1v, y1m]]))

plt.show()
