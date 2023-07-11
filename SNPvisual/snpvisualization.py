import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import matplotlib.colors as mcl
import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap


df = pd.read_excel(r"C:\Users\user\PycharmProjects\pythonProject\Biotools\SNPvisual\230711_ASFV APQA.alignment.variation2.xlsx",header=1,sheet_name="ASFV")
labels = pd.read_excel(r"C:\Users\user\PycharmProjects\pythonProject\Biotools\SNPvisual\230711_ASFV APQA.alignment.variation2.xlsx",sheet_name="accessionID")
labels = labels['Sequence ID']
dfAnnotPre = df.copy()
Annotseq = ['NC_044959','1st', '2nd', '3rd', '4th', '5th', '6th', '7th', '8th', '9th', '10th', '11th', '12th', '13th',
       '14th', '15th','16th', '17th', '18th', '19th', '20th', '21st','Korea/YC1/2019','Korea/HC224/2020']
dfForAnnot = dfAnnotPre[Annotseq]
dfForAnnotArray = np.array(dfForAnnot)

seq = ['1st', '2nd', '3rd', '4th', '5th', '6th', '7th', '8th', '9th', '10th', '11th', '12th', '13th',
       '14th', '15th','16th', '17th', '18th', '19th', '20th', '21st','Korea/YC1/2019','Korea/HC224/2020']


# Classify the data to 3 category (REF, ALT, InDel) 
for i in df.index:
    for j in seq:
        if df.loc[i, "NC_044959"] == df.loc[i, j]:
            df.loc[i, j] = 0
        elif df.loc[i, "NC_044959"] not in ["A", "C", "T", "G"] and df.loc[i, "NC_044959"] != df.loc[i, j]:
            if df.loc[i, "NC_044959"] in ["AGG", "CA", "TC"]:
                df.loc[i, j] = 12
            else: df.loc[i, j] = 7
        elif df.loc[i, "NC_044959"] == '-' and df.loc[i,j] != '-':
            df.loc[i, j] = 7
        elif df.loc[i, "NC_044959"] != '-' and df.loc[i, j] == '-':
            df.loc[i, j] = 7
        else: df.loc[i, j] = 12


# SET the color palette
h = 24
s = 0.99
v = 1

colors = [
    mcl.hsv_to_rgb((h / 360, 1, 0.7)),
    mcl.hsv_to_rgb((h / 360, 0.3, v)),
    mcl.hsv_to_rgb((h / 360, 0, v)),
]

# This line makes my own color set. To to this work, we call LinearSegmenetedColormap module.
cmap = LinearSegmentedColormap.from_list('my_cmap', colors, gamma=3)

# Drop the columns which is not necessary to draw graph.
df_sorted = df.drop(['Type', 'Variation', 'CDS gene', 'Mutation','NC_044959'],axis='columns')

# set the Nucleotide number to index (It is helpful to set the y-axis tick and label.)
df_sorted.index = df_sorted['Position']

df_sorted.rename(columns={'Position':'NC_044959'},inplace=True)

# To draw REF heatmap. set the value to 0
df_sorted['NC_044959'] = 0

# heatmap data need to convert np.array (Object type is not acceptable)
df_sorted[seq] = df_sorted[seq].astype(float)


# Build figure and axes
fig, axs = plt.subplots(2, 2, sharex="col", sharey="row", figsize=(30,50),
    gridspec_kw=dict(height_ratios=[1,5],width_ratios=[4, 1],wspace=0, hspace=0))

axs[0, 1].set_visible(False)

axs[0, 0].set_box_aspect(1/11)
axs[1, 1].set_box_aspect(3/1)

# set the palette
green = sns.light_palette("seagreen", reverse=True, as_cmap=True)
green.set_over('tomato')
sns.set(font_scale=0.6)
sns.color_palette("crest",3)

ax = sns.heatmap(df_sorted, square=False, linewidths=0.6, annot=dfForAnnotArray, cmap=cmap, cbar=False,fmt='',
                 annot_kws={'fontsize': 5, 'fontstyle': 'italic', 'color':"black",
                       'rotation':0, 'verticalalignment': 'center'},ax=axs[1,0])

# get ytick label and set the label manually to avoid automatically rotation of labels when adjusting the figure size
ticklabel = []
ticklabel.append("NC_044959")
ticklabel = ticklabel + list(labels) + seq[-2:]
ax.set_yticklabels(ax.get_yticklabels(), rotation=0, va='center')
ax.set_xticklabels(ticklabel, rotation=350, va='top', ha='left')

legend_handles = [Patch(color=colors[0], label='REF.'),
                  Patch(color=colors[1], label='InDel.'),
                  Patch(color=colors[2], label='ALT')]

ax.legend(handles=legend_handles, ncol=3, bbox_to_anchor=[0.5, -0.1], loc='lower center', fontsize=8, handlelength=.8)

# Make sums by axis: astype(bool) method yields boolean result (if a cell value is 0, the result of astype is False.)
ax0_sum = df_sorted[seq].astype(bool).sum(axis=0)
ax1_sum = df_sorted[seq].astype(bool).sum(axis=1)
ax0_mean = ax0_sum.mean()
ax1_mean = ax1_sum.mean()

ax2 = axs[1, 1].barh(y=np.array([i+0.5 for i in range(42)]), width=ax1_sum.values,alpha=0.5,color='dimgrey')
ax1 = axs[0, 0].bar(x=np.array([i+0.5 for i in range(24)]), height=np.insert(ax0_sum.values,0,0), width=0.5, align='center',alpha=0.4,color='dimgrey')
axs[0, 0].axis(xmin=0, xmax=24, ymax=27)
axs[0, 0].set_xticks([i+0.5 for i in range(24)])
axs[0, 0].hlines(y=ax0_mean,xmin=0,xmax=24,linestyles="--",colors="grey",alpha=0.6)
axs[0, 0].text(1, 16,  'mean SNP:%.3f' % ax0_mean, ha='center', va='bottom', size = 10)



## 화살표 작성을 위한 코드
arrow_props = dict(facecolor='black', arrowstyle='->')

import matplotlib.patches as patches

# 화살표를 추가할 x 범위 설정
arrow_x_range = [22, 24]

# y 값 설정
arrow_y = 20

# 양방향 화살표를 그리기 위한 화살표 패치 생성
arrow_patch = patches.FancyArrowPatch(
    (arrow_x_range[0], arrow_y), (arrow_x_range[1], arrow_y),
    arrowstyle='<->', color='blue', mutation_scale=5
)

# axs[0, 0]에 화살표 패치 추가
axs[0, 0].add_patch(arrow_patch)
axs[0, 0].text(23, 20,  "wild boar", ha='center', va='bottom', size = 10, color='black')


# write tick values by for-looping
for rect in ax1:
    height = rect.get_height()
    axs[0,0].text(rect.get_x()+rect.get_width()/2, height, '%.0f' % height, ha='center', va='bottom', size = 10)

for rect in ax2:
    width = rect.get_width()
    axs[1,1].text(2.0,rect.get_y()+rect.get_height(), '%.1f' % width, ha='center', va='bottom', size = 10)

## subplot들의 margin을 일치시키기 위한 코드

from matplotlib.transforms import Bbox

(x0m, y0m), (x1m, y1m) = axs[1, 0].get_position().get_points()  # main heatmap
(x0h, y0h), (x1h, y1h) = axs[0, 0].get_position().get_points()  # horizontal histogram
axs[0, 0].set_position(Bbox([[x0m, y0h], [x1m, y1h]]))
(x0v, y0v), (x1v, y1v) = axs[1, 1].get_position().get_points()  # vertical histogram
axs[1, 1].set_position(Bbox([[x0v, y0m], [x1v, y1m]]))

plt.tight_layout()

# Show graph!
plt.show()

