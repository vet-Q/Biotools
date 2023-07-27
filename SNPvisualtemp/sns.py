import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import math

df = pd.read_csv(r"C:\Users\user\Desktop\LAO_ASF\variantTable\LAO2\variant_general.csv")
df_targeted = pd.read_csv(r"C:\Users\user\Desktop\LAO_ASF\variantTable\LAO2\variant.csv")
df['group'] = df['number'].map(lambda x:math.ceil(round(x,-2)))
df['Type'] = "general"
df_targeted['group'] = df_targeted['number'].map(lambda x:math.ceil(round(x,-2)))
df_targeted["Type"] = 'Targeted'

want_mapping= df.groupby('group')['nt_numbers'].sum()
want_mapDf = pd.DataFrame(want_mapping)
want_mapDf['Type'] = "General"
want_mapDf.insert(0,'number',want_mapDf.index)
want_mapDf.reset_index(drop=True,inplace=True)

want_Targetmapping= df_targeted.groupby('group')['nt_numbers'].sum()
want_TargetmapDf = pd.DataFrame(want_Targetmapping)
want_TargetmapDf['Type'] = "Targeted"
want_TargetmapDf.insert(0,'number',want_TargetmapDf.index)
want_TargetmapDf.reset_index(drop=True)

Total = pd.concat([want_mapDf, want_TargetmapDf], axis=0)

sns.set_theme(style="whitegrid")

 # Load the NGS_depth dataset
# planets = sns.load_dataset("planets")

cmap = sns.cubehelix_palette(rot=-.2, as_cmap=True)
g = sns.relplot(
    data=Total,
    x='number', y='nt_numbers',
    size='nt_numbers',
    hue='Type',
    sizes=(10, 200),
)
g.set(xscale="log",yscale="log")
g.ax.xaxis.grid(True, "minor", linewidth=.25)
g.ax.yaxis.grid(True, "minor", linewidth=.25)
g.despine(left=True, bottom=True)
plt.title("LAO2 Depth scatter plot")
plt.show()

