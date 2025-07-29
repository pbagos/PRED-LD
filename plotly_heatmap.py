import pandas as pd
import plotly.graph_objects as go

# 1. Specify which chromosome to plot
chromosome = '1'  # ← set this to the chromosome you want, e.g. '2', 'X', etc.

# 2. Read in the full LD file, making sure CHR comes in as a string
df = pd.read_csv(
    'LD_info_chr_all.txt',
    sep='\t',
    dtype={'CHR': str, 'SNP_A': str, 'SNP_B': str}
)

# 3. Filter for the desired chromosome
df_chr = df[df['CHR'] == chromosome].copy()

# 4. Pivot into a square matrix of R² values
pivot = df_chr.pivot(index='SNP_A', columns='SNP_B', values='R2')

# 5. Prepare labels and matrix
x_labels = pivot.columns.astype(str)
y_labels = pivot.index.astype(str)
z_matrix = pivot.values

# 6. Define available colorscales
colorscales = ['Viridis', 'Plasma', 'Cividis', 'Inferno', 'Magma']

# 7. Build the heatmap trace
heatmap = go.Heatmap(
    z=z_matrix,
    x=x_labels,
    y=y_labels,
    colorscale=colorscales[0],
    colorbar=dict(title='R²')
)

fig = go.Figure(data=[heatmap])

# 8. Add a dropdown so you can switch colorscales on the fly
buttons = [
    dict(
        method='restyle',
        label=scale,
        args=[{'colorscale': scale}, [0]]
    ) for scale in colorscales
]

# 9. Final layout tweaks
fig.update_layout(
    title=f'LD Heatmap (R²) for Chromosome {chromosome}',
    xaxis_title='SNP_B',
    yaxis_title='SNP_A',
    updatemenus=[dict(
        active=0,
        buttons=buttons,
        x=1.15,
        y=1.1,
        xanchor='left',
        yanchor='top',
        direction='down',
        pad={'r': 10, 't': 10}
    )],
    margin=dict(l=80, r=200, t=80, b=80)
)

# 10. Show it!
fig.show()
