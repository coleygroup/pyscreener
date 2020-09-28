from collections import Counter
from typing import Dict, List, Optional, Iterable

import numpy as np
import plotly.graph_objects as go

def clustered_hist(yss: Iterable[Iterable[float]], name: str):
    fig = go.Figure()
    for i, ys in enumerate(yss):
        fig.add_trace(go.Histogram(
            x=ys,
            name=f'cluster {i}',
            xbins=dict(size=0.1),
            # opacity=0.5
        ))
    fig.update_layout(
        title_text=f'{name} Score Distribution',
        xaxis_title_text='Docking Score',
        yaxis_title_text='Count',
        barmode='stack'
    )
    fig.write_image(f'{name}_scores_histogram.pdf')

def hist(ys: Iterable[float], name: str):
    fig = go.Figure()
    fig.add_trace(go.Histogram(
        x=ys,
        xbins=dict(size=0.1)
    ))
    fig.update_layout(
        title_text=f'{name} Score Distribution',
        xaxis_title_text='Docking Score',
        yaxis_title_text='Count'
    )
    return fig.write_image(f'{name}_clusters_scores_histogram.pdf')

def text(ys: Iterable[float]):
    counts = Counter(round(y, 1) for y in ys)

    base_width = max(len(str(y)) for y in counts)
    width = 80 - (base_width + 2)   # to account for the ': '

    max_counts = counts.most_common(1)[0][1]
    for y, count in counts.items():
        counts[y] = count / max_counts

    for y in np.arange(min(counts), max(counts), 0.1):
        y=round(y, 1)
        bar = '*'*int(width*counts.get(y, 0))
        print(f'{y: >{base_width}}: {bar}')
    print()

def visualize(viz_mode: str,
              d_smi_score: Dict[str, Optional[float]], name: str,
              d_smi_score_clusters: Optional[List[Dict]] = None, **kwargs):
    if viz_mode == 'histogram':
        hist(ys=d_smi_score.values(), name=kwargs['name'])
        if d_smi_score_clusters:
            yss = (cluster.values() for cluster in  d_smi_score_clusters)
            clustered_hist(yss, name=name)
    
    if viz_mode == 'text':
        pass

    return None