from collections import Counter
import os
from pathlib import Path
from typing import Dict, List, Optional, Iterable

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from tqdm import tqdm

sns.set_theme(style='white', context='paper')
BINWIDTH = 0.1

def visualize(viz_mode: str,
              d_smi_score: Dict[str, Optional[float]], 
              d_smi_score_clusters: Optional[List[Dict]] = None,
              name: str = 'distribution', path: str = '.', **kwargs):
    if viz_mode == 'histogram':
        make_hist(ys=d_smi_score.values(), name=name, path=path)

        if d_smi_score_clusters:
            yss = (cluster.values() for cluster in  d_smi_score_clusters)
            make_clustered_hist(yss, name=name, path=path)
    
    if viz_mode == 'text':
        make_text_hist(d_smi_score.values())

    return None

def make_clustered_hist(yss: Iterable[Iterable[float]],
                        name: str = 'distribution', path: str = '.'):
    fig = plt.figure(figsize=(10, 4))

    ax1 = plt.subplot(121)
    ax2 = plt.subplot(122, sharex=ax1)
    for i, ys in enumerate(yss):
        ys = [y for y in ys if y is not None]
        bins = np.arange(min(ys), max(ys)+BINWIDTH, BINWIDTH)
        for ax in (ax1, ax2):
            ax.hist(ys, color='b', edgecolor='none', alpha=0.5,
                    bins=bins, label=f'cluster_{i}')
            ax.set_ylabel('Count')
            ax.grid(True, linewidth=1, color='whitesmoke')
    
    ax2.set_yscale('log')
    ax1.legend()
    ax2.legend()
    # add global x-axis label
    fig.add_subplot(111, frameon=False)
    plt.tick_params(labelcolor='none',
                    top=False, bottom=False, left=False, right=False)
    plt.xlabel('Score')

    plt.tight_layout()
    plt.savefig(str(Path(path)/f'{name}_scores_histogram_clusters.pdf'))
    plt.clf()

def make_hist(ys: Iterable[float], name: str = 'distribution', path: str = '.'):
    ys = [y for y in ys if y is not None]
    bins = np.arange(min(ys), max(ys)+BINWIDTH, BINWIDTH)
    
    fig = plt.figure(figsize=(10, 4))

    ax1 = plt.subplot(121)
    ax2 = plt.subplot(122, sharex=ax1)
    for ax in (ax1, ax2):
        ax.hist(ys, color='b', edgecolor='none', bins=bins)
        ax.set_ylabel('Count')
        ax.grid(True, linewidth=1, color='whitesmoke')
    ax2.set_yscale('log')
    ax2.set_ylabel('Count')
    
    # add global x-axis label
    fig.add_subplot(111, frameon=False)
    plt.tick_params(labelcolor='none',
                    top=False, bottom=False, left=False, right=False)
    plt.xlabel('Score')

    plt.tight_layout()
    plt.savefig(str(Path(path)/f'{name}_scores_histogram_clusters.pdf'))
    plt.clf()

def make_text_hist(ys: Iterable[float]):
    counts = Counter(round(y, 1) for y in ys if y is not None)

    base_width = max(len(str(y)) for y in counts)
    try:
        terminal_size = os.get_terminal_size()[0]
    except OSError:
        terminal_size = 80
    width = terminal_size - (base_width + 2)

    max_counts = counts.most_common(1)[0][1]
    for y, count in counts.items():
        counts[y] = count / max_counts

    for y in np.arange(min(counts), max(counts), 0.1):
        y=round(y, 1)
        bar = '*'*int(width*counts.get(y, 0))
        print(f'{y: >{base_width}}: {bar}')
    print()
    