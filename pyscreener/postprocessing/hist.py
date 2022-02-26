from math import ceil
import os
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

sns.set_theme(style="white", context="paper")

BINWIDTH = 0.1


def histogram(hist_mode: str, Y: np.ndarray, path: str = ".", name: str = "distribution.png"):
    if hist_mode == "image":
        plot_hist(Y, path, name)
    elif hist_mode == "text":
        print_hist(Y)

    return


def plot_hist(Y: np.ndarray, path: str = ".", name: str = "scores_distribution.png"):
    bins = np.arange(Y.min(), Y.max() + BINWIDTH, BINWIDTH)

    fig, axs = plt.subplots(1, 2, sharex=True, figsize=(10, 4))

    for ax in axs:
        ax.hist(Y, color="b", edgecolor="none", bins=bins)
        ax.grid(True, linewidth=1, color="whitesmoke")
    axs[0].set_ylabel("Count")
    axs[1].set_yscale("log")

    ax = fig.add_subplot(111, frameon=False)
    ax.tick_params(labelcolor="none", top=False, bottom=False, left=False, right=False)
    ax.set_xlabel("Score")

    fig.tight_layout()
    filepath = Path(path) / name

    fig.savefig(str(filepath), dpi=300)
    print(f"Histogram saved to {filepath}!")


def print_hist(Y: np.ndarray):
    hist, bin_edges = np.histogram(
        Y[~np.isnan(Y)], np.arange(np.nanmin(Y), np.nanmax(Y) + BINWIDTH, BINWIDTH)
    )

    try:
        terminal_size = os.get_terminal_size()[0]
    except OSError:
        terminal_size = 80
    base_width = max(len(str(bin_edges.max().round(1))), len(str(bin_edges.min().round(1))))
    width = terminal_size - (base_width + 2)

    counts = hist / hist.max()
    unit = ceil(hist.max() / width)

    header_str = f"| Score Distribution (* = {unit} counts) |"
    border_str = f"+{'-'*(len(header_str)-2)}+"

    print(f"{border_str: ^0{terminal_size}}")
    print(f"{header_str: ^0{terminal_size}}")
    print(f"{border_str: ^0{terminal_size}}")
    print("-" * terminal_size)
    for count, edge in zip(counts, bin_edges):
        bar = "*" * ceil(width * count)
        print(f"{edge:>{base_width}.1f}: {bar}")
    print("-" * terminal_size)
