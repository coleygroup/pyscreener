from typing import List

from .cluster import cluster
from .rescoring import rescore
from .visualization import visualize

def postprocess(postprocessing_options: List[str], **kwargs):
    if 'cluster' in postprocessing_options:
        d_smi_score_clusters = cluster(**kwargs)
    else:
        d_smi_score_clusters = None

    if 'visualize' in postprocessing_options:
        visualize(d_smi_score_clusters=d_smi_score_clusters, **kwargs)
