from typing import List

from .hist import histogram


def postprocess(postprocessing_options: List[str], **kwargs):
    if "none" in postprocessing_options:
        return

    if "hist" in postprocessing_options:
        histogram(**kwargs)
