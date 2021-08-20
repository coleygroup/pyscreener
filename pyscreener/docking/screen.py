from pyscreener.base import VirtualScreen

import ray

class DockingVirtualScreen(VirtualScreen):
    def __init__(self) -> None:
        super().__init__()

        if not ray.is_initialized():
            try:
                ray.init('auto')
            except ConnectionError:
                ray.init()
        
    def prepare(self):
        pass

    def run(self):
        pass
    