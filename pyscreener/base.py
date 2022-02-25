from abc import ABC
import datetime
from pathlib import Path
import re
import shutil
import tarfile
import tempfile
from typing import Optional, Union

import ray


class VirtualScreen(ABC):
    def __init__(self, path: Union[str, Path]):
        self.path = path
        self.tmp = tempfile.gettempdir()

    @property
    def path(self) -> Path:
        return self.__path

    @path.setter
    def path(self, path):
        path = Path(path)
        path.mkdir(parents=True, exist_ok=True)
        self.__path = path

    @property
    def tmp_dir(self) -> Path:
        """the Screener's temp directory"""
        return self.__tmp_dir

    @tmp_dir.setter
    def tmp_dir(self, path: Union[str, Path]):
        timestamp = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
        path = Path(path) / "pyscreener" / f"session_{timestamp}"
        path.mkdir(parents=True, exist_ok=True)

        self.__tmp_dir = path
        self.tmp_in = path / "inputs"
        self.tmp_out = path / "outputs"

    @property
    def tmp_in(self) -> Path:
        return self.__tmp_in

    @tmp_in.setter
    def tmp_in(self, path: Union[str, Path]):
        path = Path(path)
        path.mkdir(parents=True, exist_ok=True)
        self.__tmp_in = path

    @property
    def tmp_out(self) -> Path:
        return self.__tmp_out

    @tmp_out.setter
    def tmp_out(self, path: Union[str, Path]):
        path = Path(path)
        path.mkdir(parents=True, exist_ok=True)
        self.__tmp_out = path

    def collect_files(self, out_path: Optional[Union[str, Path]] = None):
        """Collect all the files from the local disks of the respective nodes

        For I/O purposes, input and output files for each simulation are
        created on the local disk of each node. If these files are desired at
        the end, they must be copied over from the node's local file system to
        the final destination.

        This is achieved by creating a gzipped tar file of the temp directory
        (the one that contains all of the input and output files for
        simulations conducted on that node) and moving these tar files under
        the desired path. Each tar file is named according the node ID from
        which it originates.

        This function should ideally only be called once during the lifetime
        of a Screener because it is slow and early calls will yield nothing
        over a single, final call.

        Parameters
        ----------
        out_path : Optional[Union[str, Path]], default=None
            the path under which the tar files should be collected to.
            If None, use self.path
        """
        out_path = Path(out_path or self.path)
        out_path.mkdir(parents=True, exist_ok=True)

        refs = []
        for node in ray.nodes():
            address = node["NodeManagerAddress"]

            @ray.remote(resources={f"node:{address}": 0.1})
            def zip_and_move_tmp():
                output_id = re.sub(r"[:,.]", "", ray.state.current_node_id())
                tmp_tar = (self.tmp_dir / output_id).with_suffix(".tar.gz")
                with tarfile.open(tmp_tar, "w:gz") as tar:
                    tar.add(self.tmp_in, arcname="inputs")
                    tar.add(self.tmp_out, arcname="outputs")
                shutil.copy(str(tmp_tar), str(out_path))

            refs.append(zip_and_move_tmp.remote())
        ray.wait(refs)
