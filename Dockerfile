FROM continuumio/miniconda3

RUN apt-get update \
    && apt-get install -y

COPY environment.yml .

RUN conda env create --file environment.yml \
    && conda clean -afy

SHELL ["conda", "run", "-n", "pyscreener", "/bin/bash", "-c"]

RUN pip install --no-input --no-cache-dir pyscreener

# install vina here
RUN mkdir vina_download \
    && cd vina_download \
    && wget https://vina.scripps.edu/wp-content/uploads/sites/55/2020/12/autodock_vina_1_1_2_linux_x86.tgz \
    && tar -xzvf autodock_vina_1_1_2_linux_x86.tgz \
    && mv autodock_vina_1_1_2_linux_x86/bin/* ../bin/ \
    && cd ..\
    && rm -rf vina_download
