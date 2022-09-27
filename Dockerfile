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

# need these for psovina install / build
RUN apt-get install make g++ libboost-all-dev -y

RUN mkdir psovina_download \
    && cd psovina_download \
    && wget -O psovina-2.0.tar.gz https://sourceforge.net/projects/psovina/files/psovina-2.0.tar.gz/download \
    && tar -xzvf psovina-2.0.tar.gz \
    && cd psovina-2.0/build/linux/release/ \
    && make \
    && mv psovina* ../../../../../bin/ \
    && cd ../../../../../ \
    && rm -rf psovina_download

# smina
RUN mkdir smina_download \
    && cd smina_download \  
    && wget -O smina https://sourceforge.net/projects/smina/files/smina.static/download \
    && chmod +x smina \
    && mv smina ../bin/ \
    && cd ../ \ 
    && rm -rf smina_download
