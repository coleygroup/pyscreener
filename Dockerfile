# ------------------------------------------------------------------------------------------------------------
FROM continuumio/miniconda3 AS base

RUN apt-get update \
    && apt-get install make g++ libboost-all-dev xutils-dev -y

COPY environment.yml .

RUN conda env create --file environment.yml \
    && conda clean -afy

SHELL ["conda", "run", "-n", "pyscreener", "/bin/bash", "-c"]

RUN pip install --no-input --no-cache-dir pyscreener


# ------------------------------------------------------------------------------------------------------------
FROM base AS vina

RUN mkdir vina_download \
    && cd vina_download \
    && wget https://vina.scripps.edu/wp-content/uploads/sites/55/2020/12/autodock_vina_1_1_2_linux_x86.tgz \
    && tar -xzvf autodock_vina_1_1_2_linux_x86.tgz \
    && mv autodock_vina_1_1_2_linux_x86/bin/* ../bin/ \
    && cd ..\
    && rm -rf vina_download


# ------------------------------------------------------------------------------------------------------------
FROM base AS psovina

RUN mkdir psovina_download \
    && cd psovina_download \
    && wget -O psovina-2.0.tar.gz https://sourceforge.net/projects/psovina/files/psovina-2.0.tar.gz/download \
    && tar -xzvf psovina-2.0.tar.gz \
    && cd psovina-2.0/build/linux/release/ \
    && make \
    && mv psovina* ../../../../../bin/ \
    && cd ../../../../../ \
    && rm -rf psovina_download


# ------------------------------------------------------------------------------------------------------------
FROM base AS smina

RUN mkdir smina_download \
    && cd smina_download \  
    && wget -O smina https://sourceforge.net/projects/smina/files/smina.static/download \
    && chmod +x smina \
    && mv smina ../bin/ \
    && cd ../ \ 
    && rm -rf smina_download


# ------------------------------------------------------------------------------------------------------------
FROM base AS qvina

SHELL ["/bin/bash", "-c"]

RUN git clone -b qvina2_1buffer --single-branch https://github.com/QVina/qvina.git \
    && cd qvina \
    && BOOST_LOC=$(whereis boost) && BOOST_PATH=${BOOST_LOC:7:19} && BOOST_VERSION=$(grep "#define BOOST_LIB_VERSION" ../../usr/include/boost/version.hpp | grep -o '".*"' | sed 's/"//g') \
    && sed -i "1s|.*|BASE=$BOOST_PATH|" Makefile && sed -i "2s|.*|BASE=$BOOST_VERSION|" Makefile \
    && make \
    && mv qvina02 ../bin/qvina \
    && cd ../ && rm -rf qvina
    # rename qvina02 to qvina otherwise pyscreener raises an error with executable name


# ------------------------------------------------------------------------------------------------------------
