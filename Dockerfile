# ------------------------------------------------------------------------------------------------------------
FROM mambaorg/micromamba:latest AS base

USER root

RUN apt-get update \
    && apt-get install wget make g++ libboost-all-dev xutils-dev libxss1 xscreensaver xscreensaver-gl-extra xvfb -y

COPY environment.yml .

ARG MAMBA_DOCKERFILE_ACTIVATE=1

RUN micromamba install -y -n base -f environment.yml \
    && pip install pyscreener \
    && micromamba clean --all --yes


# ------------------------------------------------------------------------------------------------------------
FROM base AS base-vina

RUN wget -O ADFRsuite.tar.gz https://ccsb.scripps.edu/adfr/download/1038/ \
    && mv ADFRsuite.tar.gz ../ && cd .. \
    && tar -xzvf ADFRsuite.tar.gz \
    && cd ADFRsuite_* \
    && echo "Y" | ./install.sh -d . -c 0 \
    && cd .. \
    && rm -rf ADFRsuite.tar.gz

ENV PATH="${PATH}:/ADFRsuite_x86_64Linux_1.0/bin:"


# ------------------------------------------------------------------------------------------------------------
FROM base-vina AS vina

RUN wget https://vina.scripps.edu/wp-content/uploads/sites/55/2020/12/autodock_vina_1_1_2_linux_x86.tgz \
    && tar -xzvf autodock_vina_1_1_2_linux_x86.tgz \
    && mv autodock_vina_1_1_2_linux_x86/bin/* ../bin/


# ------------------------------------------------------------------------------------------------------------
FROM base-vina AS psovina

RUN wget -O psovina-2.0.tar.gz https://sourceforge.net/projects/psovina/files/psovina-2.0.tar.gz/download \
    && tar -xzvf psovina-2.0.tar.gz \
    && cd psovina-2.0/build/linux/release/ \
    && make \
    && mv psovina* ../../../../../bin/


# ------------------------------------------------------------------------------------------------------------
FROM base-vina AS smina

RUN wget -O smina https://sourceforge.net/projects/smina/files/smina.static/download \
    && chmod +x smina \
    && mv smina ../bin/


# ------------------------------------------------------------------------------------------------------------
FROM base-vina AS qvina

SHELL ["/bin/bash", "-c"]

RUN apt-get install git -y \
    && git clone -b qvina2_1buffer --single-branch https://github.com/QVina/qvina.git \
    && cd qvina \
    && BOOST_LOC=$(whereis boost) && BOOST_PATH=${BOOST_LOC:7:19} && BOOST_VERSION=$(grep "#define BOOST_LIB_VERSION" ../../usr/include/boost/version.hpp | grep -o '".*"' | sed 's/"//g') \
    && sed -i "1s|.*|BASE=$BOOST_PATH|" Makefile && sed -i "2s|.*|BASE=$BOOST_VERSION|" Makefile \
    && make \
    && mv qvina02 ../../bin/qvina \
    && cd ../ && rm -rf qvina
    # rename qvina02 to qvina otherwise pyscreener raises an error with executable name


# ------------------------------------------------------------------------------------------------------------
