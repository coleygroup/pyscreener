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
FROM base as base-dock6

RUN apt-get install libxss1 xscreensaver xscreensaver-gl-extra xvfb -y \
    && mkdir dock6_utils

RUN cd dock6_utils \
    && mkdir sphgen_cpp_download && cd sphgen_cpp_download \
    && wget http://dock.compbio.ucsf.edu/Contributed_Code/code/sphgen_cpp.1.2.tar.gz \
    && tar -xzvf sphgen_cpp.1.2.tar.gz \
    && cd sphgen_cpp.1.2 \
    && make \
    && mv sphgen_cpp ../../ \
    && cd ../../ && rm -rf sphgen_cpp_download

RUN cd dock6_utils \
    && wget -O chimera-installer.bin "https://www.cgl.ucsf.edu/chimera/cgi-bin/secure/chimera-get.py?ident=OHeQer2WSaxn%2FepyqXlK%2BP9gulBVQ9j%2B1xtw0AnmnvIv&file=linux_x86_64%2Fchimera-1.16-linux_x86_64.bin&choice=Notified" \
    && chmod +x chimera-installer.bin


# ------------------------------------------------------------------------------------------------------------
FROM base AS base-vina

RUN wget -O ADFRsuite.tar.gz https://ccsb.scripps.edu/adfr/download/1038/ \
    && tar -xzvf ADFRsuite.tar.gz \
    && cd ADFRsuite_* \
    && echo "Y" | ./install.sh -d . -c 0 \
    && cd ..\
    && rm ADFRsuite.tar.gz

ENV PATH="${PATH}:/ADFRsuite_x86_64Linux_1.0/bin:"


# ------------------------------------------------------------------------------------------------------------
FROM base-vina AS vina

RUN mkdir vina_download \
    && cd vina_download \
    && wget https://vina.scripps.edu/wp-content/uploads/sites/55/2020/12/autodock_vina_1_1_2_linux_x86.tgz \
    && tar -xzvf autodock_vina_1_1_2_linux_x86.tgz \
    && mv autodock_vina_1_1_2_linux_x86/bin/* ../bin/ \
    && cd ..\
    && rm -rf vina_download


# ------------------------------------------------------------------------------------------------------------
FROM base-vina AS psovina

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
FROM base-vina AS smina

RUN mkdir smina_download \
    && cd smina_download \  
    && wget -O smina https://sourceforge.net/projects/smina/files/smina.static/download \
    && chmod +x smina \
    && mv smina ../bin/ \
    && cd ../ \ 
    && rm -rf smina_download


# ------------------------------------------------------------------------------------------------------------
FROM base-vina AS qvina

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
