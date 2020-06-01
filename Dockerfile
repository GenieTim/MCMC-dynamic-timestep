# Docker image file
# based on https://github.com/edoardob90/mmm2020/blob/master/Dockerfile
# Pulling from jupyter/scipy-notebook
FROM jupyter/base-notebook

USER root

# Install Ubuntu packages
RUN apt-get -y update && apt-get install -yq --no-install-recommends --fix-missing \
  apt-utils \
  build-essential \
  git \
  nano \
  vim \
  clang \
  cmake \
  libboost-test-dev \
  && apt-get clean \
  && rm -rf /var/lib/apt/lists/*

# Swtich back to normal user
USER $NB_UID

WORKDIR /home/${NB_USER}

# Create a directory for the project source
RUN mkdir -p project7-src

# Create a container volume to make data persistent
VOLUME /home/${NB_USER}/project7-src

# Download and compile librascal
RUN git clone https://github.com/cosmo-epfl/librascal .librascal

RUN mkdir -p .librascal/build && \
  cd .librascal/build && \
  cmake -DUSER=ON -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTS=ON -DBUILD_BINDINGS=ON .. && \
  make install

# Install Python requirements with Conda
RUN conda install --quiet --yes \
  'numpy' \
  'scipy=1.4*' \
  'seaborn=0.9.*' \
  'pandas=0.25.*' \
  'matplotlib-base=3.1.*' \
  'scikit-learn=0.22.*'

# Install dependencies from another channel
RUN conda config --add channels omnia --add channels conda-forge \
  && conda install -y -q openforcefield openmm openmmtools

RUN conda install -y -c conda-forge \
  'ipywidgets' \
  'ase' \
  'nglview=2.7.1' \
  'spglib' && \
  conda clean --all -f -y && \
  jupyter nbextension enable --py --sys-prefix nglview && \
  jupyter nbextension enable --py --sys-prefix widgetsnbextension && \
  jupyter labextension install @jupyter-widgets/jupyterlab-manager && \
  jupyter labextension install nglview-js-widgets@2.7.1 && \
  rm -rf $CONDA_DIR/share/jupyter/lab/staging && \
  rm -rf /home/$NB_USER/.cache/yarn && \
  rm -rf /home/$NB_USER/.node-gyp && \
  fix-permissions $CONDA_DIR && \
  fix-permissions /home/$NB_USER

# Switch to normal user
USER $NB_UID
