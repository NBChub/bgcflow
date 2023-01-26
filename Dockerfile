FROM condaforge/mambaforge:latest
LABEL io.github.snakemake.containerized="true"
LABEL io.github.snakemake.conda_env_hash="d9c34f75f4249d9c5bf03295059c6f6008db227e1418912c0c09cbbfbaa3007c"

# install app dependencies
RUN apt-get update && apt install build-essential -y

# Step 1: Retrieve conda environments

# Conda environment:
#   source: workflow/envs/antismash.yaml
#   prefix: /conda-envs/f353ad0df200ba9e83b18857a7144468
#   name: antismash_env
#   channels:
#     - conda-forge
#     - bioconda
#     - defaults
#   dependencies:
#     - antismash=6.1.1
#     - jinja2
RUN mkdir -p /conda-envs/f353ad0df200ba9e83b18857a7144468
COPY workflow/envs/antismash.yaml /conda-envs/f353ad0df200ba9e83b18857a7144468/environment.yaml

# Conda environment:
#   source: workflow/envs/arts.yaml
#   prefix: /conda-envs/2f035517b36190380bc226958f272eab
#   name: arts_env
#   channels:
#     - conda-forge
#     - defaults
#     - bioconda
#   dependencies:
#     - python==3.8
#     - pip
#     - matplotlib
#     - antismash==6.1.1
#     - jinja2
#     - mafft
#     - pip:
#       - aniso8601==9.0.1
#       - backports-lzma==0.0.14
#       - biopython==1.79
#       - blinker==1.4
#       - click==8.0.1
#       - cssselect==1.1.0
#       - cycler==0.10.0
#       - et-xmlfile==1.1.0
#       - ete3==3.1.2
#       - flask==2.0.1
#       - flask-mail==0.9.1
#       - flask-restful==0.3.9
#       - flask-sockets==0.2.1
#       - gevent==21.8.0
#       - gevent-websocket==0.10.1
#       - greenlet==1.1.1
#       - helperlibs==0.2.1
#       - itsdangerous==2.0.1
#       - jdcal==1.4.1
#       - lxml==4.6.3
#       - markupsafe==2.0.1
#       - networkx==2.6.2
#       - numpy==1.20.2
#       - openpyxl==3.0.7
#       - pillow==8.3.1
#       - pyexcelerator==0.6.4.1
#       - pyparsing==2.4.7
#       - pypng==0.0.20
#       - pyquery==1.4.3
#       - pysvg-py3==0.2.2.post3
#       - python-dateutil==2.8.2
#       - pytz==2021.1
#       - redis==3.5.3
#       - six==1.16.0
#       - straight-plugin==1.5.0
#       - werkzeug==2.0.1
RUN mkdir -p /conda-envs/2f035517b36190380bc226958f272eab
COPY workflow/envs/arts.yaml /conda-envs/2f035517b36190380bc226958f272eab/environment.yaml

# Conda environment:
#   source: workflow/envs/automlst_wrapper.yaml
#   prefix: /conda-envs/3880506184fe2a40becaec4274401f91
#   name: automlst_wrapper
#   channels:
#     - bioconda
#     - conda-forge
#     - anaconda
#   dependencies:
#     - python=2.7
#     - Flask=0.12.2
#     - Flask-Mail=0.9.1
#     - Jinja2=2.10
#     - MarkupSafe=1.0
#     - SQLAlchemy=1.2.5
#     - Werkzeug=0.14.1
#     - biopython=1.70
#     - blinker=1.4
#     - click=6.7
#     - ete3=3.1.1
#     - itsdangerous=0.24
#     - numpy=1.14.2
#     - redis-py=2.10.6
#     - scipy=1.0.0
#     - six=1.11.0
#     - uWSGI=2.0.17
#     - wheel=0.30.0
#     - wsgiref=0.1.2
#     - mafft=7.407
#     - trimal=1.4.1
#     - hmmer=3.1b2
#     - iqtree=1.6.12
#     - prodigal=2.6.3
#     - mash=2.2.1
#     - raxml=8.2.12
#     - wget
#     - unzip
#     - biopython
RUN mkdir -p /conda-envs/3880506184fe2a40becaec4274401f91
COPY workflow/envs/automlst_wrapper.yaml /conda-envs/3880506184fe2a40becaec4274401f91/environment.yaml

# Conda environment:
#   source: workflow/envs/bgc_analytics.yaml
#   prefix: /conda-envs/065809e432cb3905b2e93fe4c29a19c2
#   name: bgc_analytics
#   channels:
#     - conda-forge
#     - bioconda
#     - defaults
#   dependencies:
#     - ncbi-genome-download=0.3.1
#     - jupyterlab>3
#     - ipywidgets>=7.6
#     - jupyter-dash
#     - seaborn
#     - pandas
#     - plotly
#     - python-kaleido
#     - bokeh
#     - openpyxl
#     - nbconvert
#     - ipykernel
#     - ete3
#     - pyarrow
#     - peppy
#     - xlrd >= 1.0.0
#     - zip
#     - altair
#     - itables
#     - scikit-learn
#     - pip
#     - pip:
#         - biopython
#         - jupyterlab-dash
#         - networkx
#         - alive_progress
#         - python-louvain
RUN mkdir -p /conda-envs/065809e432cb3905b2e93fe4c29a19c2
COPY workflow/envs/bgc_analytics.yaml /conda-envs/065809e432cb3905b2e93fe4c29a19c2/environment.yaml

# Conda environment:
#   source: workflow/envs/bigscape.yaml
#   prefix: /conda-envs/6f6f3ea39414ee0158a347375aae97f2
#   name: bigscape
#   channels:
#     - conda-forge
#     - bioconda
#   dependencies:
#     - python=3.7
#     - hmmer
#     - biopython
#     - fasttree
#     - numpy
#     - scipy
#     - networkx
#     - scikit-learn=0.19.2
#     - decorator
#     - wget
#     - unzip
RUN mkdir -p /conda-envs/6f6f3ea39414ee0158a347375aae97f2
COPY workflow/envs/bigscape.yaml /conda-envs/6f6f3ea39414ee0158a347375aae97f2/environment.yaml

# Conda environment:
#   source: workflow/envs/bigslice.yaml
#   prefix: /conda-envs/f75de02f9ad013d32c409fba5c706da3
#   name: bigslice
#   channels:
#     - bioconda
#     - defaults
#   dependencies:
#     - _libgcc_mutex=0.1=main
#     - _openmp_mutex=5.1=1_gnu
#     - bzip2=1.0.8=h7b6447c_0
#     - ca-certificates=2022.10.11=h06a4308_0
#     - certifi=2022.9.24=py310h06a4308_0
#     - hmmer=3.3.2=h87f3376_2
#     - ld_impl_linux-64=2.38=h1181459_1
#     - libffi=3.4.2=h6a678d5_6
#     - libgcc-ng=11.2.0=h1234567_1
#     - libgomp=11.2.0=h1234567_1
#     - libstdcxx-ng=11.2.0=h1234567_1
#     - libuuid=1.41.5=h5eee18b_0
#     - ncurses=6.3=h5eee18b_3
#     - openssl=1.1.1s=h7f8727e_0
#     - pip=22.3.1=py310h06a4308_0
#     - python=3.10.8=h7a1cb2a_1
#     - readline=8.2=h5eee18b_0
#     - setuptools=65.5.0=py310h06a4308_0
#     - sqlite=3.40.0=h5082296_0
#     - tk=8.6.12=h1ccaba5_0
#     - tzdata=2022g=h04d1e81_0
#     - wheel=0.37.1=pyhd3eb1b0_0
#     - xz=5.2.8=h5eee18b_0
#     - zlib=1.2.13=h5eee18b_0
#     - wget
#     - unzip
#     - pip:
#       - biopython==1.76
#       - joblib==1.2.0
#       - numpy==1.23.5
#       - pandas==1.5.2
#       - pysqlite3==0.5.0
#       - python-dateutil==2.8.2
#       - pytz==2022.6
#       - scikit-learn==1.2.0
#       - scipy==1.9.3
#       - six==1.16.0
#       - threadpoolctl==3.1.0
#       - git+https://github.com/medema-group/bigslice.git@v1.1.1
RUN mkdir -p /conda-envs/f75de02f9ad013d32c409fba5c706da3
COPY workflow/envs/bigslice.yaml /conda-envs/f75de02f9ad013d32c409fba5c706da3/environment.yaml

# Conda environment:
#   source: workflow/envs/cblaster.yaml
#   prefix: /conda-envs/28be8691e392102d9d4ce62c2aff767c
#   name: cblaster
#   channels:
#     - conda-forge
#     - default
#     - bioconda
#   dependencies:
#     - diamond==2.0.15
#     - pip
#     - pip:
#       - cblaster==1.3.12
RUN mkdir -p /conda-envs/28be8691e392102d9d4ce62c2aff767c
COPY workflow/envs/cblaster.yaml /conda-envs/28be8691e392102d9d4ce62c2aff767c/environment.yaml

# Conda environment:
#   source: workflow/envs/checkm.yaml
#   prefix: /conda-envs/4671b1666196e23c958b1c6bd7e38440
#   name: checkm
#   channels:
#     - conda-forge
#     - bioconda
#     - defaults
#   dependencies:
#     - checkm-genome=1.1.3
#     - wget
#     - tar
RUN mkdir -p /conda-envs/4671b1666196e23c958b1c6bd7e38440
COPY workflow/envs/checkm.yaml /conda-envs/4671b1666196e23c958b1c6bd7e38440/environment.yaml

# Conda environment:
#   source: workflow/envs/deeptfactor.yaml
#   prefix: /conda-envs/49345ab22685dc2161f5a5e9fb34ef8e
#   name: deeptfactor
#   channels:
#     - pytorch
#     - gurobi
#     - defaults
#     - conda-forge
#     - bioconda
#   dependencies:
#     - _libgcc_mutex=0.1=main
#     - _pytorch_select=0.1=cpu_0
#     - backcall=0.2.0=py_0
#     - blas=1.0=mkl
#     - ca-certificates=2020.7.22=0
#     - certifi=2020.6.20=py36_0
#     - cffi=1.14.2=py36he30daa8_0
#     - cudatoolkit=10.0.130=0
#     - cycler=0.10.0=py36_0
#     - dbus=1.13.16=hb2f20db_0
#     - decorator=4.4.2=py_0
#     - expat=2.2.9=he6710b0_2
#     - fontconfig=2.13.0=h9420a91_0
#     - freetype=2.10.2=h5ab3b9f_0
#     - glib=2.65.0=h3eb4bd4_0
#     - gst-plugins-base=1.14.0=hbbd80ab_1
#     - gstreamer=1.14.0=hb31296c_0
#     - icu=58.2=he6710b0_3
#     - intel-openmp=2020.2=254
#     - ipykernel=5.3.4=py36h5ca1d4c_0
#     - ipython=7.16.1=py36h5ca1d4c_0
#     - ipython_genutils=0.2.0=py36_0
#     - jedi=0.17.2=py36_0
#     - joblib=0.16.0=py_0
#     - jpeg=9b=h024ee3a_2
#     - jupyter_client=6.1.6=py_0
#     - jupyter_core=4.6.3=py36_0
#     - kiwisolver=1.2.0=py36hfd86e86_0
#     - ld_impl_linux-64=2.33.1=h53a641e_7
#     - libedit=3.1.20191231=h14c3975_1
#     - libffi=3.3=he6710b0_2
#     - libgcc-ng=9.1.0=hdf63c60_0
#     - libgfortran-ng=7.3.0=hdf63c60_0
#     - libpng=1.6.37=hbc83047_0
#     - libsodium=1.0.18=h7b6447c_0
#     - libstdcxx-ng=9.1.0=hdf63c60_0
#     - libuuid=1.0.3=h1bed415_2
#     - libxcb=1.14=h7b6447c_0
#     - libxml2=2.9.10=he19cac6_1
#     - matplotlib=3.1.1=py36h5429711_0
#     - mkl=2020.2=256
#     - mkl-service=2.3.0=py36he904b0f_0
#     - mkl_fft=1.1.0=py36h23d657b_0
#     - mkl_random=1.1.1=py36h0573a6f_0
#     - ncurses=6.2=he6710b0_1
#     - ninja=1.10.1=py36hfd86e86_0
#     - numpy=1.19.1=py36hbc911f0_0
#     - numpy-base=1.19.1=py36hfa32c7d_0
#     - openssl=1.1.1g=h7b6447c_0
#     - parso=0.7.0=py_0
#     - pcre=8.44=he6710b0_0
#     - pexpect=4.8.0=py36_0
#     - pickleshare=0.7.5=py36_0
#     - pip=20.2.2=py36_0
#     - prompt-toolkit=3.0.7=py_0
#     - ptyprocess=0.6.0=py36_0
#     - pycparser=2.20=py_2
#     - pygments=2.6.1=py_0
#     - pyparsing=2.4.7=py_0
#     - pyqt=5.9.2=py36h05f1152_2
#     - python=3.6.12=hcff3b4d_2
#     - python-dateutil=2.8.1=py_0
#     - pytz=2020.1=py_0
#     - pyzmq=19.0.1=py36he6710b0_1
#     - qt=5.9.7=h5867ecd_1
#     - readline=8.0=h7b6447c_0
#     - scikit-learn=0.21.3=py36hd81dba3_0
#     - scipy=1.5.2=py36h0b6359f_0
#     - setuptools=49.6.0=py36_0
#     - sip=4.19.8=py36hf484d3e_0
#     - six=1.15.0=py_0
#     - sqlite=3.33.0=h62c20be_0
#     - tk=8.6.10=hbc83047_0
#     - tornado=6.0.4=py36h7b6447c_1
#     - traitlets=4.3.3=py36_0
#     - wcwidth=0.2.5=py_0
#     - wheel=0.35.1=py_0
#     - xz=5.2.5=h7b6447c_0
#     - zeromq=4.3.2=he6710b0_3
#     - zlib=1.2.11=h7b6447c_3
#     - pytorch=1.2.0=py3.6_cuda10.0.130_cudnn7.6.2_0
#     - pip:
#       - biopython==1.78
#       - torch==1.2.0
RUN mkdir -p /conda-envs/49345ab22685dc2161f5a5e9fb34ef8e
COPY workflow/envs/deeptfactor.yaml /conda-envs/49345ab22685dc2161f5a5e9fb34ef8e/environment.yaml

# Conda environment:
#   source: workflow/envs/eggnog.yaml
#   prefix: /conda-envs/107fa35f7f5cdef7817b23e821a46f6f
#   name: eggnog
#   channels:
#     - conda-forge
#     - bioconda
#   dependencies:
#     - python>=3.7
#     - biopython=1.76
#     - psutil=5.7.0
#     - xlsxwriter
#     - sqlite>=3.8.2
#     - eggnog-mapper=2.1.6
RUN mkdir -p /conda-envs/107fa35f7f5cdef7817b23e821a46f6f
COPY workflow/envs/eggnog.yaml /conda-envs/107fa35f7f5cdef7817b23e821a46f6f/environment.yaml

# Conda environment:
#   source: workflow/envs/fastani.yaml
#   prefix: /conda-envs/9e62724cd5cfd7be1542e4d81119a48c
#   name: fastani
#   channels:
#     - conda-forge
#     - bioconda
#     - defaults
#   dependencies:
#     - fastani=1.33
#     - gsl=2.7=he838d99_0
RUN mkdir -p /conda-envs/9e62724cd5cfd7be1542e4d81119a48c
COPY workflow/envs/fastani.yaml /conda-envs/9e62724cd5cfd7be1542e4d81119a48c/environment.yaml

# Conda environment:
#   source: workflow/envs/gtdbtk.yaml
#   prefix: /conda-envs/8f484babe423b120651944f71dbc0dd2
#   name: gtdbtk
#   channels:
#     - conda-forge
#     - bioconda
#     - defaults
#   dependencies:
#     - gtdbtk=2.1.0
#     - wget
#     - tar
RUN mkdir -p /conda-envs/8f484babe423b120651944f71dbc0dd2
COPY workflow/envs/gtdbtk.yaml /conda-envs/8f484babe423b120651944f71dbc0dd2/environment.yaml

# Conda environment:
#   source: workflow/envs/mash.yaml
#   prefix: /conda-envs/cfb3246b4c051b971e47d79628582a19
#   name: mash
#   channels:
#     - conda-forge
#     - bioconda
#     - defaults
#   dependencies:
#     - mash=2.3
RUN mkdir -p /conda-envs/cfb3246b4c051b971e47d79628582a19
COPY workflow/envs/mash.yaml /conda-envs/cfb3246b4c051b971e47d79628582a19/environment.yaml

# Conda environment:
#   source: workflow/envs/prokka.yaml
#   prefix: /conda-envs/62a64d6cebaaeccd8f5371f02ca8e9c8
#   name: prokka_env
#   channels:
#     - conda-forge
#     - bioconda
#     - defaults
#   dependencies:
#     - prokka=1.14.6
#     - perl-bioperl=1.7.2
#     - pandas
#     - biopython
#     - wget
RUN mkdir -p /conda-envs/62a64d6cebaaeccd8f5371f02ca8e9c8
COPY workflow/envs/prokka.yaml /conda-envs/62a64d6cebaaeccd8f5371f02ca8e9c8/environment.yaml

# Conda environment:
#   source: workflow/envs/roary.yaml
#   prefix: /conda-envs/4627313edbb641d43e1f032932d0bff4
#   name: roary_env
#   channels:
#     - bioconda
#     - conda-forge
#     - defaults
#     - r
#   dependencies:
#     - roary=3.13.0
#     - r-base
#     - r-ggplot2=3.3.3
RUN mkdir -p /conda-envs/4627313edbb641d43e1f032932d0bff4
COPY workflow/envs/roary.yaml /conda-envs/4627313edbb641d43e1f032932d0bff4/environment.yaml

# Conda environment:
#   source: workflow/envs/seqfu.yaml
#   prefix: /conda-envs/5d12220de6e6e9699ca6c76c54401981
#   name: seqfu
#   channels:
#     - conda-forge
#     - bioconda
#     - defaults
#   dependencies:
#     - seqfu=1.15.3
RUN mkdir -p /conda-envs/5d12220de6e6e9699ca6c76c54401981
COPY workflow/envs/seqfu.yaml /conda-envs/5d12220de6e6e9699ca6c76c54401981/environment.yaml

# Conda environment:
#   source: workflow/envs/dbt-duckdb.yaml
#   prefix: /conda-envs/0ed8477ceb477f29491cb36459d0db3a
#   name: dbt-duckdb
#   channels:
#     - conda-forge
#     - defaults
#   dependencies:
#     - python<=3.9
#     - duckcli>=0.2.1
#     - dbt-core>=1.3.1
#     - sqlfluff~=1.4.3
#     - pip
#     - pip:
#       - git+https://github.com/jwills/dbt-duckdb.git@v1.3.1
RUN mkdir -p /conda-envs/0ed8477ceb477f29491cb36459d0db3a
COPY workflow/envs/dbt-duckdb.yaml /conda-envs/0ed8477ceb477f29491cb36459d0db3a/environment.yaml

# Step 2: Generate conda environments

RUN mamba env create --prefix /conda-envs/f353ad0df200ba9e83b18857a7144468 --file /conda-envs/f353ad0df200ba9e83b18857a7144468/environment.yaml && \
    mamba env create --prefix /conda-envs/2f035517b36190380bc226958f272eab --file /conda-envs/2f035517b36190380bc226958f272eab/environment.yaml && \
    mamba env create --prefix /conda-envs/3880506184fe2a40becaec4274401f91 --file /conda-envs/3880506184fe2a40becaec4274401f91/environment.yaml && \
    mamba env create --prefix /conda-envs/065809e432cb3905b2e93fe4c29a19c2 --file /conda-envs/065809e432cb3905b2e93fe4c29a19c2/environment.yaml && \
    mamba env create --prefix /conda-envs/6f6f3ea39414ee0158a347375aae97f2 --file /conda-envs/6f6f3ea39414ee0158a347375aae97f2/environment.yaml && \
    mamba env create --prefix /conda-envs/f75de02f9ad013d32c409fba5c706da3 --file /conda-envs/f75de02f9ad013d32c409fba5c706da3/environment.yaml && \
    mamba env create --prefix /conda-envs/28be8691e392102d9d4ce62c2aff767c --file /conda-envs/28be8691e392102d9d4ce62c2aff767c/environment.yaml && \
    mamba env create --prefix /conda-envs/4671b1666196e23c958b1c6bd7e38440 --file /conda-envs/4671b1666196e23c958b1c6bd7e38440/environment.yaml && \
    mamba env create --prefix /conda-envs/49345ab22685dc2161f5a5e9fb34ef8e --file /conda-envs/49345ab22685dc2161f5a5e9fb34ef8e/environment.yaml && \
    mamba env create --prefix /conda-envs/107fa35f7f5cdef7817b23e821a46f6f --file /conda-envs/107fa35f7f5cdef7817b23e821a46f6f/environment.yaml && \
    mamba env create --prefix /conda-envs/9e62724cd5cfd7be1542e4d81119a48c --file /conda-envs/9e62724cd5cfd7be1542e4d81119a48c/environment.yaml && \
    mamba env create --prefix /conda-envs/8f484babe423b120651944f71dbc0dd2 --file /conda-envs/8f484babe423b120651944f71dbc0dd2/environment.yaml && \
    mamba env create --prefix /conda-envs/cfb3246b4c051b971e47d79628582a19 --file /conda-envs/cfb3246b4c051b971e47d79628582a19/environment.yaml && \
    mamba env create --prefix /conda-envs/62a64d6cebaaeccd8f5371f02ca8e9c8 --file /conda-envs/62a64d6cebaaeccd8f5371f02ca8e9c8/environment.yaml && \
    mamba env create --prefix /conda-envs/4627313edbb641d43e1f032932d0bff4 --file /conda-envs/4627313edbb641d43e1f032932d0bff4/environment.yaml && \
    mamba env create --prefix /conda-envs/5d12220de6e6e9699ca6c76c54401981 --file /conda-envs/5d12220de6e6e9699ca6c76c54401981/environment.yaml && \
    mamba clean --all -y
