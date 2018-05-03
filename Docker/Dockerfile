FROM ubuntu:18.04

LABEL software="CharGer"
LABEL description="Characterization of Germline variants"

WORKDIR /opt
ENV PY_VERSION 2.7.14

RUN  apt-get update \
     \
     # Install python2
     && apt-get install --no-install-recommends -y gcc make libreadline-dev libncursesw5-dev libssl-dev libsqlite3-dev libgdbm-dev libc6-dev libbz2-dev zlib1g-dev \
     && apt-get install --no-install-recommends -y wget git ca-certificates \
     \
     && wget -O - https://www.python.org/ftp/python/${PY_VERSION}/Python-${PY_VERSION}.tgz | tar xzf - \
     && cd Python-${PY_VERSION} \
     && ./configure \
     && make install \
     && cd .. \
     && rm -rf Python-${PY_VERSION} \
     \
     # Install pip
     && wget -O get-pip.py https://bootstrap.pypa.io/get-pip.py \
     && python ./get-pip.py \
     && rm -f ./get-pip.py \
     \
     # Install CharGer
     && git clone https://github.com/ding-lab/charger.git \
     && cd charger \
     && rm -rf .git* \
     && pip install . \
     \
     # Cleanup
     && apt-get remove -y '.*-dev' \
     && apt-get autoremove -y \
     && rm -rf /var/lib/apt/lists/* \

ENV PATH=/opt/charger/bin:${PATH}
