# Get a python3 environment
FROM python:3

RUN apt-get update \
 && apt-get install -y --no-install-recommends \
   libfftw3-dev \
   gcc && apt-get clean \
 && rm -rf /var/lib/apt/lists/*

# Get metabonet and install
RUN git clone https://github.com/tcameronwaller/metabonet.git
RUN cd metabonet; python3 setup.py install

# Run metabonet
ENTRYPOINT ["metabonet"]
