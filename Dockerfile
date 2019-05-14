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

# Get recon database
RUN wget https://zenodo.org/record/583326/files/Recon2M.2_MNX_Entrez_Gene.xml?download=1
RUN mv Recon2M.2_MNX_Entrez_Gene.xml?download=1 recon2m2.xml; mv recon2m2.xml dock/source/

# Get hmdb database
RUN wget http://www.hmdb.ca/system/downloads/current/hmdb_metabolites.zip
RUN unzip hmdb_metabolites.zip; mv hmdb_metabolites.xml dock/source/; rm hmdb_metabolites.zip



# Run metabonet
ENTRYPOINT ["metabonet"]
