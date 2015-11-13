FROM iplantc/dnasub_base
WORKDIR /home
MAINTAINER Matthew Vaughn <vaughn@tacc.utexas.edu>

# Install Samtools 0.1.19
RUN cd /tmp && wget --no-check-certificate "http://sourceforge.net/projects/samtools/files/samtools/0.1.19/samtools-0.1.19.tar.bz2" && tar -jxvf samtools-0.1.19.tar.bz2 && cd samtools-0.1.19 && make && cp samtools /opt/bin/ && cp bcftools/bcftools /opt/bin/ && cp misc/* /opt/bin/ && cd /tmp && rm -rf /tmp/*

# Install FASTX_toolkit
RUN cd /tmp && wget --no-check-certificate "http://hannonlab.cshl.edu/fastx_toolkit/fastx_toolkit_0.0.13_binaries_Linux_2.6_amd64.tar.bz2" && tar -jxvf fastx_toolkit_0.0.13_binaries_Linux_2.6_amd64.tar.bz2 && mv bin/* /opt/bin/ && cd /tmp && rm -rf /tmp/*

# Install Picard 1.107
# Set TACC_PICARD_DIR for compatibility
RUN cd /tmp && wget --no-check-certificate -O picard-tools-1.120.zip "https://github.com/broadinstitute/picard/releases/download/1.120/picard-tools-1.120.zip" && unzip picard-tools-1.120.zip && cp -R picard-tools-1.120 /opt && echo "export TACC_PICARD_DIR=/opt/picard-tools-1.120" >> ~/.bashrc && cd /tmp && rm -rf /tmp/*

# Install TopHat 2.0.11.1
# https://ccb.jhu.edu/software/tophat/downloads/tophat-2.0.11.Linux_x86_64.tar.gz
RUN cd /tmp && wget --no-check-certificate "https://ccb.jhu.edu/software/tophat/downloads/tophat-2.0.11.Linux_x86_64.tar.gz" && tar -zxvf tophat-2.0.11.Linux_x86_64.tar.gz && mv tophat-2.0.11.Linux_x86_64/* /opt/bin/ && cd /tmp && rm -rf /tmp/*

# Install Bowtie 2.1.0.0
RUN cd /tmp && wget --no-check-certificate "http://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.1.0/bowtie2-2.1.0-linux-x86_64.zip" && unzip bowtie2-2.1.0-linux-x86_64.zip && mv bowtie2-2.1.0/* /opt/bin/ && cd /tmp && rm -rf /tmp/*

# Install Cufflinks 2.1.1
RUN cd /tmp && wget --no-check-certificate "http://cole-trapnell-lab.github.io/cufflinks/assets/downloads/cufflinks-2.1.0.Linux_x86_64.tar.gz" && tar -zxvf cufflinks-2.1.0.Linux_x86_64.tar.gz && mv cufflinks-2.1.0.Linux_x86_64/* /opt/bin/ && cd /tmp && rm -rf /tmp/*
