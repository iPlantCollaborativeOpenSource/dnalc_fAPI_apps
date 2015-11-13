
FROM centos:centos5
WORKDIR /home
MAINTAINER Matthew Vaughn <vaughn@tacc.utexas.edu>

ENV _APP DNALC Lonestar4 runtime base
ENV _AUTHOR Matthew Vaughn
ENV _LICENSE MIT
ENV _VERSION 0.0.1

# Set up /opt and permantly set PATH
RUN mkdir -p /opt/bin && echo "PATH=\$PATH:/opt/bin" >> /etc/profile

# Basics
RUN yum -y install wget which bzip2 java-1.7.0-openjdk.x86_64

# Install git as per https://gist.github.com/eddarmitage/2001099
RUN wget --no-check-certificate http://packages.sw.be/rpmforge-release/rpmforge-release-0.5.2-2.el5.rf.i386.rpm
RUN rpm --import http://apt.sw.be/RPM-GPG-KEY.dag.txt && rpm -K rpmforge-release-0.5.2-2.el5.rf.i386.rpm && rpm -i rpmforge-release-0.5.2-2.el5.rf.i386.rpm && rm rpmforge-release-0.5.2-2.el5.rf.i386.rpm
RUN yum -y install git-gui

# Python 2.7.1 into /opt as per https://gist.github.com/stantonk/4082043
# Must be done this way to avoid breaking Centos5's yum, which
# depends on Python 2.4.x
RUN yum -y groupinstall "Development tools"
RUN yum -y install zlib-devel bzip2-devel openssl-devel ncurses-devel
RUN wget --no-check-certificate http://www.python.org/ftp/python/2.7.1/Python-2.7.1.tar.bz2 && tar xf Python-2.7.1.tar.bz2
RUN cd Python-2.7.1 && ./configure --prefix=/opt && make && make altinstall /opt/bin/python2.7
RUN cd /tmp && wget --no-check-certificate http://pypi.python.org/packages/source/d/distribute/distribute-0.6.27.tar.gz && tar xf distribute-0.6.27.tar.gz && cd distribute-0.6.27 && /opt/bin/python2.7 setup.py install && /opt/bin/easy_install-2.7 virtualenv

# Install R 3.2.1
RUN yum -y install readline-devel sqlite
RUN cd /tmp && wget --no-check-certificate "https://cran.r-project.org/src/base/R-3/R-3.2.1.tar.gz" && tar -zxvf R-3.2.1.tar.gz && cd R-3.2.1 && ./configure --prefix=/opt --with-x=no && make && make install && cd /tmp && rm -rf /tmp/*
# Install BioC and cummerbund
# These are not identical to the ones on Lonestar or in the bundle
RUN echo -e "source(\"http://bioconductor.org/biocLite.R\");\nbiocLite();\nbiocLite(c(\"cummeRbund\"));" > bio.R
RUN /opt/bin/Rscript bio.R

