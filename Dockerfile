# Init
FROM centos:centos7
MAINTAINER Keoni Gandall "koeng101@gmail.com"

# Update
#RUN yum install -y build-essential 
RUN yum install -y libssl-dev
RUN yum install -y python3-dev
RUN yum install -y python3-pip
RUN yum install -y curl
RUN yum install -y samtools

# Install requirements
COPY . /app
WORKDIR /app
RUN pip3 install -r requirements.txt

# Install minimap2
RUN curl -L https://github.com/lh3/minimap2/releases/download/v2.17/minimap2-2.17_x64-linux.tar.bz2 | tar -jxvf -
RUN ls

# Start app
CMD gunicorn wsgi:app --workers=4
