# Init
FROM ubuntu:16.04
MAINTAINER Keoni Gandall "koeng101@gmail.com"

# Update
RUN apt-get update 
RUN apt-get install -y build-essential libssl-dev libffi-dev python3-dev apt-transport-https python3-pip curl samtools

# Install requirements
COPY . /app
WORKDIR /app
RUN pip3 install -r requirements.txt

# Install minimap2
#RUN curl -L https://github.com/lh3/minimap2/releases/download/v2.17/minimap2-2.17_x64-linux.tar.bz2 | tar -jxvf -
RUN ls

# Start app
CMD gunicorn wsgi:app --workers=4

