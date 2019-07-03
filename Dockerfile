# Init
FROM ubuntu:16.04
MAINTAINER Keoni Gandall "koeng101@gmail.com"

# Update
RUN apt-get update && apt-get install -y apt-transport-https python3-pip python3-dev curl samtools

# Install requirements
COPY . /app
WORKDIR /app
RUN pip3 install -r requirements.txt

# Install minimap2
RUN mkdir storage
RUN cd ./storage
RUN curl -L https://github.com/lh3/minimap2/releases/download/v2.17/minimap2-2.17_x64-linux.tar.bz2 | tar -jxvf -
RUN ls
RUN cd ..
 
# Start app
ENTRYPOINT ["python3"]
CMD ["wsgi.py"]
