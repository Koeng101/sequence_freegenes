# Init
FROM ubuntu:16.04
MAINTAINER Keoni Gandall "koeng101@gmail.com"

# Update
RUN apt-get update && apt-get install -y apt-transport-https python-pip python-dev curl

# Install requirements
COPY . /app
WORKDIR /app
RUN pip install -r requirements.txt

# Install minimap2
RUN mkdir storage
RUN cd ./storage
RUN curl -L https://github.com/lh3/minimap2/releases/download/v2.17/minimap2-2.17_x64-linux.tar.bz2 | tar -jxvf -
RUN cd ..
 
# Start app
ENTRYPOINT ["python"]
CMD ["wsgi.py"]
