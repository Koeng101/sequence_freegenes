# Init
FROM ubuntu:16.04
MAINTAINER Keoni Gandall "koeng101@gmail.com"

# Update
RUN apt-get update && apt-get install -y apt-transport-https python-pip python-dev

# Install requirements
COPY . /app
WORKDIR /app
RUN pip install -r requirements.txt

# Start app
ENTRYPOINT ["python"]
CMD ["wsgi.py"]
