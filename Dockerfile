# Init
FROM fnndsc/ubuntu-python3
MAINTAINER Keoni Gandall "koeng101@gmail.com"


COPY . /app
WORKDIR /app

RUN pip install -r requirements.txt


# Install minimap2
#RUN curl -L https://github.com/lh3/minimap2/releases/download/v2.17/minimap2-2.17_x64-linux.tar.bz2 | tar -jxvf -
RUN ls

# Start app
CMD gunicorn wsgi:app --workers=4
