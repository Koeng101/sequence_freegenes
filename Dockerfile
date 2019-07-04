# Init
FROM python:3.7-alpine
MAINTAINER Keoni Gandall "koeng101@gmail.com"


COPY . /app
WORKDIR /app

# Update
#RUN yum install -y build-essential 
RUN apk add --no-cache --virtual=.build_dependencies musl-dev gcc python3-dev libffi-dev linux-headers && \
    pip install -r requirements.txt && \
    invoke app.dependencies.install && \
    ( \
        if [ "$INCLUDE_POSTGRESQL" = 'true' ]; then \
            apk add --no-cache libpq && \
            apk add --no-cache --virtual=.build_dependencies postgresql-dev && \
            pip install psycopg2 ; \
        fi \
    ) && \
    ( if [ "$INCLUDE_UWSGI" = 'true' ]; then pip install uwsgi ; fi ) && \
    rm -rf ~/.cache/pip && \
    apk del .build_dependencies


# Install minimap2
#RUN curl -L https://github.com/lh3/minimap2/releases/download/v2.17/minimap2-2.17_x64-linux.tar.bz2 | tar -jxvf -
RUN ls

# Start app
CMD gunicorn wsgi:app --workers=4
