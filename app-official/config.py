import boto3
import os
from botocore.client import ClientError

# Setup description
API_TITLE = os.environ['API_TITLE']
API_DESCRIPTION = os.environ['API_DESCRIPTION']

# Setup Postgres
URL = os.environ['URL']

# Setup Bucket
SPACES = boto3.session.Session().client('s3',
                        region_name=os.environ['REGION_NAME'],
                        endpoint_url=os.environ['ENDPOINT_URL'],
                        aws_access_key_id=os.environ['AWS_ACCESS_KEY_ID'],
                        aws_secret_access_key=os.environ['AWS_SECRET_ACCESS_KEY'])
BUCKET = os.environ['BUCKET']

# Setup Redis
REDIS_HOST = os.environ['REDIS_HOST']
REDIS_PORT = os.environ['REDIS_PORT']
REDIS_PASSWORD = os.environ['REDIS_PASSWORD']

