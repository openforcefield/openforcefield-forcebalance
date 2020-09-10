import sys
from pymongo import MongoClient

client = MongoClient('mongodb://localhost')
client.drop_database(sys.argv[1])
