import pymongo
import logging
import octopod_config as ocfg

class Database:

    def __init__(self):
        self.client = pymongo.MongoClient()
        
