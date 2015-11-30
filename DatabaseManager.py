import pymongo
import logging
logging.basicConfig(level=logging.INFO)

import octopod_config as ocfg

class Database:

    def __init__(self):
        self.client = pymongo.MongoClient()
        
