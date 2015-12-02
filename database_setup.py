# sqlite> .schema volumes
# CREATE TABLE volumes(vid text, sid text, eid text, date text, eye text, info text, processed int);
# sqlite> .schema experiments
# CREATE TABLE experiments(eid text, sid text, date text);
# sqlite> .schema subjects
# CREATE TABLE subjects(sid text);
# sqlite> 

import sys,os
import octopod_config as ocfg
import sqlite3
import logging
logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)

db_location = ocfg.db_path
logger.info('Setting up database in %s.'%db_location)

conn = sqlite3.connect(os.path.join(db_location,'octopod.db'))
c = conn.cursor()

c.execute('''CREATE TABLE volumes(vid text, sid text, eid text, date text, eye text, info text)''') 

