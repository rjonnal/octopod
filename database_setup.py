# sqlite> .schema volumes
# CREATE TABLE volumes(vid text, sid text, eid text, date text, eye text, info text, processed int);
# sqlite> .schema experiments
# CREATE TABLE experiments(eid text, sid text, date text);
# sqlite> .schema subjects
# CREATE TABLE subjects(sid text);
# sqlite> 

import octopod_config as ocfg

db_location = ocfg.db_path
