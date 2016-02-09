import os,sys
import fnmatch
import logging
import octopod_config as ocfg

class FileManager:

    def __init__(self):
        self.logger = logging.getLogger(__name__)
        self.logger.info('Creating FileManager.')
        paths = ocfg.data_paths
        self.keys = paths.keys()
        self.files = {}
        for key in self.keys:
            self.files[key] = []
            for root,dirnames,filenames in os.walk(paths[key]) :
                for filename in fnmatch.filter(filenames,'*.%s'%(ocfg.raw_data_extension)):
                    full_fn = os.path.join(root,filename)
                    self.files[key].append(full_fn)
                    self.logger.info('Adding %s to set %s.'%(full_fn,key))

    def get(self,key,terms=[]):
        try:
            files = self.files[key]
            if len(terms):
                out = []
                for f in files:
                    for term in terms:
                        if f.lower().find(term.lower())>-1:
                            out.append(f)
            else:
                out = files
            return out
                        
        except KeyError as e:
            self.logger.info('KeyError: "%s". Returning [].'%(e.message))
            return []

    def get_one(self,key,terms=[]):
        files = self.get(key)
        if len(files):
            matches = 0
            winner = 0
            for idx,f in enumerate(files):
                temp_matches = 0
                for term in terms:
                    if f.lower().find(term.lower())>-1:
                        temp_matches = temp_matches + 1
                if temp_matches>matches:
                    matches = temp_matches
                    winner = idx
            return files[winner]
