import os
import shutil
import datetime

def nowStr():
    return datetime.datetime.now().strftime("%Y%m%d%H%M%S")

def backupDirectory(dirname):
    if not os.path.isdir(dirname):
        sys.exit('OCTTools.misc.backupDirectory: %s is not a directory')

    [head,tail] = os.path.split(dirname)
    newdirname = os.path.join(head,os.path.join('bak','%s_%s'%(tail,nowStr())))
    #newdirname = os.path.join(head,'bak_%s_%s'%(tail,nowStr()))
    print 'Copying %s to %s'%(dirname,newdirname)
    try:
        shutil.copytree(dirname,newdirname)
    except Exception as e:
        # Copytree throws an exception, after performing the copy,
        # due to an inability to preserve file times/dates. That's
        # why we do this silly thing.
        pass
        #print 'OCTTools.misc.backupDirectory error:'
        #print e
        #print '\t%s'%e

def deleteDirectoryContents(dirname):
    if not os.path.isdir(dirname):
        sys.exit('OCTTools.misc.deleteDirectory: %s is not a directory')

    print 'Deleting %s/*'%dirname
    try:
        shutil.rmtree(dirname)
    except Exception as e:
        print 'OCTTools.misc.deleteDirectory error:'
        print '\t%s'%e

    os.makedirs(dirname)
