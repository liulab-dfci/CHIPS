#!/usr/bin/env python
import os
import sys
from optparse import OptionParser

# def _check_md5(pathDir):
#     """check the completeness of transferring
#     """
#     ID = os.path.basename(pathDir).replace('dataset', '')
#     md5file = os.path.join(pathDir, ID+'.md5')
#     if not os.path.exists(md5file):
#         return None
#     os.chdir(pathDir)
#     check_res = commands.getoutput('md5sum -c %s'%md5file).replace('%s.md5: FAILED'%ID, '')
#     os.chdir(pathDir.replace('final_chilin_result/dataset%s'%ID, ''))
#     if ('FAILED' in check_res) or ('No such file' in check_res):
#         return None
#     else:
#         return True
			
def _again_md5(options):
    """ remove the old md5 and get a new
    """
    ID = options.ID
    pathDir = options.dir
    # ID = os.path.basename(pathDir).replace('dataset', '')
    md5file = os.path.join(pathDir, ID+'.md5')
    # os.system('rm -fr %s' % md5file)
    os.chdir(pathDir)
    os.system('find ./ -type f -print0 | xargs -0 md5sum > %s.md5'%ID)
    # os.chdir(pathDir.replace('final_chilin_result/dataset%s'%ID, 'mainPlace'))

def main():
    USAGE=""
    optparser = OptionParser(usage=USAGE)
    optparser.add_option("-d", "--dir", help="direction")
    optparser.add_option("-I", "--ID", help="ID")
    (options, args) = optparser.parse_args(sys.argv)
    _again_md5(options)



if __name__ == '__main__':
    main()

