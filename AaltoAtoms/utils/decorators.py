from datetime import datetime

def timed_log(log_msg):
    def time_added(*args, **kwargs):
        return f'[{datetime.now()}]:\n\t {log_msg(*args, **kwargs)}'
    return time_added

from functools import wraps
from time import time


class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKCYAN = '\033[96m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

def timing(f):
    @wraps(f)
    def wrap(*args, **kw):
        ts = time()
        result = f(*args, **kw)
        te = time()

        print('func: %r\n\targs: [%r, %r] \n\t\ttook: %2.4f sec' % \
          (f.__name__, args, kw, te-ts))
        return result
    return wrap
