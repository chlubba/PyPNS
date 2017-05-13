from contextlib import contextmanager
import time

@contextmanager
def takeTime(action):
    takeTime.t0 = time.time()
    print 'Elapsed time to %s...'% action,
    yield
    print '%.2f s' % (time.time() - takeTime.t0)

