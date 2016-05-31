import time

class takeTime:
    def __enter__(self):
        self.t0 = time.time()
    def __exit__(self):
        print str(time.time() - self.t0)

def countALittle(maximum):
    for i in range(maximum):
        print i

# with