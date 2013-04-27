
import numpy as np

dt = np.dtype([('x', np.float64), ('y', np.float64),
               ('cell', np.int32)])
print dt

arr = np.zeros([4, 4], dtype=dt)
print arr
print arr['x']
print arr[1,1]['x']

x = [1, 2, 3, 4]
y = [10, 20, 30, 40]
cell = [1, 1, 1, 1]
arr[1, :]['x'] = x
print arr

def f(u1, u2, v1, v2):
    print 'u1', u1

f(1, 2, 3, 4)
x = [1, 2, 3, 4]
f(*x)
