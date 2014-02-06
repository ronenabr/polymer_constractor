import numpy
from numpy.linalg import det, eig
from numpy import array
from numpy import * 

vs = [ [1,0,0],
        [0,1,0],
        [0,0,1],
        [1,1,0],
        [1,0,1],
        [0,1,1],
        [-1,1,0],
        [-1,0,1],
        [0,-1,1],
        [1,1,1],
        [1,1,-1],
        [1,-1,1],
        [-1,1,1]]

vs2 = [array(v, dtype=int) for v in vs]
vs2.extend([(-1)*v for v in vs2])



g = [0]*9


        
lt = []
def find_g(g,idx):
    if idx==9:
        a = array(g, dtype=int).reshape([3,3])        
        d = det(a)
        if d==1 or d==-1:
            lt.append(a)
        return
    g[idx] = 0
    find_g(g,idx+1)
    g[idx] = 1
    find_g(g,idx+1)
    g[idx] = -1
    find_g(g,idx+1)
find_g(g,0)
print len(lt)



l2 = []
i = 0
j = 0
for a in lt:
    i+=1
    flag = True
    for u in vs2:
        sub_flag = True  in [array_equal(a.dot(u),v)  for v in vs2]
        if not sub_flag :
            j+=1
            flag = False
            break  
    if flag:
        print i, j
        l2.append(a)
print "Len,", len(l2)

numpy.save("matrices.npy",l2)




