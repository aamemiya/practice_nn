import sys
import os 
import numpy as np
import matplotlib.pyplot as plt
  
dtype=np.dtype('>f') ### big endian
fillvalue=-9.99e30

xdim=1
ydim=1
nsmp=10
obs_error_std=np.array([0.0])

#tmp=np.zeros((4,1))
#tmp.fill(fillvalue)
#print(tmp)
#quit()

def function_truth(x_in) :
  y_out = x_in[:,0:1] ** 2 
  return y_out

def out_binary(array,fname):
   nsmp=array.shape[0]
   buf=np.zeros((nsmp,1))
   buf.fill(np.nan)
   arrayout=array
#   arrayout=np.append(buf,array,axis=-1)
#   arrayout=np.append(arrayout,buf,axis=-1)
   arrayout=arrayout.astype(dtype)
   arrayout.tofile(fname)

def in_binary(fname,ndim):
   vread=np.fromfile(fname,dtype=dtype)
#   ndim=np.where(np.isnan(vread))[0][1]-1
#   nsmp=int(len(vread)/(ndim+2))
#   array=vread.reshape((nsmp,ndim+2))[:,1:-1]
   nsmp=int(len(vread)/ndim)
   array=vread.reshape((nsmp,ndim))
   return array

x=-5+10*np.random.random((nsmp,xdim))
y=function_truth(x) + obs_error_std * np.random.randn(nsmp,ydim)

out_binary(x,"x.dat")
out_binary(y,"y.dat")

xread=in_binary("x.dat",1)
yread=in_binary("y.dat",1)

plt.plot(xread,yread,"bo")
fout="sample.png"
print("plot "+ fout)
plt.savefig(fout)
plt.clf
