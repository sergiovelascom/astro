from numpy import *
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import matplotlib.cm as cm

import  pymc

data=genfromtxt('result.txt', names=True, delimiter='\t', dtype=None)


star="a_And"

starflux= where(data['object'] == star)
datastar=data[starflux]

ifilter= where(datastar['filt'] =='I')
istar=datastar[ifilter]
vfilter= where(data['filt'][starflux] =='V')
vstar=datastar[vfilter]

gain2_v= where(vstar['gain'] == 2)
gain50_v= where(vstar['gain'] == 50)
gain100_v= where(vstar['gain'] == 100)
gain200_v= where(vstar['gain'] == 200)
gain300_v= where(vstar['gain'] == 300)

gain2_i= where(istar['gain'] == 2)
gain50_i= where(istar['gain'] == 50)
gain100_i= where(istar['gain'] == 100)
gain200_i= where(istar['gain'] == 200)
gain300_i= where(istar['gain'] == 300)

gain2_i=istar[gain2_i]
gain50_i=istar[gain50_i]
gain100_i=istar[gain100_i]
gain200_i=istar[gain200_i]
gain300_i=istar[gain300_i]

gain2_v=vstar[gain2_v]
gain50_v=vstar[gain50_v]
gain100_v=vstar[gain100_v]
gain200_v=vstar[gain200_v]
gain300_v=vstar[gain300_v]

ganv=[gain2_v,gain50_v,gain100_v,gain200_v,gain300_v]
gani=[gain2_i,gain50_i,gain100_i,gain200_i,gain300_i]

cmap = cm.Blues	
cmap2= cm.Reds

color=["blue","red","yellow","black","orange"]
a=[2,50,100,200,300]

def func(x, a, b, c,d,e):
    return a*x**4 + b*x**3+ c*x**2 + d*x +e

for i in arange(0,5,1):
  vflux=ganv[i]
  iflux=gani[i]
  if iflux.size > 0:
    #plt.scatter(data['time'][vflux],(data['maximo'][vflux]),color=cmap(i / float(5)))
    #plt.scatter(data['time'][iflux],(data['maximo'][iflux]),color=cmap2(i / float(5)))
    popt_i, pcov_i = curve_fit(func, iflux['time'],iflux['maximo'])
    label="I band Gain "+str(a[i])
    plt.plot(sort(iflux['time']), sort( polyval(popt_i, iflux['time'])),color=cmap2(i +2/ float(5)),label=label)
  if vflux.size > 0:
    popt_v, pcov_v = curve_fit(func, vflux['time'],vflux['maximo'])
    label="V band Gain "+str(a[i])
    plt.plot(sort(vflux['time']), sort( polyval(popt_v, vflux['time'])),color=cmap(i +2/ float(5)),label=label)

plt.xlabel('time (ms)')
plt.ylabel('peak max(counts)')
plt.title('a_And V=2.19 I=2.06')


plt.legend(loc="upper right")
plt.savefig('a_And.png')
plt.show()

#def funcl(x, a, b):
    #return  a*x + b 

  

#print 'The chi-square result: ',  popt


    
#priors
#sig = pymc.Uniform('sig', 0.0, 100.0, value=1.)

#a = pymc.Uniform('a', -10.0, 10.0, value= 0.0)
#b = pymc.Uniform('b', -10.0, 10.0, value= 0.0)
#c = pymc.Uniform('c', -10.0, 10.0, value= 0.0)


#model
#@pymc.deterministic(plot=False)
#def mod_quadratic(x=data['time'][bs107], a=a, b=b, c=c):
      #return a*x**2 + b*x + c

#likelihood
#y = pymc.Normal('y', mu=mod_quadratic, tau=1.0/sig**2, value=data['data'][bs107], observed=True)




#plt.plot(sort(data['time'][bs107]), sort( polyval(popt, data['time'][bs107])),  color='red',label='Gain 2')


