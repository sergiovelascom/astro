import pyfits
import scipy.ndimage as snd
from numpy import *
import sys
import re
import matplotlib.pyplot as plt
plt.ion()
pathraw="./RAW_files/"
pathproc="./PROCESS/"
m = sys.argv[1]
filt = sys.argv[2]


data = [line.strip() for line in open("./listam_"+str(m)+"_"+filt, 'r')]
bias=pyfits.getdata("Bias_master_g1m"+str(m)+".fits")
flat=pyfits.getdata("FlatSky_master_Norm_"+str(m)+"_"+filt+".fits")


psf=7#poner 1/2
ran=2.5*psf #Rango busqueda LI
sizeimg = shape(bias)[1]
inner=2.5*psf
outer=5*psf

n=len(data)
xc=()
yc=()


f = open("output_"+str(m)+"_"+filt+".txt",'w')
a= ("FILE", "object",  "i_flux", "i_flux_ex", "v_flux", "v_flux_ex","psf", "inner","outer","gain", "filt", "air", "temp", "nimag", "time", "imag_number","yc_obj", "xc_obj", "xc", "yc","maximo","counts_data","npix_data","data","counts_noise","npix_noise","noise")
print "\t".join(a)

print >> f,"\t".join(a)


i_obj = { 'bsc_107': 5, 'bsc205': 4.05, 'bsc182': 0.54 , 'HD250092':9.5, 'a_And':2.19, 'GSC02500':14.5, 'HR928': 6.8, 'TYC1927':11, 'Orion': "nan" , 'Neptuno': "nan", 'Jupiter': "nan", 'M15': "nan", 'Urano': "nan",'Crab': "nan",'NGC4988':"nan",'NGC4535':"nan",'NGC4147':"nan",'NGC2419':"nan",'NGC2835':"nan"}
v_obj = { 'bsc_107': 6.05, 'bsc205': 4.0, 'bsc182': 2.01 , 'HD250092':10.01, 'a_And':2.06, 'GSC02500':15, 'HR928': 7.0, 'TYC1927':11.97, 'Orion': "nan" , 'Neptuno': "nan", 'Jupiter': "nan", 'M15': "nan", 'Urano': "nan",'Crab': "nan",'NGC4988':"nan",'NGC4535':"nan",'NGC4147':"nan",'NGC2419':"nan",'NGC2835':"nan"}
star_xc = { 'bsc_107': 412, 'bsc205': 406, 'bsc182': 425 , 'HD250092':565, 'a_And':565, 'GSC02500':335, 'HR928': 685, 'TYC1927':328, 'Orion': "nan" , 'Neptuno': "nan", 'Jupiter': "nan", 'M15': "nan" , 'Urano': "nan",'Crab': "nan",'NGC4988':"nan",'NGC4535':"nan",'NGC4147':"nan",'NGC2419':"nan",'NGC2835':"nan"}
star_yc = { 'bsc_107': 457, 'bsc205': 517, 'bsc182': 582 , 'HD250092':528, 'a_And':515, 'GSC02500':526, 'HR928': 537, 'TYC1927':425, 'Orion': "nan" , 'Neptuno': "nan", 'Jupiter': "nan", 'M15': "nan", 'Urano': "nan",'Crab': "nan" ,'NGC4988':"nan",'NGC4535':"nan",'NGC4147':"nan",'NGC2419':"nan",'NGC2835':"nan"}

i_ext=0.062
v_ext=0.124

for i in range (0,n):
  hdulist = pyfits.open(pathraw+str(data[i]))
  raw= hdulist[0].data
  prihdr = hdulist[0].header
  temp=prihdr['HIERARCH TEMPMEASURED']
  air=prihdr['AIRMASS']
  nimag=prihdr['HIERARCH CUBELENGTH']
  objeto=prihdr['OBJNAME']
  i_flux = i_obj[objeto]
  v_flux = v_obj[objeto]
  yc_obj = star_xc[objeto]
  xc_obj = star_yc[objeto]
  time=prihdr['EXPTIME']
  #time=re.sub(r'.*_(.*)g1m.*', r'\1', data[i])
  subname=re.sub(r'.*_(.*).fits', r'\1', data[i])
  #gain=prihdr['GAINMUL']
  reduc=(raw-bias)/flat
  hdu = pyfits.PrimaryHDU(reduc)
  hdulist = pyfits.HDUList([hdu])
  #hdulist.writeto(pathproc+"red_"+objeto+"_"+str(int(time))+"_m_"+str(m)+"_"+str(subname)+".fits")
  hdulist.close()
  npix_data = zeros(nimag)
  counts_data = zeros(nimag)
  npix_noise = zeros(nimag)
  counts_noise = zeros(nimag)
  if i_flux != "nan": #Quitamos los que no son estrellas (Orion....)
      i_flux_ex=i_flux+(air*i_ext) #Suma extincion
      v_flux_ex=v_flux+(air*v_ext)
      for j in range(0,nimag):
	img=reduc[j]
	maxval=img[xc_obj-ran:xc_obj+ran,yc_obj-ran:yc_obj+ran].max()
	[yc,xc]=where(reduc[j] == maxval)
	xi=xc-(outer*1.1)
	xf=xc+(outer*1.1)
	yi=yc-(outer*1.1)
	yf=yc+(outer*1.1)
	if len(yc)==1: #Evito dos maximos iguales en el rango LI
	  if xi>=0 and xf<=1024 and yi>=0 and yf<=1024: #Evito estrellas en el borde
	    for xx in range(xi,xf):
	      for yy in range(yi,yf):
		r = sqrt((xx-xc)**2+(yy-yc)**2)
		if r <=  inner:
		  counts_data[j]+=img[yy,xx]
		  npix_data[j]+=1
		else:
		  if inner<r<outer:
		    counts_noise[j]+=img[yy,xx]
		    npix_noise[j]+=1
	    if maxval>5*mean(img) and (counts_data[j]-(npix_data[j]*(counts_noise[j]/npix_noise[j])))>=0: #  sigma S/N = 5
	      lines=(data[i], objeto, i_flux, i_flux_ex, v_flux, v_flux_ex, psf,inner, outer, m, filt, air, temp, nimag, int(time), str(j+1)+"/"+str(nimag), yc_obj, xc_obj,xc[0], yc[0],maxval,counts_data[j],npix_data[j],counts_data[j]-(npix_data[j]*(counts_noise[j]/npix_noise[j])),counts_noise[j],npix_noise[j],counts_noise[j]/npix_noise[j])
	      print "\t".join([str(k) for k in lines])
	      print >> f,"\t".join([str(k) for k in lines])
	    else:
	      lines=(data[i], objeto, i_flux,i_flux_ex, v_flux, v_flux_ex, psf,inner, outer, m, filt, air, temp, nimag, int(time), str(j+1)+"/"+str(nimag), yc_obj, xc_obj,xc[0], yc[0],maxval,"nan","nan","nan","nan","nan","nan")
	      print "\t".join([str(k) for k in lines])
	      print >> f,"\t".join([str(k) for k in lines])
	  else:
	      lines=(data[i], objeto, i_flux,i_flux_ex, v_flux, v_flux_ex, psf,inner, outer, m, filt, air, temp, nimag, int(time), str(j+1)+"/"+str(nimag), yc_obj, xc_obj,xc[0], yc[0],maxval,"nan","nan","nan","nan","nan","nan")
	      print "\t".join([str(k) for k in lines])
	      print >> f,"\t".join([str(k) for k in lines])
	else:
	  lines=(data[i], objeto, i_flux, i_flux_ex, v_flux, v_flux_ex,psf,inner, outer, m, filt, air, temp, nimag, int(time), str(j+1)+"/"+str(nimag), yc_obj, xc_obj,xc[0], yc[0],maxval,"nan","nan","nan","nan","nan","nan")
	  print "\t".join([str(k) for k in lines])
	  print >> f,"\t".join([str(k) for k in lines])
      
      
      
      
    
        
      
      
      #plt.subplot(121)
      #plt.title(objeto+"  "+" m="+str(m)+" "+str(int(time))+"s "+str(filt)+" band ")
      #imgplot = plt.imshow(reduc[0])
      #imgplot.set_cmap('hot')
      #circle=plt.Circle((xc,yc),inner,color='w',fill=False)
      #circle2=plt.Circle((xc,yc),outer,color='w',fill=False)
      #fig = plt.gcf()
      #fig.gca().add_artist(circle)
      #fig.gca().add_artist(circle2)
      #plt.subplot(122)
      #imgplot = plt.imshow(reduc[0])
      #imgplot.set_cmap('hot')
      #plt.axis([xc-(outer*2),xc+(outer*2),yc-(outer*2),yc+(outer*2)])
      #circle=plt.Circle((xc,yc),inner,color='w',fill=False)
      #circle2=plt.Circle((xc,yc),outer,color='w',fill=False)
      #circle3=plt.Circle((xc,yc),1,color='b',fill=True)
      #fig2 = plt.gcf()
      #fig2.gca().add_artist(circle)
      #fig2.gca().add_artist(circle2)
      #fig2.gca().add_artist(circle3)
      #plt.show()
      #plt.title("["+str(xc[0])+","+str(yc[0])+"]"+" "+str(j+1)+"/"+str(nimag)+" "+str(maxval))
      #plt.figure()
      #plt.close()
      #plt.pause(1/1E9)
      #fig2.clear()
      #fig.clear()
