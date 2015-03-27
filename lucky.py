#!/usr/bin/env python

'''Suma las imagenes de un cubo FITS eligiendo las mejores imagenes
'''


from numpy import *
import numpy
import pyfits, os
from scipy import ndimage
from jappLib import *
from string import zfill

# Configuracion
in_path = "/media/sda2/scratch/Fastcam/observaciones/CAHA-Enero2009/20090114/"
out_path = "/media/sda2/scratch/Fastcam/observaciones/CAHA-Enero2009/20090114/"
filename = 'Gl473_cube.fits'  #Cubo de imagenes de entrada
fname, fext = filename.rsplit('.',1)  # name an extesion of input file
dark_filename = 'Master_Dark_010g1m150_I.fits'
imagesPercent = 100  # Porcentaje de imagenes con el que nos quedamos
extract_frames = False    # save also individual shifted frames
extract_dark = False
#-------------------

images = pyfits.getdata(in_path + filename)
if extract_dark:
	dark = pyfits.getdata(in_path + dark_filename)
lista = array([])
lista.shape = (0,4)
NimagesTotal, sizeX, sizeY = images.shape  # Numero total de imagenes y tamanho de las imagenes

for i in range(0, NimagesTotal):
   img = images[i]
   maxVal = ndimage.maximum(img)
   posY, posX = ndimage.maximum_position(img)
   #posX, posY = myfits.centerMass(img)
   # Anhado a la lista los valores no. de imagen, valor de px max y su posicion
   newRow = array([[i,maxVal, posX, posY]])
   lista = concatenate((lista,newRow))

# Ordeno por indices los pixel mas brillantes 
indices = argsort(lista[:,1])
Nimages = int(NimagesTotal * (imagesPercent/100.))
index = indices[-Nimages:]

posY0, posX0 = ndimage.maximum_position(images[index[0]])

mask = boxMask(images[0].shape, c=(256,256), side=500)

C = zeros_like(images[0])      # imagen de inicio en blanco para cada frame
Ctotal = zeros_like(images[0]) # imagen de inicio en blanco para el promedio

Nimage = 0
for i in index:
	posY, posX = ndimage.maximum_position(images[i]*mask)
	dx, dy = int(posX0 - posX), int(posY0 - posY)
	if extract_dark:
		images[i] = images[i] - dark
	C = shift_array(images[i], dy, dx)
	Ctotal = Ctotal + C
	if extract_frames:
		out_file = fname[4:] + '-' + zfill(Nimage,4) + '-' + zfill(i,4) + '.' + fext
		pyfits.writeto(out_path + out_file, C)
		print "Guardando imagen %s recentrada: %s " %  (i, out_file)
	Nimage = Nimage+1
	
out_file = fname[4:] + '-mean.' + fext
pyfits.writeto(out_path + fname[4:] + '-mean.' + fext, Ctotal/Nimages) 
print "Guardando imagen promediada: ", out_file
