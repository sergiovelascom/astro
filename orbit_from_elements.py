#!/usr/bin/env python
# -*- coding: utf-8 -*-


__doc__ = "Pinta la orbita de una binaria a partir de sus elementos orbitales"

from numpy import *
import matplotlib.pyplot as plt
from matplotlib import rc


d2rad = pi/180.
rad2d = 180./pi
   
def rt2xy(rho, theta):
    """Pasa de coordenadas rho y theta a (x,y) y rota 90º para tener N arriba"""
    x = rho*cos((theta+90)*(pi/180.))
    y = rho*sin((theta+90)*(pi/180.))
    return x, y

def draw_origin():
    # pinta el origen de coordenadas (posicion de la primaria)
    y0_line = plt.axhline(0, color='#cccccc', alpha=0.5, label='_nolegend_')
    x0_line = plt.axvline(0, color='#cccccc', alpha=0.5, label='_nolegend_')

def binary_position(t, P, T, e, a, i, O_node, o_peri):
   # Calcula la posicion de la compañera en coordenadas polares
   # Grados a radianes
   i = i*d2rad
   O_node = (O_node*d2rad)%(2*pi)
   o_peri = (o_peri*d2rad)%(2*pi)
 
   # Anomalia media
   M = ((2.0*pi)/P)*(t - T)  # radianes
    
   if M >2*pi: M = M - 2*pi 
   M=M%(2*pi)

   # Anomalia excentrica (1ra aproximacion)
   E0 = M  + e*sin(M) + (e**2/M) * sin(2.0*M)

   # Itero para calcular la anomalia excentrica
   for itera in range(15):
      M0 = E0 - e*sin(E0)
      E0 = E0 + (M-M0)/(1-e*cos(E0))

   # Anomalia verdadera
   true_anom = 2.0*arctan(sqrt((1+e)/(1-e))*tan(E0/2.0))
 
   #radius = (a*(1-e**2))/(1+e*cos(true_anom))
   radius = a*(1-e*cos(E0))

   theta = arctan( tan(true_anom + o_peri)*cos(i) ) + O_node
   rho = radius * (cos(true_anom + o_peri)/cos(theta - O_node))
   
   # revuelve rho ("), theta (grad)
   return rho, (theta*rad2d)%360.



# --- Calculo para nuestros elementos 

## New elements (srh)
P = 15.184596
T = 1992.297
a = 0.94600318    
e =  0.32667310
i =  76.950024  
O_node  = 322.77072 
o_peri= 191.80015 

# Calcula la posicion de la compañera a lo largo de un periodo
pos1 = [binary_position(time, P, T, e, a, i, O_node, o_peri) for time in arange(T,T+P+0.19,0.1)]
pos1 = array(pos1)

# coodenadas cartesinas
x1,y1 = rt2xy(pos1[:,0], pos1[:,1])


# --- Calculo para los elementos de Torres et al

# Torres elements
P = 15.643                  # Periodo (años)
T = 1992.297                # Paso por el periastro (años)
e = 0.295                  # Ecentricidad
a = 0.9257                  # Semieje mayor (seg de arco)
i = 103.0                  # Inclinacion de la orbita (grados)
O_node = 143.48              # Posicion del nodo ascendente (Omega, grados)
o_peri = 347.2            # Argumento del periastro (omega, grados)

# Calcula la posicion de la compañera a lo largo de un periodo
pos2 = [binary_position(time, P, T, e, a, i, O_node, o_peri) for time in arange(T,T+P+0.19,0.1)]
pos2 = array(pos2)

# coodenadas cartesinas
x2,y2 = rt2xy(pos2[:,0], pos2[:,1])


# --- Calculo para los elementos de Branham et al

# Branham elements
P = 15.86                  # Periodo (años)
T = 1994.67                # Paso por el periastro (años)
e = 0.138                  # Ecentricidad
a = 1.067                  # Semieje mayor (seg de arco)
i = 81.92                  # Inclinacion de la orbita (grados)
O_node = 146.296            # Posicion del nodo ascendente (Omega, grados)
o_peri = 276.577            # Argumento del periastro (omega, grados)

# Calcula la posicion de la compañera a lo largo de un periodo
pos3 = [binary_position(time, P, T, e, a, i, O_node, o_peri) for time in arange(T,T+P+0.19,0.1)]
pos3 = array(pos3)

# coodenadas cartesinas
x3,y3 = rt2xy(pos3[:,0], pos3[:,1])


# --- Puntos para las medidas

# Medidas Lucky imaging (Fastcam & Astralux)
fc_data = array([
#[0.45, 160.9],    # TCS  Feb 2007
[0.424, 160.6],    # TCS  Feb 2007  - Nueva calibracion
[0.631, 145.9],  # TCS  Ene 2008 
[0.638, 143.08],  # CAHA Abr 2008
#[0.6218, 143.83],  # CAHA Abr 2008
[0.648, 142.6],  # NOT May  2008
#[0.654, 142.64],  # NOT May  2008 (4x)
[0.604, 136.1],   # CAHA Ene 2009 (4x) 
[0.579, 132.6],   # CAHA Marzo 2009 (x4)
#[0.5600, 132.483], # CAHA Mayo 2009 
[0.548, 131.5], # CAHA Mayo 2009 (x4)
#[0.540, 134.30] # WHT Sep  2009
[0.231, 78.7],   # TCS 24 Junio 2010 - 2010.478953


]) 

# Medidas speckle bidimensional
spk_obs = array([
[1.098, 322.7],
[0.220, 278.0],
[0.179, 248.],
[0.177, 246.],
[0.177, 226.],
[0.186, 210.5],
[0.43, 160.],
[0.343, 164.],
[0.374, 164.],
[0.436, 161.6],
[0.602, 135.8],
])

# Medidas con HST (FGS y FOS)
hst_obs = array([
[0.2359, 26.04, 1995.5733],
[0.416, 353.2, 1996.2923],
[0.4639, 350.53, 1996.4606],
[0.7344,338.44, 1997.4241],
[0.9245, 333.46, 1998.2510]
])


# Medidas visuales y fotomentricas
visual_obs = array([
  [   0.72, 321.7 ,  1938.34],  
  [   1.0,  312.  , 1938.38 ], 
  [   1.10, 321.3 , 1938.42 ], 
  [   0.81, 316.1 , 1939.18 ], 
  [   0.73, 311.4 , 1940.16 ], 
  [   0.29, 300.3 , 1941.18 ], 
  [   0.4 , 239.  , 1941.44 ], 
  [   0.15, 260.  , 1941.51 ], 
  [  0.45 , 126.4 , 1946.20 ], 
  [  0.70 , 335.6 , 1950.13 ], 
  [  0.80 , 325.0 , 1952.29 ], 
  [  0.69 , 326.3 , 1952.37 ], 
  [  0.51 , 294.2 , 1952.53 ], 
  [  0.68 , 323.8 , 1953.21 ],
  [  0.73 , 316.7 , 1954.15 ], 
  [  0.73 , 338.7 , 1955.38 ], 
  [  0.43 , 316.9 , 1956.17 ], 
  [  0.44 , 315.8 , 1956.33 ], 
  [  0.8  , 169.3 , 1960.30 ], 
  [  1.5  , 140.  , 1962.22 ], 
  [  0.37 , 353.1 , 1965.176], 
  [  0.2  ,  0.   , 1965.40],
  [  0.56 , 343.5 , 1966.29 ], 
  [  0.72 , 329.2 , 1967.18 ], 
  [  0.73 , 326.9 , 1968.20 ], 
  [  0.75 , 336.6 , 1968.34 ],
  [  0.77 , 328.6 , 1969.16 ],
  [  0.84 , 324.5 , 1970.32 ],
  [  0.73 , 323.1 , 1970.46 ],
  [  0.68 , 320.2 , 1971.30 ],
  [  0.35 , 125.8 , 1979.18 ],
  [  0.87 , 324.0 , 1984.50 ],
  [  0.86 , 330.8 , 1985.32 ],
  [  0.81 , 312.8 , 1987.72 ],
  [  0.30 , 304.9 , 1989.39 ],
  [  0.18 , 197.4 , 1991.17 ]
])

fcx, fcy = rt2xy(fc_data[:,0], fc_data[:,1])
spkx, spky = rt2xy(spk_obs[:,0], spk_obs[:,1])
hstx, hsty = rt2xy(hst_obs[:,0], hst_obs[:,1])
visx, visy = rt2xy(visual_obs[:,0], visual_obs[:,1])

plt.rc("font", size=12)

draw_origin() # Marca la posicion de la primaria

# Lineas de las orbitas
plt.plot(x1, y1, 'k-', lw=0.8, label='This work')
plt.plot(x2, y2, 'r:', lw=1.0, label='Torres et al.')
#plt.plot(x3, y3, 'g:', lw=0.8, label='Branham')

# Puntos de las observaciones
plt.plot(fcx, fcy, 'r.', label='Lucky Imaging')
plt.plot(spkx, spky, 'b^', label='Speckle')
plt.plot(hstx, hsty, 'gs', label='HST')
#plt.plot(visx, visy, 'm.', label='Visual/Photo')


plt.xlabel(r'$\Delta\alpha$ ($^{\prime\prime}$)', fontsize=18)
plt.ylabel(r'$\Delta\delta$ ($^{\prime\prime}$)', fontsize=18)

#plt.legend(loc='upper left')

plt.show()
plt.grid()
