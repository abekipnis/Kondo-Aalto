#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jun 26 18:44:06 2021

@author: Whales
"""

#Rename the filename to your file and load the file

import numpy as np

import matplotlib as plt

import plotly.express as px
import pdb

import plotly.offline
# %pylab inline
import plotly.graph_objects as go
import pandas as pd
import gif
filename=r"S:\Createc_new\STMDATA\Ag\Small Kondo corrals\Ag 2021-08-13 2p5 nm radius\grid\Createc2_210814.214635.specgrid"
filename=r'S:\Createc_new\STMDATA\Ag\Small Kondo corrals\Ag 2021-08-16 2p5 nm radius empty\Createc2_210816.223358.specgrid'




# In this section the specgrid parameters are loaded and printed below

f = open(filename, "rb")

a = np.fromfile(f, dtype=np.uint32,count=256)

f.close

f = open(filename, "rb")

b = np.fromfile(f, dtype=np.float32,count=256)

f.close

version = a[0]

nx=a[1]

ny=a[2]

dx=a[3]

dy=a[4]

specxgrid=a[5]

specygrid=a[6]

vertpoints=a[7]
    
vertmandelay=a[8]

vertmangain=a[9]
# data type is different here, and if not int32 it would read wrong.
biasvoltage=b[10]

tunnelcurrent=b[11]

imagedatasize=a[12]

specgriddatasize=a[13]

specgridchan=a[14]

specgridchannelselectval=a[15]

specgriddatasize64=np.int64(a[17])

specgriddatasize64=(specgriddatasize64 << 32) + a[16]

xstart=a[18]

xend=a[19]

ystart=a[20]

yend=a[21]

specgridchannelselectval2=a[22]

specgridnx=a[23]

specgridny=a[24]

specgriddx=a[25]

specgriddy=a[26]

specgridcenterx=a[27]

specgridcentery=a[28]







# print(version,nx,ny,dx,dy,specxgrid,specygrid)

# print(vertpoints,vertmandelay,vertmangain)

# print(biasvoltage,tunnelcurrent)

# print(imagedatasize)

# print(specgriddatasize,specgridchan,specgridchannelselectval,specgriddatasize64)

# print(xstart,xend,ystart,yend)

# print(specgridchannelselectval2)

# print(specgridnx,specgridny);

# print(specgriddx,specgriddy);

# print(specgridcenterx,specgridcentery);
 



'''
# Load the Spectrum Data
# load the Spectrum data

# specvz contains the data of the biasvoltage-ramp and the z-ramp. 
#It is 2d array with dimensions specvz[vertpoints,specgridchan]

# the index runs from 0 .. x-1

# specdata contains the individual spectra. 

# the dimensions are specdata[specxpoints,specypoints,vertpoints,specgridchan]

'''


# f = open(filename, "r")

# b = np.fromfile(f, dtype=np.float32,count=256)
'位移256位？offset？fromfile 与之前读取过的位置有关系，之前读取过float32现在继续float32就会从之前的位置继续。'


count3=vertpoints*3

specvz = np.fromfile(f, dtype=np.float32,count=count3)
# print(specvz)
specvz3 = specvz.reshape(vertpoints,3)

''' data is the raw data of spectra, without any reshape or 
'''
data = np.fromfile(f, dtype=np.float32)


f.close

#specdata=np.arange(vertpoints*specgridchan*specxpoints*specypoints).reshape(specxgrid,specygrid,specgridchan,vertpoints)

# print(data.size)

' why? specgridny is not good ? because there might be some points out of scan frame?'
specgridnynew=int(data.size/vertpoints/specgridchan/specgridnx)

print(specgridnynew==specgridny)

# data=data.reshape(int(data.size/2),2)
      
# data2=np.resize(data,specgridnx*specgridnynew*vertpoints*specgridchan)


'need to check the spectra data if it matches the postion' 
specdata=data.reshape(specgridnx,specgridnynew,vertpoints,specgridchan)
'data reshape for each spectra point'

# print(specdata[0,1,:,:])
specdata_norm=np.empty_like(specdata)
'normalization of each single didv'
for i in range(xstart-1,xend): #python range，arange里重要的一点是不包括尾部，左闭右开
    for j in range(ystart-1,yend):
        specdata_norm[i,j,:,1] = specdata[i,j,:,1]/np.mean(specdata[i,j,:,1])
       

# i in range ((xstart-1),(xend-1)), j in range ((ystart-1),(yend-1)



# Display Spectroscopy Maps
# Display Spectroscopy maps 


# =============================================================================
# # '创建一个二维矩阵，对应每一个spectra点，每一个点赋值相对应能量的di值，然后再画出来这个二维矩阵'
# # image1=np.arange(specgridnx*specgridnynew).reshape(specgridnx,specgridnynew)
# 
# # for k in range(1,10):
# 
# #     image1=np.arange(specgridnx*specgridnynew).reshape(specgridnx,specgridnynew)
# 
# #     for i in range(0, specgridnx-1):
# 
# #        for j in range(0, specgridnynew-1):
# 
# #             image1[i,j]=specdata[i,j,100*k,1]
# 
# #     figure(figsize=(2,2), dpi=80)
# 
# #     plt.imshow(image1,cmap = cm.Greys_r)
# 
# # #    plt.savefig(filename+'.png',dpi=200)
# 
# =============================================================================

# '输出'
# import os
# filename = 'pseudo-nanonis.3ds'
# with open(filename, 'w') as output:
#     output.write('Grid dim='+ str(specgridnx) + 'x' + str(specgridny) + os.linesep)
#     output.write('Grid settings=0.000000E+0;0.000000E+0;1E-9;1E-9;00.000000E+0'+ os.linesep)
#     output.write('Sweep Signal="Bias (V)"'+ os.linesep)
#     output.write('Fixed parameters="Sweep Start;Sweep End"'+ os.linesep)
#     output.write('Experiment parameters='+ os.linesep)
#     output.write('# Parameters (4 byte)=2'+ os.linesep)
#     output.write('Experiment size (bytes)='+str(specgriddatasize)+ os.linesep)
#     output.write('Points='+ str(vertpoints) + os.linesep)
#     output.write('Channels="Lock-in Current"'+ os.linesep)
  
#     output.write(':HEADER_END:'+os.linesep)
    
# with open(filename, 'ab') as output:
#     output.write(bytes(13))
#     output.write(bytes(10))
#     for i in np.arange(xstart-1,xend): 
#         for j in np.arange(ystart-1,yend):

#             output.write(specvz3[0,0]+specvz3[0,-1])
#             output.write(np.ascontiguousarray(specdata_norm[i,j,:,1]))
            
            
# =============================================================================
# output = open('pseudo-nanonis.3ds', 'w')
# line = '\n'
#   
# #write()
# output.write()
# output.close()
# =============================================================================


# # '使用matplotlib来画图'
# k=300 #same as the k in plotly
# specmap=np.empty((xend,yend))
# for i in range(xstart-1,xend):
#     for j in range(ystart-1,yend):
#         specmap[i,j] = specdata[i,j,k,1]

# import matplotlib.pyplot as plt

# fig= plt.figure(figsize=(6,6))
# fig=plt.pcolormesh(specmap,cmap='YlOrRd_r')


# plt.contour(specmap)

# fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
# surf = ax.plot_surface(np.arange(xstart-1,xend-1),np.arange(ystart-1,yend-1), specmap[:,:], linewidth=0, antialiased=False)

# =============================================================================
# fig = px.imshow(specmap,zmin=0,zmax=50)
# 
# 
# import plotly.offline
# 
# plotly.offline.plot(fig)
# fig.write_html("file.html")
# 
# =============================================================================
 

# 3D dimension animation
# Import data
import time
import numpy as np

from skimage import io



specmapvol=np.empty((xend,yend,vertpoints))
for i in range(xstart-1,xend):
    for j in range(ystart-1,yend):
        specmapvol[i,j,:] = specdata[i,j,:,1]
        
vol = specmapvol
volume = vol.T
r, c = volume[0].shape

# Define frames
import plotly.graph_objects as go
nb_frames = vertpoints-1

fig = go.Figure(frames=[go.Frame(data=go.Surface(
    z=(k * 0.1) * np.ones((r, c)), #z=(k * 0.1) is from bottom to top, z=(nb_frames/10 - k * 0.1) top to bottom- depend on where the energy starts
    surfacecolor=np.flipud(volume[k]),#if z using k*0.1, then here is volume[k],rather than volume[nb_frames - k]
    #cmin=0.025*np.percentile(volume[k].flatten(),2.5), cmax=0.90*np.percentile(volume[k].flatten(),90)
    ),
    name=str(k) # you need to name the frame for the animation to behave properly
    )
    for k in range(vertpoints)])#nb_frames+1 or vertpoints because range doesn't include the right side

    
    
    
# Add data to be displayed before animation starts
fig.add_trace(go.Surface(
    z=nb_frames/10 * np.ones((r, c)),
    surfacecolor=np.flipud(volume[nb_frames]),
    colorscale='ylorrd_r',
    #cmin=0.025*np.percentile(volume[nb_frames].flatten(),2.5), cmax=0.90*np.percentile(volume[nb_frames].flatten(),90),
    colorbar=dict(thickness=20, ticklen=4)
    ))

# 每一帧之间的过渡？
def frame_args(duration):
    return {
            "frame": {"duration": duration},
            "mode": "immediate",
            "fromcurrent": True,
            "transition": {"duration": duration, "easing": "linear"},
        }

# 滑块？
sliders = [
            {
                "pad": {"b": 10, "t": 60},
                "len": 0.9,
                "x": 0.1,
                "y": 0,
                "steps": [
                    {
                        "args": [[f.name], frame_args(0)],
                        "label": str(k)+': '+str(specvz3[k,0])+' meV',#str(nb_frames-k)if top to bottom
                        "method": "animate",
                    }
                    for k, f in enumerate(fig.frames)
                ],
            }
        ]

# Layout 滑动条与页面显示？
fig.update_layout(
          title='dI/dV map from '+str(specvz3[0,0])+' meV to '+str(specvz3[-1,0])+' meV', 
          width=600,
          height=600,
          scene=dict(
                    zaxis=dict(range=[-0.1, (nb_frames/10+0.1)], autorange=False),
                    aspectratio=dict(x=1, y=1, z=1),
                    ),
          updatemenus = [
            {
                "buttons": [
                    {
                        "args": [None, frame_args(50)],
                        "label": "&#9654;", # play symbol
                        "method": "animate",
                    },
                    {
                        "args": [[None], frame_args(0)],
                        "label": "&#9724;", # pause symbol
                        "method": "animate",
                    },
                ],
                "direction": "left",
                "pad": {"r": 10, "t": 70},
                "type": "buttons",
                "x": 0.1,
                "y": 0,
            }
          ],
          sliders=sliders
          
)
import plotly.offline

plotly.offline.plot(fig)
# fig.write_html("file2.html")

# pdb.set_trace()
# #%%

# # Display individual Spectra
# #Display individual spectra

# c=specdata[specgridnx-1,specgridny-1,:,1]*10.0/512000


# for i in range(len(specvz3[0,:])):
    # d=specvz3[:,i];
    # plt.plot(d,c)
    # plt.savefig(r'S:\Createc_new\STMDATA\Ag\Small Kondo corrals\Ag 2021-08-13 2p5 nm radius\grid\%d.png' %(i),dpi=200)
# Pandas DataFrame with random data

# Gif function definition
# @gif.frame
# def plot(i):
    # fig = go.Figure()
    # fig.add_trace(go.Surface(
        # z=nb_frames/10 * np.ones((r, c)),
        # surfacecolor=np.flipud(volume[i]),
        # colorscale='ylorrd_r',
        # #cmin=0.025*np.percentile(volume[nb_frames].flatten(),2.5), cmax=0.90*np.percentile(volume[nb_frames].flatten(),90),
        # colorbar=dict(thickness=20, ticklen=4)
    # ))
    # fig.update_layout(width=500, height=300)
    # return fig

# # Construct list of frames
# frames = []
# for i in range(len(volume)):
    # frame = plot(i)
    # frames.append(frame)

# # Save gif from frames with a specific duration for each frame in ms
# gif.save(frames, r'S:\Createc_new\STMDATA\Ag\Small Kondo corrals\Ag 2021-08-13 2p5 nm radius\grid\2p5nm pm 20mV grid.gif', duration=10)