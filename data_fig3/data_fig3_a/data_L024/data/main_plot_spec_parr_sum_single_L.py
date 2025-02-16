#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 11 20:21:42 2021

@author: maobinbin

Plot the spectrum of the system calculated by ED

compare the results of rho and |rho|
"""
import numpy as np
import warnings
import scipy.sparse as sp
from scipy.sparse.linalg import eigsh
import matplotlib.pyplot as plt
warnings.filterwarnings('ignore')
from matplotlib.ticker import MultipleLocator,LinearLocator, FormatStrFormatter,FuncFormatter

ii=0+1j

para_Jr=173200
data_in=np.loadtxt('input.in')
L_all = [int(data_in[0])]
para_parr_dir=1
para_parr_max=50
para_nbins=10
para_isteps=10000
para_isteps=100000

import matplotlib.pyplot as plt
from matplotlib import rc
plt.cla
plt.rcParams['xtick.major.pad']='5'
plt.rcParams['ytick.major.pad']='3'
plt.rcParams['axes.linewidth'] = '1'
rc('text', usetex=True)
value_fontsize=20
value_marker_size=7
value_line_width=1


color_list=['purple','olive','green','blue','magenta','olive']
marker_list=['o','s','<','^','>','*','p','d']
line_list=['-',':','-',':','-.']
line_list=['^r','^b','^k','^c','^g','^m','or','ob','ok','oc','og','om','*r','*b','*k','*c','*g','*m','dr','db','dk','dc','dg','dm','sr','sb','sk','sc','sg','sm']
#line_list=['sr','ob','dc','<g','xm','-sr','-or','-ob','-ok','-oc','-og','-om']

  




limit_x_min=-1.01*np.pi
limit_x_max=1.01*np.pi


limit_y_min=2.8
limit_y_max=3.4

nrows =1
ncols = 1
fig, axs = plt.subplots(nrows,ncols,sharex=True,figsize=(5,5))
plt.subplots_adjust(hspace=0.01,wspace=0.4)


#***************************************
# plot the error of numerical results
#***************************************
ax2 = plt.subplot(1,1,1)
value_markeredgewidth=1
plt.cla()

#entropy_data=np.zeros([np.size(L_all),para_Nlel])

max_val_all=[]
for iL, L in enumerate(L_all):
    if L==18:
        para_parr_dir=2
    else:
        para_parr_dir=1    
#    for iSz in range(int(L/2)+1,int(L/2)+2):
    spec_data=[]
    spec_data_all=[]
#    for iSz in range(12):
    for iSz in range(int(L/2+1)):
        for k in range(int(L/2+1)):
            file_name = './eigs_Sz%d_k%d.dat' %(iSz,k)
            spec_data=np.loadtxt(file_name)
            spec_data_all=np.hstack((spec_data_all,spec_data))
            
            
        
    spec_data_all = np.reshape(spec_data_all,-1)
#    spec_data_all =np.sort(spec_data_all)
    max_val=max(spec_data_all)
#    max_val=spec_data_all[-3]
    max_val_all.append(max_val)
#print('max_val',max_val)
print('ED L=%d, Jr=%d' %(L,para_Jr))
#print(np.sort(spec_data_all))


max_val = 1
for iL, L in enumerate(L_all):
    if L==18:
        para_parr_dir=2
    else:
        para_parr_dir=1    
        
    max_val=max_val_all[iL]

#    for iSz in range(int(L/2)+1,int(L/2)+2):
#    for iSz in range(int(L/2),int(L/2)+1):
#    ax2.plot([0],[10],line_list[iL],linewidth=value_line_width,markersize=value_marker_size,markeredgewidth=value_markeredgewidth,fillstyle='none',label=r'$L=%d$' %(L))

#    for iSz in range(15,16):
#    for iSz in range(L+1):
    count_line=0
#    for iSz in range(1):
    for iSz in range(int(L/2+1)):
#    for iiSz,iSz in enumerate([0,1,2,3]):
#        for k_aux in range(1):
        for k_aux  in range(int(L/2+1)):
#        for k_aux in range(-L,L,1):    
            k=np.mod(k_aux,L)
            
            spec_data=[]
            file_name = './eigs_Sz%d_k%d.dat' %(iSz,k)           
            spec_data=np.loadtxt(file_name)
            print('(%d, %d)' %(k_aux,iSz),-np.log(spec_data))
#            count_line = iSz
            if k_aux == 0:
                ax2.plot([0],[-1],line_list[count_line],linewidth=value_line_width,markersize=value_marker_size,markeredgewidth=value_markeredgewidth,fillstyle='none',label=r'$Sz=%d$' %(iSz))#
                
            if (np.size(spec_data)==1 and spec_data !=0):
                ax2.plot((k_aux)/L*2*np.pi,-np.log(spec_data)+np.log(max_val),line_list[count_line],linewidth=value_line_width,markersize=value_marker_size,markeredgewidth=value_markeredgewidth,fillstyle='none')#
                
                if k_aux!=0:
                    ax2.plot((-k_aux)/L*2*np.pi,-np.log(spec_data)+np.log(max_val),line_list[count_line],linewidth=value_line_width,markersize=value_marker_size,markeredgewidth=value_markeredgewidth,fillstyle='none')#

                '''
                if iSz==int(L/2):
                    ax2.plot((k_aux)/L*2,-np.log(abs(spec_data))+np.log(max_val),'xr',linewidth=value_line_width,markersize=value_marker_size,markeredgewidth=value_markeredgewidth,fillstyle='none')#
                if iSz==int(L/2)+1:
                    ax2.plot((k_aux)/L*2,-np.log(abs(spec_data))+np.log(max_val),'+b',linewidth=value_line_width,markersize=value_marker_size,markeredgewidth=value_markeredgewidth,fillstyle='none')#
                '''
            elif (np.size(spec_data)>1 and np.ceil(sum(abs(spec_data)))>0):
                spec_data=np.reshape(spec_data,[1,np.size(spec_data)])
                data_x=[]
                data_x=np.ones([1,np.size(spec_data)])*(k_aux)/L*2*np.pi
                ax2.plot(data_x[0,:],-np.log(spec_data[0,:])+np.log(max_val),line_list[count_line],linewidth=value_line_width,markersize=value_marker_size,markeredgewidth=value_markeredgewidth,fillstyle='none')#
                if k_aux!=0:
                    data_x=[]
                    data_x=np.ones([1,np.size(spec_data)])*(-k_aux)/L*2*np.pi
                    ax2.plot(data_x[0,:],-np.log(spec_data[0,:])+np.log(max_val),line_list[count_line],linewidth=value_line_width,markersize=value_marker_size,markeredgewidth=value_markeredgewidth,fillstyle='none')#
                '''
                if iSz==int(L/2):
                    ax2.plot(data_x[0,:],-np.log(abs(spec_data[0,:]))+np.log(max_val),'xr',linewidth=value_line_width,markersize=value_marker_size,markeredgewidth=value_markeredgewidth,fillstyle='none')#

#                    ax2.plot(data_x[0,0],min(-np.log(abs(spec_data[0,:]))+np.log(max_val)),'xr',linewidth=value_line_width,markersize=value_marker_size,markeredgewidth=value_markeredgewidth,fillstyle='none')#
                
                if iSz==int(L/2)+1:
                    ax2.plot(data_x[0,:],-np.log(abs(spec_data[0,:]))+np.log(max_val),'+b',linewidth=value_line_width,markersize=value_marker_size,markeredgewidth=value_markeredgewidth,fillstyle='none')#)#
                '''
#                    ax2.plot(data_x[0,0],min(-np.log(abs(spec_data[0,:]))+np.log(max_val)),'+k',linewidth=value_line_width,markersize=value_marker_size,markeredgewidth=value_markeredgewidth,fillstyle='none')#
        count_line=count_line+1
#for iL, L in enumerate(L_all):
#    ax2.plot([0],[10],line_list[iL],linewidth=value_line_width,markersize=value_marker_size,markeredgewidth=value_markeredgewidth,fillstyle='none',label=r'$L=%d$' %(L))



#y_coeff = 0.25
#x=[(100-i)*np.pi/100 for i in range(200)]
#y=[]
#for ix,x_val in enumerate(x):
#    y.append(y_coeff -y_coeff *np.cos(x_val/1.2))

#for iy,y_val in enumerate(y):
#    y[iy]=y[iy]+0.08
#for ix,x_val in enumerate(x):
#    x[ix]=x[ix]-0.07


#ax2.plot(x,y,'-k',linewidth=value_line_width,markersize=value_marker_size,markeredgewidth=value_markeredgewidth,fillstyle='none')



ax2.set_xlabel(r'$k$', fontsize=value_fontsize)
ax2.set_ylabel(r'$\xi_{\alpha}-\xi_0$', fontsize=value_fontsize)#,fontproperties=font_name
for label in ax2.xaxis.get_ticklabels():
    label.set_fontsize(value_fontsize)
for label in ax2.yaxis.get_ticklabels():
    label.set_fontsize(value_fontsize)

ax2.set_xlim([-np.pi,np.pi])


ax2.set_ylim([0,2])
#ax2.set_xticks(np.arange(-1,1.5,0.5))
#ax2.set_yticks(np.arange(0,2.5,0.5))

#ax2.set_xticks = [-2*np.pi, -3*np.pi/2, -np.pi, -np.pi/2, 0, np.pi/2, np.pi, 3*np.pi/2, 2*np.pi]
#ax2.set_set_xticklabes = ['-2π', '-3π/2', '-π', '-π/2', 0, 'π/2', 'π', '3π/2', '2π' ]
#ax2.set_xticks(np.arange(-1,1.5,0.5))
#ax2.xticklabes = ['-2π', '-3π/2', '-π', '-π/2', 0, 'π/2', 'π', '3π/2', '2π' ]

'''
def pi_formatter(x, pos):
    """
    比较罗嗦地将数值转换为以pi/4为单位的刻度文本
    """
    m = np.round(x / (np.pi/4))
    n = 4
    if m%2==0: m, n = m/2, n/2
    if m%2==0: m, n = m/2, n/2
    if m == 0:
        return "0"
    if m == 1 and n == 1:
        return "$\pi$"
    if n == 1:
        return r"$%d \pi$" % m
    if m == 1:
        return r"$\frac{\pi}{%d}$" % n
    return r"$\frac{%d \pi}{%d}$" % (m,n)
'''
def pi_formatter(x, pos):
    """
    比较罗嗦地将数值转换为以pi/4为单位的刻度文本
    """
#    print ('x=',x)
#    print ('pos=',pos)    
    m = np.round(x/ (np.pi/8) )
    n = 8
    if m%2==0: m, n = m/2, n/2
    if m%2==0: m, n = m/2, n/2
    if m == 0:
        return "0"
    if m == 1 and n == 1:
        return "$\pi$"
    if n == 1:
#        print ('188, m=%d, n=%d' %(m,n))    
        return r"$%d \pi$" % m

    if m == 1:
#        print ('191, m=%d, n=%d' %(m,n))
        return r"$\frac{\pi}{%d}$" % n
    if m == -1:
#        print ('191, m=%d, n=%d' %(m,n))
        return r"$-\frac{\pi}{%d}$" % n

    if m == -2:
#        print ('191, m=%d, n=%d' %(m,n))
        return r"$-\pi$"
    if m == 2:
#        print ('191, m=%d, n=%d' %(m,n))
        return r"$\pi$"


#    print ('194, m=%d, n=%d' %(m,n))    
    return r"$\frac{%d \pi}{%d}$" % (m,n)


ax2.xaxis.set_major_locator( MultipleLocator(np.pi/2) )
ax2.xaxis.set_major_formatter( FuncFormatter( pi_formatter ) )

#ax2.yaxis.set_major_locator( MultipleLocator(np.pi/4) )
#ax2.yaxis.set_major_formatter( FuncFormatter( pi_formatter ) )

#ax2.xaxis.set_major_locator( [-1,-0.5,0,0.5,1])
#ax2.xaxis.set_major_formatter( [r'$-\pi$',r'$\frac{-\pi}{2}$',r'$0$',r'$\frac{\pi}{2}$',r'$\pi$'])




plt.legend(loc='center right',ncol=3,fancybox=True,shadow=False,numpoints=1,frameon=False,bbox_to_anchor=(1,1.1))#,
leg=plt.gca().get_legend()
ltext=leg.get_texts()
plt.setp(ltext,fontsize=value_fontsize-2)


pdf_file_name = 'fig_AFM_L%03d_Jr%06d_nbins_%05d.pdf' %(L,para_Jr,para_nbins)
plt.savefig(pdf_file_name, dpi=None, facecolor='w', edgecolor='w', orientation='portrait', format="pdf", bbox_inches='tight')

