#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 19 19:57:51 2019

@author: binbin
"""
import numpy as np
from scipy import sparse
import scipy.sparse.linalg as sparse_la
import copy
import warnings
from scipy.sparse.linalg import eigsh
from matplotlib.ticker import MultipleLocator,LinearLocator, FormatStrFormatter,FuncFormatter
import datetime
warnings.filterwarnings('ignore')
import warnings
import scipy.sparse as sp
from scipy.sparse.linalg import eigsh
import os
import subprocess
import time



class model_para:
    """

    DESCRIPTION: use "model_para" class to define all necessary global parameters

    """
    def __init__(self, L=4, Jx=1.0, Jz=1.0):
        self.L = L
        self.Jx = Jx
        self.Jz = Jz
        self.n_eig = 20
        print('L = %d, Jx = %.2f, Jz = %.2f' %(L, Jx, Jz))


ii=0+1j

data_in=np.loadtxt('input.in')

para_Jr=173200
L_in = int(data_in[0])
L_all=[L_in]
L=L_in
para_parr_start=1
para_parr_max=100
para_nbins=10
#L_all=[14]
para_msteps=10000
para_isteps=100000


starttime = datetime.datetime.now()


eta = 0.1
deg_err = 1e-10
para_aux = 1
Jx_in = para_aux 
Jz_in = para_aux 

site_dim = 2

sp_den='s'




eqs = 'wu'
eqs = 'OGk'
eqs = 'OGInvk'

    



##***************************
##  plot
import matplotlib.pyplot as plt
from matplotlib import rc
plt.cla
plt.rcParams['xtick.major.pad']='5'
plt.rcParams['ytick.major.pad']='3'
plt.rcParams['axes.linewidth'] = '1'
rc('text', usetex=True)
value_fontsize=20
value_marker_size=7
value_line_width=1.2


color_list=['purple','olive','green','blue','magenta','olive']
marker_list=['o','s','<','^','>','*','p','d']
line_list=['-',':','-',':','-.']
line_list=['-or','-ob','-ok','-oc','-og','-om','-sr','-or','-ob','-ok','-oc','-og','-om']
nrows =1
ncols = 1
fig, axs = plt.subplots(nrows,ncols,sharex=True,figsize=(5,4))
plt.subplots_adjust(hspace=0.01,wspace=0.4)

#***************************************
# plot the error of numerical results
#***************************************
ax2 = plt.subplot(1,1,1)
value_markeredgewidth=2
#print ('E0 =',np.min(eigs_all))
#print ('eigs_plot =',eigs_plot-E0)



#data_x=data[:,0]-energy
#data_y=-data[:,3]
#ax2.plot(data_x,data_y,'-r',linewidth=value_line_width,markersize=value_marker_size,markeredgewidth=value_markeredgewidth,fillstyle='none')#


x=[]
file_name='data_plot_k_all_L%d.dat' %(L)
x_aux=np.loadtxt(file_name)
#x=x-L/2


x=[]
for ix,x_val in enumerate(x_aux):
    x.append((x_val-L/2)/L*2*np.pi)


y=[]
file_name='data_plot_omega_all_L%d.dat' %(L)
y=np.loadtxt(file_name)



z=[]
file_name='data_plot_spec_func_L%d.dat' %(L)
z=np.loadtxt(file_name)



print(np.shape(z))

[row,col] = np.shape(z)
omeg_val_all=[]
for i in range(row):
    max_value = np.max(z[i,:])
    max_ind = np.where(z[i,:]==max_value)
    omeg_val_all.append(y[max_ind][0])
    print('---- max_value',max_value)
    print('max_ind',max_ind)
    print('max_val',z[i,max_ind])

y = y-omeg_val_all[int(L/2)]
omeg_val_all=omeg_val_all-omeg_val_all[int(L/2)]


xx,yy = np.meshgrid(x,y)

#fig, ax = plt.subplots()
c = ax2.pcolormesh(xx,yy,z.T,cmap='jet')#hot_r,jet
fig.colorbar(c, ax=ax2)


#ax2.plot(x,omeg_val_all,'or',linewidth=value_line_width,markersize=value_marker_size,markeredgewidth=value_markeredgewidth,fillstyle='none')#


x_fit = x[int(L/2):int(L/2)+2]
x=copy.deepcopy(x_fit)
y = copy.deepcopy(omeg_val_all[int(L/2):int(L/2)+2])

func = np.polyfit(x, y, 1)
xn = np.linspace(0, 3, 100)
yn = np.poly1d(func)
#plt.plot(xn, yn(xn), x, y, 'ow')
plt.plot(xn, yn(xn), '-w')
print('polyfit ', func)



plt.text(-np.pi/4*3, 2, r'$v \approx %.4f$' %(func[0]), fontsize=value_fontsize,color='white')

#plt.text(-np.pi/4*3, 2, r'$d_1 \approx %.4f$' %(func[0]*np.pi/6.0), fontsize=value_fontsize,color='white')

file_name = '../data_v_spec_func.dat'
np.savetxt(file_name,[func[0]])


#ax2.text(2,11,r'$%d$ sites' %(L),fontsize=20)

ax2.set_xlabel(r'$k$', fontsize=value_fontsize)
ax2.set_ylabel(r'$\omega$', fontsize=value_fontsize)#,fontproperties=font_name

for label in ax2.xaxis.get_ticklabels():
    label.set_fontsize(value_fontsize)
for label in ax2.yaxis.get_ticklabels():
    label.set_fontsize(value_fontsize)

'''
plt.legend(loc='center right',ncol=1,fancybox=True,shadow=False,numpoints=1,frameon=False)#bbox_to_anchor=(0.8,0.97),
leg=plt.gca().get_legend()
ltext=leg.get_texts()
plt.setp(ltext,fontsize=value_fontsize-2)
'''


#ax2.set_xscale("log")
ax2.set_xlim(-np.pi, np.pi)
ax2.set_ylim([0,3])


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
        return "$0$"
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


#plt.xticks(x, calendar.month_name[1:13],color='blue',rotation=60)#此处locs参数与X值数组相同

#plt.xticks([])
#xx=range(0, int(L+L/4), int(L/4))
#plt.xticks(xx,[r'$0$',r'$\frac{\pi}{2}$',r'$\pi$',r'$\frac{3\pi}{2}$',r'$2\pi$'])


pdf_file_name = 'fig_SpecFunc_ladder_L_%03d_J_%.2f_%s_slope.pdf' %(L,para_aux,eqs)
plt.savefig(pdf_file_name, dpi=None, facecolor='w', edgecolor='w', orientation='portrait', format="pdf", bbox_inches='tight')

    
    
endtime = datetime.datetime.now()
print('time total:',endtime -starttime)

