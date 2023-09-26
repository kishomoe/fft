import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.pylab import mpl
 
mpl.rcParams['font.sans-serif'] = ['SimSun']
mpl.rcParams['axes.unicode_minus']=False
#plt.rcParams['figure.dpi'] = 300
plt.rcParams['figure.constrained_layout.use'] = True

fig = plt.figure()
fixed_width_mm = 340
fixed_height_mm = 170
fixed_width_inch = fixed_width_mm / 25.4
fixed_height_inch = fixed_height_mm / 25.4
fig.set_size_inches(fixed_width_inch, fixed_height_inch)

ax1 = fig.add_subplot(121)
ax2 = fig.add_subplot(122)

data = pd.read_csv('Bx.csv')
x = data['Distance [mm]']
y = data['bx []']
ax1.plot(x,y,color='r',label='orgin',ls='--')

fft_y=np.fft.fft(y)
N=fft_y.shape[0]

orders=11
order=orders+1

fundamental=np.abs(fft_y[1])/N*2
harmonics=np.abs(fft_y[2:order])/N*2
thd=np.sqrt(np.sum(harmonics**2))/fundamental
thd_per = round(thd*100, 2)

for i in range(order):
	new_fft_y=np.zeros_like(fft_y)
	new_fft_y[i]=fft_y[i]
	ifft_y=np.fft.ifft(new_fft_y).real*2
	ax1.plot(x,ifft_y,label='order='+str(i))
	normalization_y=np.abs([fft_y[i]])/N*2
	ax2.bar(i,normalization_y)

ax1.legend()
ax1.set_title('各次谐波分解示意图')
ax1.set_xlabel('长度(mm)')
ax1.set_ylabel('各次谐波幅值大小')
ax1.set_xlim(x[0],x[len(x)-1])

ax2.set_title('各次谐波FFT分解柱状图 (THD='+str(thd_per)+'%)')
ax2.set_xlabel('谐波次数')
ax2.set_ylabel('各次谐波幅值大小')
ax2.set_xlim(0,order)
ax2.set_xticks(np.arange(0,order,1))

plt.show()