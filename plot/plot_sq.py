import numpy as np
import matplotlib.pyplot as plt
# setting
plt.style.use([
    'yeplot.yeplot',      # overall style
    'yeplot.onecolumn',   # one-column format
    'yeplot.times',       # use times font
])
fsize = 8
Size = 2
plt.xticks(fontsize=fsize)
plt.yticks(fontsize=fsize)
plt.grid(linestyle = '--', linewidth = 0.5)
plt.ylabel(r'S(G)',fontsize=fsize)
plt.xlabel(r'$|q|$ ($\AA$)',fontsize=fsize)
#plt.xlim(0,400)

#data
mp2 = np.loadtxt('Li444.mp2',skiprows=1)
bmp2 = np.loadtxt('Li444.bmp2',skiprows=1)
r = np.transpose(bmp2)[0]
sqd = np.transpose(bmp2)[1]
sqx = np.transpose(bmp2)[2]
sq = sqd + sqx

# plot
plt.scatter(r,sqd,label='direct',marker="o",color='Red',s=Size)
plt.scatter(r,sqx,label='exchange',marker="o",color='Blue',s=Size)
plt.scatter(r,sq,label='S(G)',marker="o",color='Green',s=Size)

plt.legend()
plt.savefig("Li444-bmp2.pdf",format='pdf')


