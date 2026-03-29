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
plt.xlim(0,0.1)

#data
data = np.loadtxt('ueg-502-1021-20.log',skiprows=1)
q = np.transpose(data)[0]
bmp2 = np.transpose(data)[1]
mp2 = np.transpose(data)[2]
ccd = np.transpose(data)[3]
rpa = np.transpose(data)[4]

# plot
plt.scatter(q,mp2,label='mp2',marker="o",color='Black',s=Size)
plt.scatter(q,bmp2,label='bmp2',marker="o",color='Red',s=Size)
plt.scatter(q,ccd,label='ccd',marker="o",color='Blue',s=Size)
plt.scatter(q,rpa,label='frpa',marker="o",color='Green',s=Size)

plt.legend()
plt.savefig("ueg-502-1021-20.pdf",format='pdf')


