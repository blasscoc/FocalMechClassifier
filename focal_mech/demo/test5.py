from numpy import array, argmin, arange, argmax, deg2rad, dot
from scipy.linalg import norm

import matplotlib.pyplot as plt
import mplstereonet

from focal_mech.lib.sph_harm import WignerD2
from focal_mech.lib.classify_mechanism import classify, translate_to_sphharm
from focal_mech.io.read_hash import read_demo, read_hash_solutions
from focal_mech.util.hash_routines import hash_to_classifier


hash_solns = read_hash_solutions("example1.out")

# we want solutions that are symetric
polarity_data = read_demo("north1.phase", "scsn.reverse", reverse=True)
inputs = hash_to_classifier(polarity_data, parity=1)

event = 3146815

result = classify(*inputs[event], kernel_degree=2)
Alm = translate_to_sphharm(*result, kernel_degree=2)

alm = array([Alm[0,0], 
             Alm[1,-1], Alm[1,0], Alm[1,1], 
             Alm[2,-2], Alm[2,-1], Alm[2,0], Alm[2,1], Alm[2,2]])

corr = []
for event in polarity_data.keys():
    
    if event == 3146815:
        continue
    
    result = classify(*inputs[event], kernel_degree=2)
    Blm = translate_to_sphharm(*result, kernel_degree=2)

    blm = array([Blm[0,0], 
                 Blm[1,-1], Blm[1,0], Blm[1,1], 
	         Blm[2,-2], Blm[2,-1], Blm[2,0], Blm[2,1], Blm[2,2]])
             
    c = norm(alm.conjugate().dot(blm))/(norm(alm)*norm(blm))

    print event, c
    corr.append(c)

blm = array([Alm[0,0], 
             Alm[1,-1], Alm[1,0], Alm[1,1], 
             Alm[2,-2], Alm[2,-1], Alm[2,0], Alm[2,1], Alm[2,2]])
                
D = WignerD2(deg2rad(90), deg2rad(0), deg2rad(0) )
blm[4:] = dot(D, blm[4:]) 

# rotate event by 90-degreees
corr.append(norm(alm.conjugate().dot(blm))/(norm(alm)*norm(blm)))

# \alpha makes sense for sandstone or granite, 4\sqrt(5) is the other term
blm = [4.7, 0,0,0, 0,0,9,0,0]

corr.append(norm(alm.conjugate().dot(blm))/(norm(alm)*norm(blm)))

events = [i for i in polarity_data.keys() if i != 3146815]

events.append("90-degree")
events.append("Tensile")

N = len(events)

fig = plt.figure(facecolor="white", figsize=(20,30))
ax = fig.add_subplot(111)
rects1 = ax.bar([i+0.5 for i in range(N)], corr, 0.99, align = 'center')

rects1[argmin(corr[:-2])].set_color('r')
rects1[argmax(corr[:-2])].set_color('r')

rects1[-1].set_color('c')
rects1[-2].set_color('g')


ax.set_yticks(arange(0,1,0.1))      
ax.set_xticks(range(0,N))     
ax.set_xticklabels(events, rotation=25, fontsize=18 )
ax.yaxis.grid(True)

for tick in ax.yaxis.get_major_ticks():
    tick.label.set_fontsize(18) 
    
plt.xlim([0,N])

plt.ylabel("Correlation Score", fontsize=24)

plt.show()
    
    