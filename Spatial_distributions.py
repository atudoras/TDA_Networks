#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import numpy as np
import scipy.stats as ss

def gaussianN_2d(mu,sigma,N):
    gpoints = np.zeros((N,2))
    for i in range(N):
        gpoints[i] = np.array([random.gauss(mu,sigma),random.gauss(mu,sigma)])
    return gpoints

def discgaussianN_2d(xmax,ymax,N,sigma):
    x = np.arange(-xmax/2,xmax/2 + 1, 1)
    xU, xL = x + 0.5, x - 0.5 
    probX = ss.norm.cdf(xU, scale = sigma) - ss.norm.cdf(xL, scale = sigma)
    probX = probX / probX.sum() # normalize the probabilities so their sum is 1
    
    
#     numsX = np.random.choice(x, size = N, p = prob)
#     numsX = numsX + round(xmax/2)
    
    y = np.arange(round(-ymax/2),round(ymax/2)+1)
    yU, yL = y + 0.5, y - 0.5 
    probY = ss.norm.cdf(yU, scale = sigma) - ss.norm.cdf(yL, scale = sigma)
    probY = probY / probY.sum() # normalize the probabilities so their sum is 1
#     numsY = np.random.choice(y, size = N, p = prob)
#     numsY = numsY + round(ymax/2)
    
    possible_points = list()
    for elx in np.arange(0,xmax+1):
        for ely in np.arange(0,ymax+1):
            possible_points.append(np.array([elx,ely]))

    Space_Position = np.array(possible_points).reshape(-1, 2)
    
    P = np.zeros(len(Space_Position))
    pos = 0
    for count, value in enumerate(y):
        for count2, value2 in enumerate(x):
            P[len(x)*count+count2] = probY[count]*probX[count2]
    
    P = P/ P.sum()
    
    return Space_Position[np.random.choice(range(len(Space_Position)), size = N, replace=False, p= P)]

#Plot configurations of points
def plot_points(points1_l,xmin,xmax,ymin,ymax,folder):
    x1_l = points1_l[:,0]
    y1_l = points1_l[:,1]
    
    plt.scatter(x1_l, y1_l, c='b')
    plt.xlim((xmin-0.5,xmax+0.5))
    plt.ylim((ymin-0.5,ymax+0.5))
    plt.xlabel('x',fontsize=12)
    plt.ylabel('y',fontsize=12)
    plt.grid()
    plt.savefig(os.path.join(folder,"PointsConfiguration.png"))
    plt.savefig(os.path.join(folder,"PointsConfiguration.svg"))
    plt.close()

