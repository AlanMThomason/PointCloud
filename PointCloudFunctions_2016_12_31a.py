import numpy as np
#from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
#import matplotlib.widgets as widgets
import time
#import FT_20160901a as FT
import FT_20170106a as FT
import pandas as pd
#import lsf_circle as lsfc


def TimeStamp(t_then, TSMessage) : 
    t_now = time.time()
    t_span = time.time()-t_then
    t_str = '%.3f' % t_span
    if (t_then > 0) : print(TSMessage, t_str)
    return(t_now)
    
def SetOnOrigin(p,q,r):
    #print("I am in the function SetOnOrigin")
    p = p - min(p)
    q = q - min(q)
    r = r - min(r)
   
    return(p,q,r)
    
def SetOnPoint(p,q,r,index):

    p = p - p[index]
    q = q - q[index]
    r = r - r[index]
    
    return(p,q,r)
    
def FitPlane(p,q,r,npoints):
# Create the necessary matrices to find the plane representing all of the points.

    if (npoints ==0) : 
        print("No points to evaluate")
        G = [0,0,0]
        Rsqd = 0
        return(G,Rsqd)
    D = np.array([[0.0,0.0,0.0],[0.0,0.0,0.0],[0.0,0.0,0.0]],np.float64)
    E = np.array([[0.0,0.0,0.0],[0.0,0.0,0.0],[0.0,0.0,0.0]],np.float64)
    F = np.array([0.0,0.0,0.0],np.float64)

    D[0][0] = sum(p*p)
    D[1][0] = sum(p*q)
    D[2][0] = sum(p)

    D[0][1] = sum(p*q)
    D[1][1] = sum(q*q)
    D[2][1] = sum(q)

    D[0][2] = sum(p)
    D[1][2] = sum(q)
    D[2][2] = npoints

    #print("D = ",D)
    E = np.linalg.inv(D)
    F[0] = sum(p*r)
    F[1] = sum(q*r)
    F[2] = sum(r)
    #print("E = ",E)
    #print("F = ",F)
    G = np.dot(E,F)
    
    #Determine R values
    SSres = np.array(npoints, np.float64)
    res = np.array(npoints, np.float64)
    tot = np.array(npoints, np.float64)
    r_bar = np.average(r)
    tot = r - r_bar
    #print("r = ", r)
    #print("r_bar = ", r_bar)
    #print("tot = ",tot)
    tot = tot * tot
    res = (G[0]*p + G[1]*q + G[2] - r)
    #print("res = ",res)
    res = res * res
    SSres = res.sum()/npoints
    SStot = tot.sum()/npoints
    #print("SStot = ",SStot)
    #print("SSres = ",SSres)
    #print("SSres.sum() = ", SSres)
    #print("SStot.sum() = ", SStot)
    Rsqd = 1 - (SSres/SStot)
    
    #print("G = ",G)

    return(G, Rsqd)
    
def CreateSlope(p,q,npoints):
    num = sum(p*q)-sum(p)*sum(q)/npoints
    den = sum(p*p)-sum(p)*sum(p)/npoints
    slope = num/den
    intercept = (sum(q) - slope * sum(p) ) / npoints
    print("Slope  = ",slope, "Intercept = ", intercept) 
    return(slope,intercept)
    
def ReduceArraySize(p,q,r,factor):
    reduced_length = int(len(p)/factor)
    #print("Reduced Length =", reduced_length)
    p2 = np.zeros(reduced_length)
    q2 = np.zeros(reduced_length)
    r2 = np.zeros(reduced_length) 
    #print("p, q, r ",p, q, r)     
    #print("p2, q2, r2 ",p2, q2, r2) 
    j = 0
    k = 0
    # p2[0] = p[0]
    # q2[0] = q[0]
    # r2[0] = r[0] 
    for i in range(len(p)):
        #print(" i, j, k , p2, q2, r2 ",i, j, k, p2, q2, r2) 
        j= j+ 1
        if j==factor:
            j = 0
            p2[k] = p[i]
            q2[k] = q[i]
            r2[k] = r[i]
            k = k + 1
    print("reduced_length", reduced_length)
    return(p2,q2,r2, reduced_length)
            
     
def InitializeArrays(i_max):
    X = np.zeros(i_max)
    Y = np.zeros(i_max)
    Z = np.zeros(i_max)
    return(X,Y,Z)
    
def FindFileLength(FilePathAndName, limit):
    i_max= 0
    # open the input file
    infile = open(FilePathAndName, "r")
    atomsno = infile.readline().rstrip('\n').split(" ")

    while len(atomsno[0]) > 0:
        atomsno = infile.readline().rstrip('\n').split(" ")
        i_max = i_max + 1
        if i_max >= limit:
            atomsno[0]=""
            i_max = np.long(limit)
            print("Hit the Limit of ", limit, " points.")

    print(i_max)
    
    return(i_max)
    
def LoadData(FilePathAndName, p,q,r,i_max):
    infile = open(FilePathAndName, "r")
    atomsno = infile.readline().rstrip('\n').split(" ")
    i = 0
    while i < i_max:
        p[i] = np.float64(atomsno[0])
        q[i] = np.float64(atomsno[1])
        r[i] = np.float64(atomsno[2])
      
        atomsno = infile.readline().rstrip('\n').split(" ")
        i = i + 1
    infile.close()
    return(p,q,r)

    
def load_known_size(FilePathAndName, i_max) :
    x = np.empty((i_max, 3), dtype = np.float64)
    with open(FilePathAndName) as f:
        for i_max, line in enumerate(f):
            x[i_max, :] = line.split()
    return x
    
def PandaLoad(FilePathAndName) :

    df = pd.read_csv(FilePathAndName, sep=' ', header=None, dtype = np.float64)
    BigD = len(df.index)
    p = np.empty(BigD, dtype = np.float64)
    q = np.empty(BigD, dtype = np.float64)
    r = np.empty(BigD, dtype = np.float64)
    p = df[0]
    q = df[1]
    r = df[2]
    return(p, q, r)
    
def PandaLoad4(FilePathAndName) :

    df = pd.read_csv(FilePathAndName, sep=' ', header=None, dtype = np.float64)
    BigD = len(df.index)
    p = np.empty(BigD, dtype = np.float64)
    q = np.empty(BigD, dtype = np.float64)
    r = np.empty(BigD, dtype = np.float64)
    s = np.empty(BigD, dtype = np.float64)
    p = df[0]
    q = df[1]
    r = df[2]
    s = df[3]
    return(p, q, r, s)


def OutputData(FilePathAndName, p,q,r,i_max):
    outfile = open(FilePathAndName, "w")
    i = 0
    while i < i_max:
        x_str = '%.4f' % p[i]
        y_str = '%.4f' % q[i]
        z_str = '%.4f' % r[i]
        LineToWrite = x_str + " " + y_str + " " + z_str + '\n'
        outfile.write(LineToWrite)
        i = i + 1
    outfile.close()
    return(p,q,r)
    
   
def CreatePlane(p,q,r,p2,q2,r2,G,n):  
    p_max = max(p)
    q_max = max(q)
    r_max = max(r)
    p_min = min(p)
    q_min = min(q)
    r_min = min(r)


    for i in range(n):
        for j in range(n):
            p2[i*n+j] = (p_max-p_min)*i/n + p_min
            q2[i*n+j] = (q_max-q_min)*j/n + q_min
            r2[i*n+j] = (p2[i*n+j]*G[0]+q2[i*n+j]*G[1])+G[2]
    return(p2,q2,r2)
  
def RotatePoints(p,q,Angle):
    p2 = p * np.cos(Angle) + q * np.sin(Angle)
    q2 =  -1*p * np.sin(Angle) + q * np.cos(Angle) 
    print("Angle = ", Angle, " np.cos(Angle) ", np.cos(Angle), "np.sin(Angle)", np.sin(Angle))
    return(p2,q2)
   

   
def PlotData(X1,Y1, X2, Y2,X3, Y3, PlotTitle):
    fig = plt.figure(figsize = [20,20],facecolor = 'white')

    x_range = max(X1)-min(X1)
    y_range = max(Y1)-min(Y1)

    
    max_range = max(x_range, y_range)
    max_range = 2*10*int(max_range/10)
    print("max_range = ", max_range)
    x_min = int(min(X1)/10)*10
    y_min = int(min(Y1)/10)*10
    
    #range=max(max(X1), max(Y1), max(Z1))
    #range = max(int(range * 1.1),int(range)+1)
    ax = fig.add_subplot(111)
    ax.plot(X1,Y1,markersize = 1,marker = 'o', color = 'blue', linestyle = 'none')
    ax.plot(X2,Y2,markersize = 2,marker = 'o', color = 'red', linestyle = 'none')
    ax.plot(X3,Y3,markersize = 4,marker = 'o', color = 'green', linestyle = 'none')
    ax.set_xlim(x_min,x_min+max_range)
    ax.set_ylim(y_min,y_min+max_range)
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_title('XY Plane')
    ax.minorticks_on()
    ax.grid(color='g',which = 'major', linestyle = '-', linewidth =0.3)
    ax.grid(color='g',which = 'minor', linestyle = '-', linewidth =0.1)



    # Show the plot to the screen
    plt.show()



def EstimateCenterPoint(Xin, Yin, IndexArrayCP, Dens_min) :    

    # X and Y are subsets of the larger arrays.  They were derived from the larger
    #    sets using IndexArrayCP.
    
    print("Estimate Center Point")
    pi = np.pi

        
    select_len = 299  #This is the 
    
    if(Density_Calc_Outdated) :
        density_big[:] = 0.0
        LittleD = IndexArrayCP.sum()
        X = np.empty(LittleD, dtype = np.float64)
        Y = np.empty(LittleD, dtype = np.float64)
        Z = np.empty(LittleD, dtype = np.float64)
        
        X[:] = Xin
        Y[:] = Yin
        Z[:] = 0.0
        BigD = len(IndexArrayCP)
        i_orig = np.array(range(BigD))
        print("i_orig = ", i_orig)
        i_circle = np.empty(LittleD, dtype = np.int64)
        i_circle_r = np.empty(LittleD, dtype = np.float64)
        print("BigD = ", BigD, "IndexArrayCP.sum() = ", IndexArrayCP.sum(), "LittleD = ", LittleD)
        FT.crunch(i_orig, i_circle_r, IndexArrayCP, BigD, LittleD)
        i_circle = i_circle_r.astype(np.int64)
        print("i_circle = ", i_circle)
        # Now self.i_circle should represent all the indices from which data was taken.
        
        samples = min(9, LittleD)
        select_len = min(select_len, np.int(LittleD/2))
        print("select_len = ", select_len, "samples = ", samples)
        #samples is how many data points around each point we will evaluate
        #select_len is the number of datapoints we are going to take unless we reach the density limits first.
        Z = np.empty(LittleD, np.float64)
        # We aren't interested in Z, but we must create an array of zeros to use the distance subroutine
        #Z[:] =0
        X_ctrd = np.empty(LittleD, np.float64)
        Y_ctrd = np.empty(LittleD, np.float64)
        Z_ctrd = np.empty(LittleD, np.float64)
        Z_ctrd[:] = 0.0
        

        X_sel = np.empty(select_len, np.float64)
        Y_sel = np.empty(select_len, np.float64)
        
        dist = np.empty(LittleD, np.float64)
        indices = np.empty(LittleD, np.int)
        density = np.empty(LittleD, np.float64)
        #i_select = np.empty(LittleD, dtype = np.bool_)
        
    
        density[:] = 0.0

        #print("length of X, Y", len(X),len(Y))

        #tl = PCF.TimeStamp(0, "Start of loop")
        #ttl = tl
        
        
        for i in range(0, LittleD-1) :
            # For the entire subset of datapoints in the circle.
            # Center all points to the current index.
            FT.setonpoint2(X,Y,Z, X_ctrd, Y_ctrd, Z_ctrd, i, LittleD)
            #FT.setonpoint(X,Y,Z, i, LittleD )

            FT.distance2(X_ctrd,Y_ctrd,Z_ctrd, dist, LittleD)
            indices = np.argsort(dist)
            #so indices represents from the smalles to the largest the indices within the subset
            # keep in mind that if we wish to find the index of the larger subset we will have to 
            # make some manipulations.  The index of a larger subset is
            # i_circle.
            
            for j in range(1,samples) :    
                density[i] += 1/max(0.0001,dist[indices[j]]**0.5)
                print("i = ", i, "j = ", "density = ", density[i], dist[indices[j]])
            density_big[i_circle[i]] = density[i]
    
    # for i in range(0, LittleD - 1) 
        
        # if (density[i] <  Dens_min) : 
            # #print("i = ", i, "i_circle = ", i_circle, "indices = ", indices[j])
            # print("i = ", i, "density = ", density[i])
            # IndexArrayCP[i_circle[i]] = False
        
        # i should be the index within the subset
        # i_circle 
    
    #IndexArrayCP = np.where(density_big > Dens_min, True, False)
    # print("Density = ", density)
    # print("Max of Density", max(density))
    # print("Min of Density", min(density))
    # print("Dens_min = ", Dens_min)

    return(IndexArrayCP)
    # print("Density = ", density)
    # print("Max of Density", max(density))
    # i_select[:] = False
    # indices = np.argsort(density)
    # #keep in mind that 
    
    # for i in range(0,select_len) :
        # i_select[indices[LittleD-i-1]] = True
        # print("i = ", i, "indices[Bid-i-1]", indices[LittleD-i-1], X[indices[LittleD-i-1]], Y[indices[LittleD-i-1]], density[indices[LittleD-i-1]])

    # i_select[:] = False
    # i = LittleD - 1
    # keep_looking = True
    # while ((keep_looking) and (i > 0)) :
        # if (density[indices[i]] > Dens_min) :
            # i_select[indices[i]] = True
            # i = i - 1
            # print("i = ", i, X[indices[LittleD-i-1]], Y[indices[LittleD-i-1]], density[indices[LittleD-i-1]])
        # else :
            # keep_looking = False
    
    # #print("i = ", i, X[indices[LittleD-i-1]], Y[indices[LittleD-i-1]], density[indices[LittleD-i-1]])



    # FT.crunch(X, X_sel, i_select, LittleD, select_len)
    # FT.crunch(Y, Y_sel, i_select, LittleD, select_len)


    # prior to going to find circle
    # X_cent, Y_cent, R = lsfc.find_circle(X_sel, Y_sel)
    # return(X_cent, Y_cent, R, i_select)
    
def EstimateDensity(Xin, Yin, IndexArrayCP) :    

    # Xin and Yin are subsets of the larger arrays.  They were derived from the larger
    #    sets using IndexArrayCP.
    
    print("Estimate Density")

    select_len = 299  #This is the 
    
    
    LittleD = IndexArrayCP.sum()
    X = np.empty(LittleD, dtype = np.float64)
    Y = np.empty(LittleD, dtype = np.float64)
    Z = np.empty(LittleD, dtype = np.float64)
    
    X[:] = Xin
    Y[:] = Yin
    Z[:] = 0.0
    BigD = len(IndexArrayCP)
    density_big = np.empty(BigD, np.float64)
    density_big[:] = 0.0
    i_orig = np.array(range(BigD))
    #print("i_orig = ", i_orig)
    i_circle = np.empty(LittleD, dtype = np.int64)
    i_circle_r = np.empty(LittleD, dtype = np.float64)
    #print("BigD = ", BigD, "IndexArrayCP.sum() = ", IndexArrayCP.sum(), "LittleD = ", LittleD)
    FT.crunch(i_orig, i_circle_r, IndexArrayCP, BigD, LittleD)
    i_circle = i_circle_r.astype(np.int64)
    #print("i_circle = ", i_circle)
    # Now self.i_circle should represent all the indices from which data was taken.
    
    samples = min(9, LittleD)
    select_len = min(select_len, np.int(LittleD/2))
    #print("select_len = ", select_len, "samples = ", samples)
    #samples is how many data points around each point we will evaluate
    #select_len is the number of datapoints we are going to take unless we reach the density limits first.
    Z = np.empty(LittleD, np.float64)
    # We aren't interested in Z, but we must create an array of zeros to use the distance subroutine
    #Z[:] =0
    X_ctrd = np.empty(LittleD, np.float64)
    Y_ctrd = np.empty(LittleD, np.float64)
    Z_ctrd = np.empty(LittleD, np.float64)
    Z_ctrd[:] = 0.0

    dist = np.empty(LittleD, np.float64)
    indices = np.empty(LittleD, np.int)
    #i_select = np.empty(LittleD, dtype = np.bool_)
    

    #density[:] = 0.0

    #print("length of X, Y", len(X),len(Y))

    #tl = PCF.TimeStamp(0, "Start of loop")
    #ttl = tl
    
    for i in range(0, LittleD-1) :
    #for i in range(0, 2) :
        # For the entire subset of datapoints in the circle.
        # Center all points to the current index.
        FT.setonpoint2(X,Y,Z, X_ctrd, Y_ctrd, Z_ctrd, i, LittleD)
        #FT.setonpoint(X,Y,Z, i, LittleD )

        FT.distance2(X_ctrd,Y_ctrd,Z_ctrd, dist, LittleD)
        indices = np.argsort(dist)
        #so indices represents from the smalles to the largest the indices within the subset
        # keep in mind that if we wish to find the index of the larger subset we will have to 
        # make some manipulations.  The index of a larger subset is
        # i_circle.
        
        density_scalar = 0.0
        for j in range(1,samples) :    
            density_scalar += 1/max(0.0001,dist[indices[j]]**0.5)
            #print("i = ", i, " j = ", j, "indices[j] ", indices[j], "dist[indices[j]], ", dist[indices[j]], "density = ", density_scalar)
        density_big[i_circle[i]] = density_scalar 
    return(density_big)
 

def AR(InputArray, Selector, OutputArray) :
    OutputArray = InputArray[Selector]
    #return (OutputArray)
            
    

def PlotData(X_disp,Y_disp,Z_disp,x_plane,y_plane,z_plane,PlotTitle):
    fig = plt.figure(figsize = [20,5],facecolor = 'white')

    ax = fig.add_subplot(141, projection = '3d')
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_zlabel("Z")
    ax.scatter(X_disp,Y_disp,Z_disp, c='red', s = 10, marker = 'o')
    ax.scatter(x_plane, y_plane, z_plane, c = 'black', s = 1, marker = 'o')
    range=max(max(X_disp), max(Y_disp), max(Z_disp))
    range = max(int(range * 1.1),int(range)+1)
    ax.set_xlim(0,range)
    ax.set_ylim(0,range)
    ax.set_zlim(0,range)
    ax.set_title(PlotTitle)
    
    ax2 = fig.add_subplot(142)
    ax2.plot(X_disp,Y_disp,markersize = 1,marker = 'o', linestyle = 'none')
    ax2.set_xlim(0,range)
    ax2.set_ylim(0,range)
    ax2.set_xlabel("X")
    ax2.set_ylabel("Y")
    ax2.set_title('XY Plane')
    ax2.minorticks_on()
    ax2.grid(color='g',which = 'major', linestyle = '-', linewidth =0.3)
    ax2.grid(color='g',which = 'minor', linestyle = '-', linewidth =0.1)
    
    ax3 = fig.add_subplot(143)
    ax3.plot(Y_disp,Z_disp,markersize = 1,marker = 'o', linestyle = 'none')
    ax3.set_xlim(0,range)
    ax3.set_ylim(0,range)
    ax3.set_xlabel("Y")
    ax3.set_ylabel("Z")
    ax3.set_title('YZ Plane')
    ax3.minorticks_on()
    ax3.grid(color='g',which = 'major', linestyle = '-', linewidth =0.3)
    ax3.grid(color='g',which = 'minor', linestyle = '-', linewidth =0.1)
    
    ax4 = fig.add_subplot(144)
    ax4.plot(X_disp,Z_disp, markersize = 1, marker = 'o', linestyle = 'none')
    ax4.set_xlim(0,range)
    ax4.set_ylim(0,range)
    ax4.set_xlabel("X")
    ax4.set_ylabel("Z")
    ax4.set_title('XZ Plane')
    ax4.minorticks_on()
    ax4.grid(color='g',which = 'major', linestyle = '-', linewidth =0.3)
    ax4.grid(color='g',which = 'minor', linestyle = '-', linewidth =0.1)
    

    # Show the plot to the screen
    plt.show()

        
def PlotData2(X_disp,Y_disp,X2, Y2, PlotTitle):
    fig = plt.figure(figsize = [5,5],facecolor = 'white')

    #range=max(max(X_disp), max(Y_disp))
    #range = max(int(range * 1.1),int(range)+1)
    
    ax2 = fig.add_subplot(111)
    ax2.plot(X_disp,Y_disp,markersize = 2,marker = 'o', color = 'r', linestyle = 'none')
    ax2.plot(X2,Y2,markersize = 2,marker = 'o', linestyle = 'none')
    #ax2.set_xlim(0,range)
    #ax2.set_ylim(0,range)
    ax2.set_xlim(-6,6)
    ax2.set_ylim(-6,6)
    ax2.set_xlabel("X")
    ax2.set_ylabel("Y")
    ax2.set_title('XY Plane')
    ax2.minorticks_on()
    ax2.grid(color='g',which = 'major', linestyle = '-', linewidth =0.3)
    ax2.grid(color='g',which = 'minor', linestyle = '-', linewidth =0.1)
    

    

    # Show the plot to the screen
    plt.show()



