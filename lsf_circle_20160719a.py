from numpy import *
from scipy import  odr
from scipy      import optimize

def calc_R(xc, yc):
    #""" calculate the distance of each 2D points from the center (xc, yc) """
    return sqrt((x-xc)**2 + (y-yc)**2)

#@countcalls
def f_2(c):
    #""" calculate the algebraic distance between the 2D points and the mean circle centered at c=(xc, yc) """
    Ri = calc_R(*c)
    return Ri - Ri.mean()
      
def find_circle (x_in, y_in) :
 
    global x
    global y
    
    method_2  = "leastsq"
    
    x = empty(len(x_in), float64)
    y = empty(len(y_in), float64)
    
    x[:] = x_in[:]
    y[:] = y_in[:]
     
    x_m = mean(x)
    y_m = mean(y)
    
    #return

    center_estimate = x_m, y_m
    center_2, ier = optimize.leastsq(f_2, center_estimate)

    xc_2, yc_2 = center_2
    Ri_2       = calc_R(xc_2, yc_2)
    R_2        = Ri_2.mean()
    residu_2   = sum((Ri_2 - R_2)**2)
    residu2_2  = sum((Ri_2**2-R_2**2)**2)/len(x_in)
    #ncalls_2   = f_2.ncalls
    
    #print  (fmt % (method_2 , xc_2 , yc_2 , R_2 , ncalls_2 , Ri_2.std() , residu_2 , residu2_2 ))
    #print  (fmt % (method_2 , xc_2 , yc_2 , R_2 , Ri_2.std() , residu_2 , residu2_2 ))
    print(method_2)
    print("x_m", x_m, "y_m", y_m)
    print ("xc_2", xc_2 , "yc_2", yc_2)
    print ("R_2", R_2 , "Ri_2.std()", Ri_2.std() , "residu_2", residu_2 , "residu2_2", residu2_2 )
    

    return(xc_2, yc_2, R_2)
 