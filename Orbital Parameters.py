
# coding: utf-8

# In[36]:


import numpy as np


# In[37]:


# This function will calculate the length of the semi-major axis of an ellipse given the distance at perigee and apogee.

def semi_major(r_p,r_a):
    a = (r_p + r_a)/2
    return a


# In[38]:


# This function will calculate the period (in seconds) of an ellipse given the length of the semi-major axis. Mu is the 
# gravitational parameter and is defined by the gravitational constant (G) times the sum of the two masses in the 
# system (G(M_1 + M_2)). However, if M_1 >> M_2, the expression can omit the value of M_2 (GM_1). The default value for 
# mu has been set to Earth's value (398,600 km^3/s^2).

def orbital_period(a, mu=None):
    if mu == None:
        mu = 398600
    period = ((2*np.pi)/((mu)**(1/2)))*(a)**(3/2)
    T = period
    return T


# In[39]:


# This function will calculate an object's velocity using the energy equation. As arguments, it requires the shape (first
# letter) of the orbit and the object's position. As optional arguments, it accepts a value for the semi-major axis
# and the gravitational parameter, mu. If mu is omitted, the default value for Earth will be used. If "a" is omitted, 
# and the shape is anything other than a circle or parabola, you will receive a "TypeError".

def velocity_E(shape, r, a=None, mu=None):
    if mu == None:
        mu = 398600
    if shape == "c":
        v = (mu/r)**(1/2)
    elif shape == "e":
        v = (mu*(2/r - 1/a))**(1/2)
    elif shape == "p":
        v = (2*mu/r)**(1/2)
    else:
        v = (mu*(1/a + 2/r))**(1/2)
    return v


# In[123]:


# This function will calculate on object's angular momentum in vector form, given vectors for position and velocity. 

def ang_mom(R, V):
    H = np.cross(R,V)
    h = (np.dot(H,H))**(1/2)
    return H, h


# In[99]:


# This function uses the distance fomula to calculate the craft's distance from the focus.

def distance(h, e, theta, mu = None):
    if mu == None:
        mu = 398600
    r = (h**2/mu)*(1/(1 + e*np.cos(theta)))
    return r


# In[100]:


# This function calculates the right ascension of the ascending node.

def cap_omega(N,n):
    if N[1] >= 0:
        o = np.arccos(N[0]/n)
    else:
        o = 2*np.pi - np.arccos(N[0]/n)
        
    O = o*180/(np.pi)
    
    return O


# In[101]:


# This function calculates the argument of perigee for an orbit.

def arg_perigee(N, e, n, e_abs):
    if e[2] >= 0:
        W = np.arccos(np.dot(N,e)/(n*e_abs))
    else:
        W = 2*np.pi - np.arccos(np.dot(N,e)/(n*e_abs))
    w = W*180/np.pi
    
    return w


# In[132]:


# This function will return all orbital elements given the position and velocity vectors of a craft.

def orbital_elements(R,V,mu=None):
    
    if mu == None:
        mu = 398600
    
    #Absolute value of R
    r = (np.dot(R,R))**(1/2)
    
    #Absolute value of V
    v = (np.dot(V,V))**(1/2)
    
    #Radial Component of Velocity
    v_rad = np.dot(R,V)/r
    
    #Tangential Componenet of Velocity
    v_tan = (v**2 - v_rad**2)**(1/2)
    
    #calling ang_mom to get vector and scalar forms of h
    H, h = ang_mom(R,V)
    
    #using H_z (H[2]) and h to find the inclination (i) of the orbit
    i_rad = np.arccos(H[2]/h)
    i = i_rad*180/np.pi
    
    #using the unit vector k and the vector form of h to calculate the node line (N) and its magnitude
    k = [0,0,1]
    N = np.cross(k,H)
    n = (np.dot(N,N)**(1/2))
    
    #using N_x (N[0]) and the scalar form of N to calculate the right ascension of the ascending node.
    O = cap_omega(N,n)
    
    #Eccentricity vector
    e = (1/mu)*((v**2 - mu/r)*R - r*v_rad*V)
    e_abs = (np.dot(e,e))**(1/2)
    
    #Argument of perigee
    w = arg_perigee(N, e, n, e_abs)
        
    #True Anolmaly
    if v_rad >= 0:
        Theta = np.arccos(np.dot(e,R)/(e_abs*r))
    else:
        Theta = 2*np.pi - np.arccos(np.dot(e,R)/(e_abs*r))
    theta = Theta*180/np.pi
    
    r_p = distance(h, e_abs, 0)
    r_a = distance(h, e_abs, np.pi)
    a = semi_major(r_p,r_a)
    b = a*(1-e**2)**(1/2)
    
    T = orbital_period(a)
    
    
    print("The distance from the focus to the craft is: " + str(r))
    print("The magnitude of the craft's current velocity is: " + str(v))
    print("The magnitude of the craft's current radial velocity is: " + str(v_rad))
    print("The magnitude of the craft's current tangential velocity is: " + str(v_tan))
    print("The angular momentum vector is: " + str(H))
    print("The magnitude of the angular momentum is: " + str(h))
    print("The inclination of the orbit is: " + str(i))
    print("The vector for the ascending node is: " + str(N))
    print("The magnitude of the ascending node is: " + str(n))
    print("The right ascension of the ascending node is: " + str(O))
    print("The eccentricity vector is: " + str(e))
    print("The eccentricity of the orbit is: " + str(round(e_abs,4)))
    print("The argument of perigee for the orbit is: " + str(round(w, 2)) + " degrees.")
    print("The true anomaly for the craft's position is: " + str(round(theta, 2)) + " degrees.")
    print("The distance at perigee is: " + str(r_p))
    print("The distance at apogee is: " + str(r_a))
    print("The semi-major axis of the orbit is: " + str(a))
    print("The semi-minor axis of the orbit is: " + str(b))
    print("The period of the orbit (in seconds) is: " + str(T))
        
    return r, v, v_rad, v_tan, H, h, i, N, n, O, e, e_abs, w, theta, r_p, r_a, a, b, T


    


# In[137]:


R = np.asarray([-6600, -1300, -5200], dtype='int64')
V = np.asarray([-0.4, -0.5, -0.6])


# In[138]:


r, v, v_rad, v_tan, H, h, i, N, n, O, e, e_abs, w, theta, r_p, r_a, a, b, T = orbital_elements(R,V)

