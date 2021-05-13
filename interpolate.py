import numpy as np
from scipy.integrate import odeint 
from scipy import interpolate
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider #for making interactive plots
# from SEQ import P_z, trans_rate_exp

plt.rc('font',family ='serif') # use Latex in matplotlib
plt.rc('text',usetex =True)

### some paramters from adapted the data on Figure 4.3, Townsend
g = 2 # g factor of electron
q = 4.803207 *10**(-10) #  esu Gaussian units
m = 1.883533* 10**(-25)  # g
c = 3.0*10**10 #cm/s
B0 = 60 #gauss
B1 = 6  #gauss
w00 = g*q/(2*m*c)*B0 
w10 = g*q/(2*m*c)*B1
T1 = 2*np.pi/w10 #12 microscond 12*10^-6, period of transition if w1 =  w10, 80kHz


## evaluate the derivatives at a specific state at time t: 
    #input: Psi, an array [aR,aI,bR,bI]; 
    #        t, the specific time; 
    #       and params, [w0,w1,w] the scale for all three frequencies

def f(Psi,t,params):
    w0,w1,w = params #w0,w1,w passed from the function input are scale numbers
    w0 = w0*w00 # w0 is a scale, times with a factor of w00 to get the actual frequency
    w1 = w1*w00
    w = w*w00
    aR = Psi[0] 
    aI = Psi[1]
    bR = Psi[2]
    bI = Psi[3]
    daR = 1/2*(w0*aI+w1*np.cos(w*t)*bI) # Real part of the derivative a
    daI = -1/2*(w0*aR+w1*np.cos(w*t)*bR) # Imaginary part of the derivative a
    dbR =  1/2*(w1*np.cos(w*t)*aI-w0*bI)
    dbI =  -1/2*(w1*np.cos(w*t)*aR-w0*bR)
    
    return np.array([daR,daI,dbR,dbI])


# find P_up, P_down and transition rate
    # input: omega_0, omega_1, omega, initial condition of the spin state r0 
    # out put: t: a array of time steps from 0 to 2T1 microseconds
    #           P_up: an array of the probability of getting +z when we measure S_z at each time step
    #           P_down: an array of the probability of getting -z when we measure S_z at each time step
    #           transition_rate: the proportion of particles that transition from +z to -z, defined as the maximum value of P_down.def P_z(w0,w1,w,Psi0,N):
def P_z(w0,w1,w,Psi0,N):
    t = np.linspace(0,2*T1,N) #time step from t=0 to t= 2T1
    params = [w0,w1,w] 
    psoln = odeint(f, Psi0, t, args=(params,)) # solve by odeint
    aR = psoln[:,0]
    aI = psoln[:,1]
    bR = psoln[:,2]
    bI = psoln[:,3]
    P_up = list(map(lambda r,i: abs(r+i*1j)**2,aR,aI))    # probablitly of measuring the spin-up: the magnitude square of the coefficent on |z+>
    P_down = list(map(lambda r,i: abs(r+i*1j)**2,bR,bI))  # 
    transition_rate = np.max(P_down)                      #defining as the maximum value of P_down as mentioned above
    return t,P_up,P_down,transition_rate

Psi0 = [1.,0.,0.,0.] # at t=0, Psi = |+z>
w0,w1 = 1.,0.42 # default value
f0,f1= np.array([w0,w1])*800 #kHZ for plot title
w_list = np.arange(0.4,1.4,0.005) # a list of applied frequency values 
rate = np.array(list(map(lambda w: P_z(w0,w1,w,Psi0,600), w_list)))[:,3] # find the transition rate at each applied freq

f = interpolate.interp1d(w_list[::10], rate[::10],kind="cubic")  #interpolate based on the approximated solutions

fig= plt.figure(figsize=(15,13))
ax = fig.add_subplot(111)
fig.subplots_adjust(bottom=0.3, top=0.8)
ax.plot(w_list,rate,label='Transition Probability from $| \langle +z | \psi \\rangle |^2$ to $| \langle -z | \psi \\rangle |^2$')
ax.plot(np.arange(0.4,1.35,0.01),f(np.arange(0.4,1.35,0.01)),label='interpolated')
ax.legend(loc='upper right', fontsize='large')
ax.set_xlabel('$\omega$', fontsize=20)
ax.xaxis.set_label_coords(1.01, -0.01)
ax.set_title('$\omega_0= %0.2f kHz$, $\omega_1 = %0.2f kHz$'%(w0,w1), fontsize=20)
ax.grid(True)
plt.savefig("C:/Users/lenovo/Desktop/PHYS302/PHYS302finalproject/PHYS302finalproject/interpolated.png")
plt.show()