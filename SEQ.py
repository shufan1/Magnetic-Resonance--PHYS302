
import numpy as np
from scipy.integrate import odeint 

import matplotlib.pyplot as plt
from matplotlib.widgets import Slider #for making interactive plots

plt.rc('font',family ='serif') # use Latex in matplotlib
plt.rc('text',usetex =True)


# In[24]:
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




### At resonance, and a weak oscilating field
# find the solution, and the probablity of +z and -z when w=w0
Psi0 = [1.0,0.,0.,0.] #initial condition, spin up: {{1+0j},{0+0j}}
w0,w1,w = 1.0,0.1,1.0 #
t,P_up,P_down,transition_rate = P_z(1.0,0.1,1.0,Psi0,600) #600 time step
 
fig= plt.figure(figsize=(15,13)) # initialize figure
ax = fig.add_subplot(111) # define the axies where we are going to plot P_zup and P_zdown
fig.subplots_adjust(bottom=0.3, top=0.8) # adjust the positions to make room for putting sliders later

# plot the default result, at resonance, and a weak oscilating field
P_zup_line, = ax.plot(t,P_up,label='Probability of spin up $| \langle +z | \psi \\rangle |^2$')
P_zdown_line, = ax.plot(t,P_down,label='Probability of spin down $| \langle -z | \psi \\rangle |^2$')

# set plot title, x, y axis labels and limits
f0,f1,f_ac = np.array([w0,w1, w])*800 #kHz for plot title, 
ax.set_title('$f_0= %0.2f kHz$, $f_1 = %0.2f kHz$(applied field strength), $f = %0.2f kHz$(applied field frequency)'%(f0,f1,f_ac), fontsize=25)
ax.set_xlim(0,2*T1)
ax.legend(loc='upper right', fontsize='x-large')
ax.set_xlabel('$t$ ($\mu$s)',fontsize=18)
ax.set_ylabel('$Probablity$ ', fontsize=19)
ax.xaxis.set_label_coords(1.08, -0.01)

#setting up the labels on the x axis to be the time in units of microseconds
xtic0, xtic1, xtic2, xtic3, xtic4 = 0, round(0.5*T1/(10**(-6)), 3), round(T1/(10**(-6)), 3), round(1.5*T1/(10**(-6)), 3),  round(2.0*T1/(10**(-6)), 3)
xtick_name = (xtic0, xtic1, xtic2, xtic3, xtic4) #convert the x-axis tick to units of microsecond 
ax.set_xticks(np.arange(0,2.5 *T1,T1/2))
ax.set_xticklabels(xtick_name,fontsize=18)
ax.grid(True)

# annotate to show the transition rate
annotation = ax.annotate("Transition rate = %0.2f"%transition_rate, (0.66,0.685),xycoords = "figure fraction", fontsize =20)

# Create axes for sliders
ax_w0 = fig.add_axes([0.3, 0.22, 0.4, 0.03]) #create slider axes for w0 / strength of the constant field
ax_w0.spines['top'].set_visible(True) 
ax_w0.spines['right'].set_visible(True)
ax_w1 = fig.add_axes([0.3, 0.16, 0.4, 0.03])  #create slider axes for w1 / strength of the alternating field
ax_w1.spines['top'].set_visible(True)
ax_w1.spines['right'].set_visible(True)
ax_w = fig.add_axes([0.3, 0.10, 0.4, 0.03]) #create slider axes for w / frequnecies of the oscillating field
ax_w.spines['top'].set_visible(True)
ax_w.spines['right'].set_visible(True)
ax_t = fig.add_axes([0.3, 0.05, 0.4, 0.03]) # create axes for the slider to adjust how many steps 
ax_t.spines['top'].set_visible(True)
ax_t.spines['right'].set_visible(True)
# put four sliders at where we defined our axes
w0_slider = Slider(ax=ax_w0,label='$\omega_0$',valmin=0.,valmax=2.,valinit=1.0,facecolor='#cc7000',dragging= True)
w1_slider = Slider(ax=ax_w1,label='$\omega_1$',valmin=0.,valmax=0.5,valinit=0.1,facecolor='#cc7000',dragging= True)
w_slider = Slider(ax=ax_w,label='$\omega$',valmin=0.,valmax=2.,valinit=1.0,facecolor='#cc7000',dragging= True)
t_slider = Slider(ax=ax_t,label='N steps',valmin=100.,valmax=1000.,valfmt='%0.0f',valinit=600,facecolor='#cc7000',dragging= True)


# Update values: this function is called whenever the value on the slider is changed
def update(val): 
    # retrive the values of w0, w1, w from the sliders
    w0 = w0_slider.val 
    w1 = w1_slider.val
    w = w_slider.val
    N = int(t_slider.val) # number of time step must be integer
    t,zup,zdown,transition_rate = P_z(w0,w1,w,Psi0,N) #calculate the probability of spin-up or spin-down with the new w0,w1,w values
    P_zup_line.set_data(t, zup) # update the new value on the old line 
    P_zdown_line.set_data(t, zdown)
    f0,f1,f_ac = np.array([w0,w1, w])*800 #kHz for plot title, 
    ax.set_title('$f_0= %0.2f kHz$, $f_1 = %0.2f kHz$(applied field strength), $f = %0.2f kHz$(applied field frequency)'%(f0,f1,f_ac), fontsize=25)
    annotation.set_text("Transition rate = %0.2f"%transition_rate) #update annotation
    fig.canvas.draw_idle() #pdate a figure that has been altered, but not automatically re-drawn

w0_slider.on_changed(update) #when the slider change, call the update function
w1_slider.on_changed(update)
w_slider.on_changed(update)
t_slider.on_changed(update)

plt.show()


#--------------------------------------------------------------------------------------
# ### varying the frequency of the oscillating field - transition rate

# In[72]:

### expected transition_rate vs applied AC field frequency from Townsend Figure 4.6
def trans_rate_exp(w0,w1,w):
    w0 = w0*w00
    w1 = w1*w00
    w = w*w00
    y = w1**2/4/((w0-w)**2+w1**2/4)
    return y

Psi0 = [1.,0.,0.,0.] # at t=0, Psi = |+z>
w0,w1 = 1.,0.1 # default value
f0,f1= np.array([w0,w1])*800 #kHZ for plot title
w_list = np.arange(0.4,1.4,0.005) # a list of applied frequency values 
rate = np.array(list(map(lambda w: P_z(w0,w1,w,Psi0,600), w_list)))[:,3] # find the transition rate at each applied freq

fig= plt.figure(figsize=(15,13))
ax = fig.add_subplot(111)
fig.subplots_adjust(bottom=0.3, top=0.8)
transition_rate_line, = ax.plot(w_list,rate,label='Transition Probability from $| \langle +z | \psi \\rangle |^2$ to $| \langle -z | \psi \\rangle |^2$')
exp_line, = ax.plot(w_list, trans_rate_exp(w0,w1,w_list),label="expectation")

ax.legend(loc='upper right', fontsize='large')
ax.set_xlabel('$\omega$', fontsize=20)
ax.xaxis.set_label_coords(1.01, -0.01)
ax.set_title('$f_0= %0.2f kHz$, $f_1 = %0.2f kHz$(applied field strength)'%(f0,f1), fontsize=20)
ax.grid(True)


#Create axes for sliders
ax_w0 = fig.add_axes([0.3, 0.18, 0.4, 0.03]) #create slider axes for w0 / strength of the constant field
ax_w0.spines['top'].set_visible(True) 
ax_w0.spines['right'].set_visible(True)
ax_w1 = fig.add_axes([0.3, 0.12, 0.4, 0.03])  #create slider axes for w1 / strength of the alternating field
ax_w1.spines['top'].set_visible(True)
ax_w1.spines['right'].set_visible(True)
ax_t = fig.add_axes([0.3, 0.06, 0.4, 0.03]) # create axes for the slider to adjust how many steps 
ax_t.spines['top'].set_visible(True)
ax_t.spines['right'].set_visible(True)

# put four sliders at where we defined our axes
w0_slider = Slider(ax=ax_w0,label='$\omega_0$',valmin=0.,valmax=2.,valinit=1.0,facecolor='#cc7000',dragging= True)
w1_slider = Slider(ax=ax_w1,label='$\omega_1$',valmin=0.,valmax=0.5,valinit=0.1,facecolor='#cc7000',dragging= True)
t_slider = Slider(ax=ax_t,label='N steps',valmin=100.,valmax=1000.,valfmt='%0.0f',valinit=600,facecolor='#cc7000',dragging= True)


# Update values:，this function is called whenever the value on the slider is changed

def update2(val):
    #retrive w,w1,N from the sliders
    w0 = w0_slider.val
    w1 = w1_slider.val
    N = t_slider.val
    rate = np.array(list(map(lambda w: P_z(w0,w1,w,Psi0,N), w_list)))[:,3] #calculate the transition rate under the new w and w1 values
    transition_rate_line.set_data(w_list,rate)
    exp_line.set_data(w_list, trans_rate_exp(w0,w1,w_list))

    f0,f1 = np.array([w0,w1])*800 #kHz for plot title, 
    ax.set_title('$f_0= %0.2f kHz$, $f_1 = %0.2f kHz$(applied field strength)'%(f0,f1), fontsize=25)
    # update plot title
    fig.canvas.draw_idle()

w0_slider.on_changed(update2)
w1_slider.on_changed(update2)
t_slider.on_changed(update2)

plt.show()




# # def RK_Solve(f,t0,tf,aR0,aI0,bR0,bI0,N,params):
# #     w0,w1,w = params
# #     h = (tf-t0)/N
# #     tpoints = np.arange(t0,tf,h)
# #     aR = []
# #     aI =[]
# #     bR= []
# #     bI = []

# #     r = np.array([aR0,aI0,bR0,bI0],complex)
# #     for t in tpoints:
# #         aR.append(r[0])
# #         aI.append(r[1])
# #         bR.append(r[2])
# #         bI.append(r[3])

# #         k1 = h*f(r,t,params)
# #         k2 = h*f(r+0.5*k1,t+0.5*h,params)
# #         k3 = h*f(r+0.5*k2,t+0.5*h,params)
# #         k4 = h*f(r+k3,t+h,params)
# #         r += 1/6*(k1+2.*k2+2.*k3+k4)
         
# #     return tpoints,np.array(aR),np.array(aI),np.array(bR),np.array(bI)
