import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from matplotlib.widgets import Slider
from magnetic_resonance import P_z

hbar = 4.135667696*10**-15
g = 2
q = 4.803207 *10**(-10) # in Gaussian units
m = 1.883533* 10**(-25)  
c = 3.0*10**10
B0 = 60
w00 = g*q/(2*m*c)*B0


def trans_rate_exp(w0,w1,w):
	w0 = w0*w00
	w1 = w1*w00
	w = w*w00
	y = w1**2/4/((w0-w)**2+w1**2/4)
	return y


r0 = [1.,0.,0.,0.]
w0,w1 = 1.,0.1
f0,f1= np.array([w0,w1])*800 #kHZ for plot title
w_list = np.arange(0.6,1.4,0.01)
rate = np.array(list(map(lambda w: P_z(w0,w1,w,r0), w_list)))[:,3]

fig= plt.figure(figsize=(15,13))
ax = fig.add_subplot(111)
fig.subplots_adjust(bottom=0.3, top=0.8)
transition_rate_line, = ax.plot(w_list,rate,label='Transition Probability from $| \langle +z | \psi \\rangle |^2$ to $| \langle -z | \psi \\rangle |^2$')
exp_line, = ax.plot(w_list, trans_rate_exp(w0,w1,w_list),label="expectation")

ax.legend(loc='upper right', fontsize='large')
ax.set_xlabel('$\\frac{\omega}{\omega_0}$', fontsize=20)
ax.xaxis.set_label_coords(1.01, -0.01)
ax.set_title('$f_0= %0.2f kHz$, $f_1 = %0.2f kHz$(applied field strength)'%(f0,f1), fontsize=20)
ax.grid(True)


ax_w0 = fig.add_axes([0.3, 0.18, 0.4, 0.03])
ax_w0.spines['top'].set_visible(True)
ax_w0.spines['right'].set_visible(True)
ax_w1 = fig.add_axes([0.3, 0.12, 0.4, 0.03])
ax_w1.spines['top'].set_visible(True)
ax_w1.spines['right'].set_visible(True)

w0_slider = Slider(ax=ax_w0,label='$\omega_0$',valmin=0.,valmax=2.,valinit=1.0, facecolor='#cc7000')
w1_slider = Slider(ax=ax_w1,label='$\omega_1$',valmin=0.,valmax=2.,valinit=0.1, facecolor='#cc7000')

# Update values
def update(val):
	w0 = w0_slider.val
	w1 = w1_slider.val
	rate = np.array(list(map(lambda w: P_z(w0,w1,w,r0), w_list)))[:,3]
	transition_rate_line.set_data(w_list,rate)
	exp_line.set_data(w_list, trans_rate_exp(w0,w1,w_list))
	fig.canvas.draw_idle()

w0_slider.on_changed(update)
w1_slider.on_changed(update)




plt.show()

