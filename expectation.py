import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button
from magnetic_resonance import P_z, plot_prob, plot_transition_rate

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


# Create main axis
fig = plt.figure(figsize=(6, 4))
ax = fig.add_subplot(111)
fig.subplots_adjust(bottom=0.2, top=0.75)
# Create axes for sliders
ax_w0 = fig.add_axes([0.3, 0.85, 0.4, 0.05])
ax_w0.spines['top'].set_visible(True)
ax_w0.spines['right'].set_visible(True)
ax_w1 = fig.add_axes([0.3, 0.92, 0.4, 0.05])
ax_w1.spines['top'].set_visible(True)
ax_w1.spines['right'].set_visible(True)
# Create sliders

w0_slider = Slider(ax=ax_w0,label='$\omega_0$',valmin=0.,valmax=2.,valinit=1.0, facecolor='#cc7000')
w1_slider = Slider(ax=ax_w1,label='$\omega_1$',valmin=0.,valmax=2.,valinit=0.1, facecolor='#cc7000')

# Plot default data

w_list = np.arange(0.6,1.4,0.01)
w0,w1 = 1.0,0.1
f_d, = ax.plot(w_list, trans_rate_exp(w0,w1,w_list),label="expectation")

 # Update values
def update(val):
    w0 = w0_slider.val
    T = w1_slider.val
    f_d.set_data(w_list, trans_rate_exp(w0,w1,w_list))
    fig.canvas.draw_idle()

w0_slider.on_changed(update)
w1_slider.on_changed(update)


plt.show()
