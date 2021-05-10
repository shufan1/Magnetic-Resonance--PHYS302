# Import slider package
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from matplotlib.widgets import Slider
from magnetic_resonance import P_z, plot_prob, plot_transition_rate



r0 = [1.,0.,0.,0.]
w0,w1,w = 1.,0.1,1.
t, P_up,P_down, transition_rate = P_z(w0,w1,w,r0)
f0,f1,fap = np.array([w0,w1,w])*800 #kHZ for plot title

fig= plt.figure(figsize=(15,7))
ax = fig.add_subplot(111)
fig.subplots_adjust(bottom=0.3, top=0.8)
# plot default solution
P_zup_line, = ax.plot(t,P_down,label='Probability of spin up $| \langle -z | \psi \\rangle |^2$')
ax.set_title('$f_0= %0.2f kHz$, $f_1 = %0.2f kHz$(applied field strength), $f = %0.2f kHz$(applied field frequency)'%(f0,f1,fap), fontsize=25)
# ax.set_xlim(0,3*T1)
ax.legend(loc='upper right', fontsize='xx-large')
# ax.set_ylabel('$Probablity$ ')
ax.set_xlabel('$t$')
# ax.set_xticks(np.arange(0,3*T1,T1/2))
ax.grid(True)
transition_rate = np.max(P_down)
# ax.annotate("Transition rate = %0.2f"%transition_rate, (0.72,0.67),xycoords = "figure fraction", fontsize = 20)


# Create axes for sliders
ax_w0 = fig.add_axes([0.3, 0.18, 0.4, 0.03])
ax_w0.spines['top'].set_visible(True)
ax_w0.spines['right'].set_visible(True)
ax_w1 = fig.add_axes([0.3, 0.12, 0.4, 0.03])
ax_w1.spines['top'].set_visible(True)
ax_w1.spines['right'].set_visible(True)
ax_w = fig.add_axes([0.3, 0.05, 0.4, 0.03])
ax_w.spines['top'].set_visible(True)
ax_w.spines['right'].set_visible(True)

w0_slider = Slider(ax=ax_w0,label='$\omega_0$',valmin=0.,valmax=2.,valinit=1.0)
w1_slider = Slider(ax=ax_w1,label='$\omega_1$',valmin=0.,valmax=2.,valinit=0.1)
w_slider = Slider(ax=ax_w,label='$\omega$',valmin=0.,valmax=2.,valinit=1.0)


# Update values
def update(val):
    w0 = w0_slider.val
    w1 = w1_slider.val
    w = w_slider.val
    t,zup,zdown,transition_rate = P_z(w0,w1,w,r0) 
    P_zup_line.set_data(t, zdown)
    fig.canvas.draw_idle()

w0_slider.on_changed(update)
w1_slider.on_changed(update)
w_slider.on_changed(update)

plt.show()
