import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import matplotlib.animation as animation


class SubplotAnimation(animation.TimedAnimation):
    def __init__(self):
        fig = plt.figure()

        ax1 = fig.add_subplot(2, 1, 1)
        ax2 = fig.add_subplot(2, 1, 2)
        fig.subplots_adjust(wspace=0, hspace=0)

        gauss_data=np.loadtxt("gauss.dat")

        sim_time=gauss_data[:,0]
        cv1=gauss_data[:,1]
        cv2=gauss_data[:,2]
        height=gauss_data[:,3]

        self.t = gauss_data[:,0]
        self.cv1 = gauss_data[:,1]
        self.cv2 = gauss_data[:,2]
        self.height = gauss_data[:,3]

        ax1.set_ylabel('CVs', fontsize=14, weight='bold', style='italic')
        self.line1 = Line2D([], [], lw=0.1, color='black')       
        ax1.add_line(self.line1)
        ax1.set_xlim(min(sim_time), max(sim_time))
        ax1.set_xticks([])
        ax1.set_ylim(min(cv1)-0.02, max(cv1)+0.02)
        ax1.set_yticks([0,0.5,1.0,1.5,2])

        self.line3 = Line2D([], [], lw=0.1, color='red')       
        ax1.add_line(self.line3)

        fig.legend(('CV1', 'CV2'), bbox_to_anchor=(0.04, 0.77, 0.85, 0.1),frameon=False) 

        ax2.set_xlabel('Simulation time (ps)', fontsize=14, weight='bold', style='italic')
        ax2.set_ylabel('Gaussian Height', fontsize=14, weight='bold', style='italic')
        self.line2 = Line2D([], [], lw=0.15, color='black')        
        ax2.add_line(self.line2)
        ax2.set_xlim(min(sim_time), max(sim_time))
        ax2.set_ylim(0, max(height)+1)
        ax2.yaxis.set_ticks([0,2,4,6,8,10])

        animation.TimedAnimation.__init__(self, fig, interval=10, blit=True)

    def _draw_frame(self, framedata):
        i = framedata
        head = i - 1
        head_slice = (self.t > self.t[i] - 1.0) & (self.t < self.t[i])

        self.line1.set_data(self.t[:i], self.cv1[:i])
        
        self.line2.set_data(self.t[:i], self.height[:i])

        self.line3.set_data(self.t[:i], self.cv2[:i])

        self._drawn_artists = [self.line1, self.line2, self.line3]

    def new_frame_seq(self):
        return iter(range(self.t.size))

    def _init_draw(self):
        lines = [self.line1, self.line2, self.line3]
        for l in lines:
            l.set_data([], [])

ani = SubplotAnimation()
ani.save('mini_animation.mp4')
plt.show()
