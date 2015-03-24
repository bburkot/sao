__author__ = 'Blazej'

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from numpy import arange, sin, pi
import threading

class Dummy():

    def plotme(self, iteration = 1):

        print "%ix plotting... " % iteration,
        t = arange(0.0, 2.0, 0.01)
        s = sin(2*pi*t)

        fig, ax = plt.subplots()
        ax.plot(t, s)
        ax.set_xlabel('time (s)')
        ax.set_ylabel('voltage (mV)')
        ax.set_title('About as simple as it gets, folks (%i)' % iteration)
        #fig.savefig("19110942_%i_test.png" % iteration)

    def threadme(self, iteration = 1):

        thread_plot = threading.Thread(target=self.plotme,
                                      args=(iteration,))
        thread_plot.start()
        thread_plot.join()

dummy = Dummy()
dummy.threadme(1)
dummy.threadme(2)
plt.show()