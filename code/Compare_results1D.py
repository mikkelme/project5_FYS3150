from dump_reader import *
from analytic_solution import *


def TwoDimPlot_compare(x,t, u_num, u_ana):
    X,T = np.meshgrid(x,t)
    colormap = "plasma"
    vmin = 0; vmax = 1
    fig = plt.figure(num=0, dpi=80, facecolor='w', edgecolor='k')
    #fig, ax = plt.subplots(ncols=1, nrows=2, dpi=80, facecolor='w', edgecolor='k')

    plt.subplot(211)
    plt.title("Numerical")
    mesh_num = plt.pcolormesh(X,T,u_num, vmin = vmin, vmax = vmax, cmap = colormap)
    #plt.xlabel("x", fontsize=14)
    plt.ylabel("time", fontsize=14)
    #plt.colorbar()

    plt.subplot(212)
    plt.title("Analytical")
    plt.pcolormesh(X,T,u_ana, vmin = vmin, vmax = vmax, cmap = colormap)
    #plt.colorbar()
    plt.xlabel("$x$", fontsize=14)
    plt.ylabel("time", fontsize=14)

    plt.tight_layout(pad=1.1, w_pad=0.7, h_pad=0.2)
    fig.subplots_adjust(right = 0.8)
    cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
    fig.colorbar(mesh_num, cax = cbar_ax, label = "$u$")
    #plt.legend(loc = "best", fontsize = 13)

    #plt.savefig("../../article/figures/phase_subplot.pdf", bbox_inches="tight")
    plt.show()


def Animation_compare(x,t, u_num, u_ana):
    import matplotlib.animation as animation
    fig = plt.figure()
    ax = plt.axes(xlim=(x[0], x[-1]), ylim = [0,1])
    numerical, = ax.plot([], [], lw=2, label = "numerical")
    analytical, = ax.plot([], [], lw=2, label = "analytical")

    # initialization function: plot the background of each frame
    def init():
        numerical.set_data([], [])
        analytical.set_data([], [])
        return numerical,

    # animation function.  This is called sequentially
    def animate(i):
        i = int(i)
        numerical.set_data(x,u_num[i])
        analytical.set_data(x,u_ana[i])
        ax.set_title(f"t = {t[i]} ({t[i]/t[-1]*100:.2f})%")
        #ax.set_ylim(0, np.max(u_ana))
        return numerical, analytical,

    # call the animator.
    time_frames = np.linspace(0,len(t)-1,len(t))
    anim = animation.FuncAnimation(fig, animate, init_func=init,
                                   frames=time_frames, interval=1, blit=False)

    plt.legend(bbox_to_anchor=(0.4, 0.80))
    plt.show()





if __name__ == "__main__":
    filename = "main.txt"
    x, t, u_num = read_dump(filename)
    u_ana = OneDim_analytic(x,t,L=1)
    TwoDimPlot_compare(x,t, u_num, u_ana)
    Animation_compare(x,t, u_num, u_ana)






#
