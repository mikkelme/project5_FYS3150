from dump_reader import *
from analytic_solution import *


def TwoDimPlot_compare(x,t, u_num, u_ana):
    X,T = np.meshgrid(x,t)
    colormap = "plasma"
    vmin = 0; vmax = 1
    fig = plt.figure(num=0, dpi=80, facecolor='w', edgecolor='k')

    plt.subplot(211)
    plt.title("Numerical")
    mesh_num = plt.pcolormesh(X,T,u_num, vmin = vmin, vmax = vmax, cmap = colormap)
    plt.ylabel("time", fontsize=14)

    plt.subplot(212)
    plt.title("Analytical")
    plt.pcolormesh(X,T,u_ana, vmin = vmin, vmax = vmax, cmap = colormap)
    plt.xlabel("$x$", fontsize=14)
    plt.ylabel("time", fontsize=14)

    plt.tight_layout(pad=1.1, w_pad=0.7, h_pad=0.2)
    fig.subplots_adjust(right = 0.8)
    cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
    fig.colorbar(mesh_num, cax = cbar_ax, label = "$u$")

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
        return numerical, analytical,

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
                                   frames=time_frames, interval=20, blit=False)

    plt.legend(bbox_to_anchor=(0.4, 0.80))
    plt.show()


def Error_compare(filenames):
    t1 = 0.1
    t2 = 0.2

    U01 = []    #dx = 0.1
    method01 = [] #order of method
    U001 = []   #dx = 0.01
    method001 = [] #order of method
    for filename in filenames:
        x, t, u, dx, dt, Time, N = read_dump(filename)
        if dx == 0.1:
            x01 = x
            t01 = t
            U01.append(u)
            method01.append(filename.split("D")[1].split("0")[0])
        elif dx == 0.01:
            x001 = x
            t001 = t
            U001.append(u)
            method001.append(filename.split("D")[1].split("0")[0])


    U01 = np.array(U01) #Shape = method, t, x
    U001 = np.array(U001) #Shape = method, t, x
    t1_idx = np.argmin(np.abs(t - t1))
    t2_idx = np.argmin(np.abs(t - t2))


    fig0 = plt.figure(num=0, dpi=80, facecolor='w', edgecolor='k')
    plt.subplots_adjust(hspace=0.35)

    plt.subplot(211)
    plt.title("dx = 0.1")
    linestyle = ["-.", "--", "--"]
    marker = "o"
    markersize = 5
    t1_idx = np.argmin(np.abs(t01 - t1))
    t2_idx = np.argmin(np.abs(t01 - t2))
    for i in range(len(U01)):
        plott1 = plt.plot(x01, U01[i,t1_idx], linestyle = linestyle[i], marker = marker, markersize = markersize, label = method01[i])
        plt.plot(x01, U01[i,t2_idx], linestyle = linestyle[i], marker = marker, markersize = markersize, color = plott1[0].get_color())
    plt.ylabel("$u$",fontsize = 14)

    plt.subplot(212)
    plt.title("dx = 0.01")
    linestyle = ["-.", "--", "--"]
    marker = None
    markersize = 5
    t1_idx = np.argmin(np.abs(t001 - t1))
    t_2idx = np.argmin(np.abs(t001 - t2))
    for i in range(len(U001)):
        plott1 = plt.plot(x001, U001[i,t1_idx], linestyle = linestyle[i], marker = marker, markersize = markersize, label = method001[i])
        plt.plot(x001, U001[i,t_2idx], linestyle = linestyle[i], marker = marker, markersize = markersize, color = plott1[0].get_color())
        plt.legend(loc = "best", fontsize = 13)
    plt.ylabel("$u$", fontsize = 14)
    plt.xlabel("x", fontsize = 14)
    plt.tight_layout(pad=1.1, w_pad=0.7, h_pad=0.2)


    fig1 = plt.figure(num=1, dpi=80, facecolor='w', edgecolor='k')
    plt.subplots_adjust(hspace=0.35)

    plt.subplot(211)
    plt.title("dx = 0.1")
    linestyle = "-"
    marker = "o"
    markersize = 5

    t1_idx = np.argmin(np.abs(t01 - t1))
    t2_idx = np.argmin(np.abs(t01 - t2))
    u_ana01 = OneDim_analytic(x01,[t01[t1_idx],t01[t2_idx]],L=1)
    for i in range(len(U01)):
        plott1 = plt.plot(x01, np.abs(U01[i,t1_idx]-u_ana01[0]), linestyle = linestyle, marker = marker, markersize = markersize, label = method01[i] + " t1")
        plt.plot(x01, np.abs(U01[i,t2_idx]-u_ana01[1])/u_ana01[1], linestyle = linestyle, marker = marker, markersize = markersize, label = method01[i] + " t2")
    plt.ylabel("Relative error",fontsize = 14)

    plt.subplot(212)
    plt.title("dx = 0.01")
    linestyle = "-"
    marker = None
    markersize = 5

    t1_idx = np.argmin(np.abs(t001 - t1))
    t2_idx = np.argmin(np.abs(t001 - t2))
    u_ana001 = OneDim_analytic(x001,[t001[t1_idx],t001[t2_idx]],L=1)
    for i in range(len(U001)):
        plott1 = plt.plot(x001, np.abs(U001[i,t1_idx]-u_ana001[0]), linestyle = linestyle, marker = marker, markersize = markersize, label = method001[i] + " t2")
        plt.plot(x001, np.abs(U001[i,t2_idx]-u_ana001[1])/u_ana001[1], linestyle = linestyle, marker = marker, markersize = markersize, label = method001[i] + " t2")
        plt.legend(loc = "best", fontsize = 13)
    plt.ylabel("Relative error",fontsize = 14)
    plt.xlabel("x",fontsize = 14)
    plt.tight_layout(pad=1.1, w_pad=0.7, h_pad=0.2)


    fig0.savefig("../article/figures/compare_show.pdf", bbox_inches="tight")
    fig1.savefig("../article/figures/compare_error.pdf", bbox_inches="tight")
    plt.show()



def Error_compare_dt(filenames):

    t1 = 0.1
    data = [[], [], []]


    method = np.array(["Explicit", "Implicit", "CN"])
    colorcycle = plt.rcParams['axes.prop_cycle'].by_key()['color']

    for filename in filenames:
        x, t, u, dx, dt, Time, N = read_dump(filename)
        t_idx = np.argmin(np.abs(t - t1))
        u_ana = OneDim_analytic(x, [t[t_idx]], L=1)

        method_idx = np.where(method == filename.split("D")[1].split("0")[0])[0][0]
        data[method_idx].append([dt, np.sum(np.abs(u - u_ana))])


    data = np.array(data)

    for i in range(len(data)):
        for j in range(len(data[i])):
            plt.plot(data[i,j,0], data[i,j,1], "o", color = colorcycle[i], label = method[i])
            print(data[i,j,0], data[i,j,1])
    plt.legend()
    plt.show()










if __name__ == "__main__":
    # filenames = ["1DExplicit0.1.txt", "1DImplicit0.1.txt", "1DCN0.1.txt", \
    #             "1DExplicit0.01.txt", "1DImplicit0.01.txt", "1DCN0.01.txt"]

    # Error_compare(filenames)

    filenames = ["1DExplicit0.1.txt", "1DImplicit0.1.txt", "1DCN0.1.txt"]
    Error_compare_dt(filenames)


    # filename = "1DExplicit0.1.txt"
    # x, t, u_num, dx, dt, Time, N = read_dump(filename)
    # u_ana = OneDim_analytic(x,t,L=1)
    # TwoDimPlot_compare(x,t, u_num, u_ana)
    # Animation_compare(x,t, u_num, u_ana)






#
