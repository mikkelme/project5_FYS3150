from dump_reader import *
from analytic_solution import *
import os

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


def error_compare(xy ,t, u, dx, dt, Time, N):
    X,Y = np.meshgrid(xy,xy)
    colormap = "plasma"
    vmin = 0; vmax = 1

    def remove_axis(i):
        if i == 0:
            ax = plt.gca()
            ax.set_xticks([])
        if i == 1:
            ax = plt.gca()
            ax.set_xticks([])
            ax.set_yticks([])
        if i == 3:
            ax = plt.gca()
            ax.set_yticks([])

    t_end = 0.1
    t_end_idx = np.argmin(np.abs(t-t_end))
    t_idx = np.linspace(0,t_end_idx, 4).astype(int)
    u_ana = TwoDim_analytic(xy, xy, t[t_idx], L=1)

    fig0 = plt.figure(num=0, dpi=80, facecolor='w', edgecolor='k')
    for i in range(4):
        plt.subplot(2,2,i+1)
        remove_axis(i)
        plt.title(f"t = {t[t_idx[i]]:.3f}")
        mesh = plt.pcolormesh(X,Y,u[t_idx[i]], vmin = vmin, vmax = vmax, cmap = colormap, linewidth=0, rasterized=True)
        mesh.set_edgecolor('face')
    plt.tight_layout(pad=1.1, w_pad=0.7, h_pad=0.2)
    fig0.subplots_adjust(right = 0.8)
    cbar_ax = fig0.add_axes([0.85, 0.15, 0.05, 0.7])
    fig0.colorbar(mesh, cax = cbar_ax, label = "u")

    fig1 = plt.figure(num=1, dpi=80, facecolor='w', edgecolor='k')
    colormap = "plasma"
    error = np.abs(u[t_idx] - u_ana)[:,1:-1,1:-1]/u_ana[:,1:-1,1:-1]
    vmin = 0; vmax = np.max(error)
    text_color = ["white", "white", "white", "black"]
    for i in range(4):
        plt.subplot(2,2,i+1)
        remove_axis(i)
        plt.title(f"t = {t[t_idx[i]]:.3f}")
        mesh = plt.pcolormesh(X[1:-1,1:-1], Y[1:-1,1:-1], error[i], vmin = vmin, vmax = vmax, cmap = colormap, linewidth=0, rasterized=True)
        mesh.set_edgecolor('face')
        plt.text(0.5,0.5,f"max = {np.max(error[i]):.2e}", fontsize = 14, ha="center", va = "center", color = text_color[i])
    plt.tight_layout(pad=1.1, w_pad=0.7, h_pad=0.2)
    plt.subplots_adjust(hspace = 0.20, wspace = 0.06)
    fig1.subplots_adjust(right = 0.8)
    cbar_ax = fig1.add_axes([0.85, 0.15, 0.05, 0.7])
    cb = fig1.colorbar(mesh, cax = cbar_ax, label = "Relativ error")
    cb.formatter.set_powerlimits((0,0))
    cb.update_ticks()
    # fig0.savefig("../article/figures/2D_sol.pdf", bbox_inches="tight")
    # fig1.savefig("../article/figures/2D_compare.pdf", bbox_inches="tight")
    plt.show()






if __name__ == "__main__":
    xy, t, u, dx, dt, Time, N = read_dump2(filename = "2DExplicit0.01.txt")
    error_compare(xy ,t, u, dx, dt, Time, N)
