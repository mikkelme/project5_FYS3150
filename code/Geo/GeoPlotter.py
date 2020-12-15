import sys
sys.path.append("..")
from dump_reader import *


def TempDist(x, t, u):

    fig = plt.figure(num=0, dpi=80, facecolor='w', edgecolor='k')

    time = np.array([[10],[50],[100],[500]])
    # time = np.array([[0],[25],[150],[1000]])

    idx = np.argmin(np.abs(t-time), axis = 1)
    # print(u[idx[1]])
    for i in range(len(idx)):
        plt.plot(u[idx[i]], x/1e3, label = f"t = {t[idx[i]]:.0f} My" )

    analytical = plt.plot(x/120000*(1300 - 8) + 8, x/1e3, "--", dashes = (5,10),  label = "Analytical")
    plt.legend(loc = "best", fontsize = 13)
    plt.gca().invert_yaxis()
    plt.xlabel("T [$^{\circ}$C]", fontsize = 14)
    plt.ylabel("z [km]",fontsize = 14)
    plt.tight_layout(pad=1.1, w_pad=0.7, h_pad=0.2)
    # plt.savefig("../../article/figures/SteadyState_BRQ0.pdf", bbox_inches="tight")
    plt.show()

def max_diff(x,t,u, dx):
    # print(len(u))
    diff_temp_evo = np.zeros(len(u))
    diff_depth_evo = np.zeros(len(u))
    for i in range(len(u)):
        diff_array = np.abs(u[0] - u[i])
        diff_idx = np.argmax(diff_array)
        diff_temp = diff_array[diff_idx]
        diff_depth = diff_idx*dx

        diff_temp_evo[i] = diff_temp
        diff_depth_evo[i] = diff_depth

    peak_diff_idx = np.argmax(diff_temp_evo)

    fig0 = plt.figure(num=0, figsize=(7,8), dpi=80, facecolor='w', edgecolor='k')
    # fig0 = plt.figure(num=0, dpi=80, facecolor='w', edgecolor='k')

    plt.subplot(211)
    plt.plot(t, diff_temp_evo)
    peak_point = plt.plot(t[peak_diff_idx], diff_temp_evo[peak_diff_idx], "o", label = "Peak:    $\Delta T $  = %.2f $^{\circ}$C\n           Time = %.2f My" %(diff_temp_evo[peak_diff_idx],t[peak_diff_idx]))
    end_point = plt.plot(t[-1], diff_temp_evo[-1], "o", label ="End point:  $\Delta T$  = %.2f $^{\circ}$C \n                 Time = 1000 My " %(diff_temp_evo[-1]))
    plt.xlabel("Time [My]",fontsize = 14)
    plt.ylabel("max $\Delta$T [$^{\circ}$C]", fontsize = 14)
    plt.legend(loc = "best", fontsize = 13)
    plt.tight_layout(pad=1.1, w_pad=0.7, h_pad=0.2)


    plt.subplot(212)
    plt.plot(t, diff_depth_evo/1e3)
    plt.plot(t[peak_diff_idx], diff_depth_evo[peak_diff_idx]/1e3, "o", color = peak_point[0].get_color(), label = "Peak:    z    = %.2f km \n         Time = %.2f My" %(diff_depth_evo[peak_diff_idx]/1e3,t[peak_diff_idx]))
    plt.plot(t[-1], diff_depth_evo[-1]/1e3, "o", color = end_point[0].get_color(), label ="End point:    z    = %.2f km \n                 Time = 1000 My " %(diff_depth_evo[-1]/1e3))

    plt.xlabel("Time [My]",fontsize = 14)
    plt.ylabel("z at max $\Delta T$ [km]", fontsize = 14)
    plt.legend(loc = "best", fontsize = 13)
    plt.tight_layout(pad=1.1, w_pad=0.7, h_pad=0.2)

    fig1 = plt.figure(num=1, dpi=80, facecolor='w', edgecolor='k')
    z = np.linspace(0,dx*len(u[0])-1, len(u[0]))
    plt.plot(z/1e3, np.abs(u[0] - u[-1]), label = "Time = 1000 My")
    plt.xlabel("z [km]", fontsize = 14)
    plt.ylabel("$\Delta$T [$^{\circ}$C]",fontsize = 14)
    plt.legend(loc = "best", fontsize = 13)
    plt.tight_layout(pad=1.1, w_pad=0.7, h_pad=0.2)

    # print(np.abs(u[0] - u[-1]))
    for z,T in zip(z,np.abs(u[0] - u[-1])):
        print(f"z = {z/1e3:.1f},   T= {T:.3f}")
        if z > 20000:
            break
        # print(a/1e3, b)

    # fig0.savefig("../../article/figures/DeltaT.pdf", bbox_inches="tight")
    # fig1.savefig("../../article/figures/DeltaT_Last.pdf", bbox_inches="tight")



    plt.show()




if __name__ == "__main__":
    x, t, u, dx, dt, Time, N = read_dump(filename = "1DTempDist_AR_decay1.txt")
    max_diff(x, t, u, dx)
    #TempDist(x, t, u)
    # Animation_show(x, t, u, 1)
