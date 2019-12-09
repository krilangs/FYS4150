import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D

def get_dim(filename):
    """
    Get the dimension of the data files.
    """
    with open(filename+"_DIM.txt", "r") as infile:
        line1 = infile.readline()
        n = int(line1.split()[0])
        m = float(line1.split()[1])
        m = int(m)
    return n, m

def read_file(filename, n, m):
    """
    Read the data files and order in the the files according to the dimensions.
    """
    a = np.fromfile(filename)
    return a.reshape(m, n)

def plot_2D_EarthSun(m):
    """
    Plot 2D figure of the motion of Earth around the Sun for varying dt
    for Euler and Verlet.
    """
    filenames = ["earth_sun_euler_", "earth_sun_verlet_"]
    i = 0
    for j in range(m):
        for filename in filenames:
            exp = 10**(1+j)
            dt = 10**(-1-j)
            dt = float(dt)
            n, m = get_dim(filename + f"{exp}")
            pos_x = np.fromfile(filename + f"{exp}"+"_x")
            pos_y = np.fromfile(filename + f"{exp}"+"_y")
            plt.figure()
            if i == 0:
                plt.title(f"Forward Euler algorithm with $\\Delta t$={dt}", size=15)
            else:
                plt.title(f"Velocity Verlet algorithm with $\\Delta t$={dt}", size=15)
            plt.plot(pos_x, pos_y, label="Earth")
            plt.plot([0,0], [0,0], "o", label="Sun")
            plt.xlabel("x [AU]", size=15)
            plt.ylabel("y [AU]", size=15)
            plt.legend(loc="upper left", fontsize=15)
            plt.tight_layout()
            plt.savefig("Figures/"+filename+f"{dt}.png")
            i += 1
        plt.show()

def plot_conservation_SE():
    """
    Plot energies and angular momentum of the Earth for Euler and Verlet.
    """
    filenames = ["earth_sun_euler_1000","earth_sun_verlet_1000"]
    for filename in filenames:
        E = np.fromfile(filename + "_energy")
        K = np.fromfile(filename + "_kinetic")
        U = np.fromfile(filename + "_potential")
        L = np.fromfile(filename + "_angular")
        x = np.linspace(0, 10, len(E))
        plt.figure()
        if filename == "earth_sun_euler_1000":
            plt.title("Energies with Euler algo", size=15)
        else:
            plt.title("Energies with Verlet algo", size=15)
        plt.plot(x, E)
        plt.plot(x, K)
        plt.plot(x, U)
        plt.ylabel(r"Energies [$M_\odot$($\frac{AU}{yr})^2$]", size=15)
        plt.xlabel("Time [Years]", size=15)
        plt.legend(["Total energy", "Kinetic energy", "Potential energy"],
                   fontsize=15)
        plt.tight_layout()
        plt.savefig("Figures/Energy_"+filename+".png")

        plt.figure()
        if filename == "earth_sun_euler_1000":
            plt.title("Angular momentum with Euler algo", size=15)
        else:
            plt.title("Angular momentum with Verlet algo\n", size=15)
        plt.plot(x, L)
        plt.ylabel(r"Angular momentum [$M_E$ $\frac{AU^2}{yr}$]", size=15)
        plt.xlabel("Time [Years]", size=15)
        plt.tight_layout()
        plt.savefig("Figures/Mom_"+filename+".png")
    plt.show()

def plotBetaMinMax(filenames):
    """
    Plot orbit of Earth for velocity just above and below escape velocity.
    """
    plt.figure(figsize=(10,6))
    for i, filename in enumerate(filenames):
        plt.subplot(121+i)
        n, m = get_dim(filename)
        pos_x = np.fromfile(filename + "_x")
        pos_y = np.fromfile(filename + "_y")
        if i == 0:
            plt.title("Velocity of Earth just\n below escape velocity\n"
                      + "$v_E$ = 8.875 AU/yr", size=15)
            plt.plot(pos_x, pos_y, label="Earth")
        elif i == 1:
            plt.title("Velocity of Earth just\n above escape velocity \n"
                      + "$v_E$ = 8.89 AU/yr", size=15)
            plt.plot(pos_x, pos_y, label="Earth")
        plt.plot([0,0], [0,0], "o", label="Sun")
        plt.xlabel("x [AU]", size=15)
        plt.ylabel("y [AU]", size=15)
        plt.legend(fontsize=15)
    plt.tight_layout()
    plt.savefig("Figures/Velocities_min_max.png")
    plt.show()

def plotBeta(filenames, label):
    """
    Plot the velocity of the Earth as beta increases.
    """
    plt.figure(figsize=(10,10))
    for i, filename in enumerate(filenames):
        plt.subplot(221+i)
        n, m = get_dim(filename)
        pos_x = np.fromfile(filename + "_x")
        pos_y = np.fromfile(filename + "_y")
        if i == 0:
            plt.title("Velocity of Earth with\n v = 6.8 AU/yr, $\\beta$="+label,
                      size=15)
        elif i == 1:
            plt.title("Velocity of Earth with\n v = 7.4 AY/yr, $\\beta$="+label,
                      size=15)
        elif i == 2:
            plt.title("Velocity of Earth with\n v = 8.0 AU/yr, $\\beta=$"+label,
                      size=15)
        elif i == 3:
            plt.title("Velocity of Earth with\n v = 8.8 AU/yr, $\\beta=$"+label,
                      size=15)
        plt.plot(pos_x, pos_y, label="Earth")
        plt.plot([0,0], [0,0], "o", label="Sun")
        plt.xlabel("x [AU]", size=15)
        plt.ylabel("y [AU]", size=15)
        plt.legend(loc="lower left", fontsize=15)
    plt.tight_layout()
    plt.savefig("Figures/Beta_change_"+label+".png")
    plt.show()

def plot_EJ(filenames, names):
    """
    Plot the motion of the Earth and Jupiter with the Sun fixed as the center
    of mass of the system.
    """
    for filename in filenames:
        plt.figure()
        nEJ, mEJ = get_dim(filename)
        EJpos_x = read_file(filename + "_x", nEJ, mEJ)
        EJpos_y = read_file(filename + "_y", nEJ, mEJ)

        for j, N in enumerate(names):
            plt.plot(EJpos_x[j], EJpos_y[j], label=N)
        plt.plot([0,0], [0,0], "o", label="Sun")
        plt.xlabel("x [AU]", size=15)
        plt.ylabel("y [AU]", size=15)
        plt.title("Motion of the Earth and Jupiter", size=15)
        if filename=="EarthJupiter1000x":
            plt.xlim([-15, 10])
            plt.ylim([-10, 10])
        plt.legend(fontsize=15)
        plt.tight_layout()
        plt.savefig("Figures/" + filename + ".png")
    plt.show()

def plot_SEJ(names, filename):
    """
    Plot the motion of the three-body problem with CM of the Solar system.
    """
    n, m = get_dim(filename)
    pos_x = read_file(filename + "_x", n, m)
    pos_y = read_file(filename + "_y", n, m)

    plt.figure()
    plt.title("Motion of the Sun, Earth and Jupiter", size=15)
    for i, n in enumerate(names):
        plt.plot(pos_x[i], pos_y[i], label=n, linewidth=0.5)
    plt.xlabel("x [AU]", size=15)
    plt.ylabel("y [AU]", size=15)
    plt.legend(fontsize=15)
    plt.tight_layout()
    #plt.savefig("Figures/" + filename + ".png")

    plt.figure()
    plt.title("Motion of the Sun zoomed in", size=15)
    plt.plot(pos_x[0], pos_y[0], label="Sun", linewidth=0.5)
    plt.xlabel("x [AU]", size=15)
    plt.ylabel("y [AU]", size=15)
    plt.legend(fontsize=15)
    plt.tight_layout()
    #plt.savefig("Figures/" + filename + "_zoomed.png")
    plt.show()

def plot_cons_SEJ(filenames):
    """
    Plot energies and angular momentum for three-body problem.
    """
    for filename in filenames:
        E = np.fromfile(filename + "_energy")
        K = np.fromfile(filename + "_kinetic")
        U = np.fromfile(filename + "_potential")
        L = np.fromfile(filename + "_angular")
        x = np.linspace(0, 10, len(E))

        plt.figure()
        plt.title("Energies for three-body problem", size=15)
        plt.plot(x, E, label="Total energy")
        plt.plot(x, K, label="Kinetic energy")
        plt.plot(x, U, label="Potential energy")
        plt.ylabel(r"Energies [$M_\odot$($\frac{AU}{yr})^2$]", size=15)
        plt.xlabel("Time [Years]", size=15)
        plt.legend(fontsize=15)
        plt.tight_layout()
        #plt.savefig("Figures/Energy_"+filename+".png")

        plt.figure()
        plt.title("Angular momentum for three-body problem\n", size=15)
        plt.plot(x, L)
        plt.ylabel(r"Angular momentum [$M_E$ $\frac{AU^2}{yr}$]", size=15)
        plt.xlabel("Time [Years]", size=15)
        plt.tight_layout()
        #plt.savefig("Figures/Mom_"+filename+".png")
    plt.show()

def InnerSystem(filename):
    """
    Plot the inner Solar system in 3D and 2D.
    """
    names = ["Sun", "Mercury", "Venus", "Earth", "Mars"]
    n, m = get_dim(filename)
    pos_x = read_file(filename + "_x", n, m)
    pos_y = read_file(filename + "_y", n, m)
    pos_z = read_file(filename + "_z", n, m)

    fig = plt.figure(figsize=(8,8))
    ax = fig.gca(projection="3d")
    ax.text2D(0.3, 0.95, "3D plot of the inner Solar system",
              transform=ax.transAxes, size=15)
    years = int(2*n/100)

    ax.set_xlim(-1.5, 1.5)
    ax.set_ylim(-1.5, 1.5)
    ax.set_zlim(-1.0, 1.0)
    ax.set_xlabel("x [AU]", fontsize=15)
    ax.set_ylabel("y [AU]", fontsize=15)
    ax.set_zlabel("z [AU]", fontsize=15)
    ax.tick_params(labelsize=15)
    ax.grid(True)
    for i, N in enumerate(names):
        ax.plot(pos_x[i][0:years], pos_y[i][0:years], pos_z[i][0:years],
                label=N, linewidth=0.5)
    ax.legend(loc="upper left", fontsize=15)
    ax.xaxis.labelpad = 10
    ax.yaxis.labelpad = 10
    ax.zaxis.labelpad = 10
    fig.tight_layout()
    plt.savefig("Figures/3D_Inner_SolarSystem.png")

    plt.figure()
    plt.title("2D plot of the inner Solar system", size=15)
    for i, N in enumerate(names):
        plt.plot(pos_x[i][0:years], pos_y[i][0:years], label=N, linewidth=0.5)
    plt.legend(loc="best")
    plt.xlabel("x [AU]", size=15)
    plt.ylabel("y [AU]", size=15)
    plt.legend(fontsize=15)
    fig.tight_layout()
    plt.savefig("Figures/2D_Inner_SolarSystem.png")
    plt.show()

def SolarSystem(names, filename):
    """
    Plot the whole Solar system in 3D and 2D.
    """
    n, m = get_dim(filename)
    pos_x = read_file(filename + "_x", n, m)
    pos_y = read_file(filename + "_y", n, m)
    pos_z = read_file(filename + "_z", n, m)
    maxX = np.max(pos_x)
    maxY = np.max(pos_y)
    max = np.max([maxX, maxY])
    fig = plt.figure(figsize=(10,10))
    ax = fig.gca(projection="3d")
    ax.text2D(0.3, 0.8, "3D plot of the Solar system", transform=ax.transAxes,
              size=15)

    ax.set_xlim(-max, max)
    ax.set_ylim(-max, max)
    ax.set_zlim(-20, 20)
    ax.set_xlabel("x [AU]", fontsize=15)
    ax.set_ylabel("y [AU]", fontsize=15)
    ax.set_zlabel("z [AU]", fontsize=15)
    ax.tick_params(labelsize=15)
    ax.grid(True)
    for i, N in enumerate(names):
        ax.plot(pos_x[i], pos_y[i], pos_z[i], label=N, linewidth=0.5)
    ax.legend(loc="center left", fontsize=15)
    ax.xaxis.labelpad = 10
    ax.yaxis.labelpad = 10
    fig.tight_layout()
    plt.savefig("Figures/3D_SolarSystem.png")

    plt.figure(figsize=(7,6))
    plt.title("2D plot of the Solar system", size=15)
    for i, N in enumerate(names):
        plt.plot(pos_x[i], pos_y[i], label=N, linewidth=0.5)
    plt.legend(loc="best")
    plt.xlabel("x [AU]", size=15)
    plt.ylabel("y [AU]", size=15)
    plt.legend(fontsize=15)
    fig.tight_layout()
    plt.savefig("Figures/2D_SolarSystem.png")
    plt.show()


if __name__=="__main__":
    """Uncomment the desired task below to run:"""
    names = ["Sun", "Mercury", "Venus", "Earth", "Mars", "Jupiter", "Saturn",
             "Uranus", "Neptune", "Pluto"]

    """5c) """
    #plot_2D_EarthSun(5)
    #plot_conservation_SE()

    """5d) """
    #plotBetaMinMax(["Earth_Sun_min", "Earth_Sun_max"])
    #plotBeta(["Beta_0_1", "Beta_0_2", "Beta_0_3", "Beta_0_4"], "3.00")
    #plotBeta(["Beta_1_1", "Beta_1_2", "Beta_1_3", "Beta_1_4"], "3.33")
    #plotBeta(["Beta_2_1", "Beta_2_2", "Beta_2_3", "Beta_2_4"], "3.67")
    #plotBeta(["Beta_3_1", "Beta_3_2", "Beta_3_3", "Beta_3_4"], "4.00")

    """5e) """
    #plot_EJ(["EarthJupiter", "EarthJupiter10x"], ["Earth", "Jupiter"])
    #plot_EJ(["EarthJupiter1000x"], ["Earth", "Jupiter"])
    #plot_cons_SEJ(["EarthJupiter", "EarthJupiter10x", "EarthJupiter1000x"])

    """5f) """
    plot_SEJ(["Sun", "Earth", "Jupiter"], "SunEarthJupiter")
    plot_cons_SEJ(["SunEarthJupiter"])
    #InnerSystem("SolarSystem")
    #SolarSystem(names, "SolarSystem")




