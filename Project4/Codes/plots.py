import matplotlib.pyplot as plt
import numpy as np

print("Choose which Task to run by typing the one of the letters (C, D, E):")
print("Task C - Equilibrium plots:")
print("Task D - probability histogram")
print("Task E & F - phase transition plots")

Task = input("Write here: ")

if Task == "C":
    MCcycles = []
    # Ordered, T=1.0
    E_mean_O1 = []; M_mean_O1 = []; Nconfigs_O1 = []
    # Unordered, T=1.0
    E_mean_U1 = []; M_mean_U1 = []; Nconfigs_U1 = []
    # Ordered, T=2.4
    E_mean_O2 = []; M_mean_O2 = []; Nconfigs_O2 = []
    # Unordered, T=2.4
    E_mean_U2 = []; M_mean_U2 = []; Nconfigs_U2 = []

    with open("Ordered_1") as file:
        lines = file.readlines()
        #Skip the first two lines
        for j in range(2,len(lines)):
            line = lines[j]
            vals = line.split()
            MCcycles.append(float(vals[0]))
            E_mean_O1.append(float(vals[1]))
            M_mean_O1.append(float(vals[2]))
            Nconfigs_O1.append(float(vals[3]))

    with open("Unordered_1") as file:
        lines = file.readlines()
        #Skip the first two lines
        for j in range(2,len(lines)):
            line = lines[j]
            vals = line.split()
            E_mean_U1.append(float(vals[1]))
            M_mean_U1.append(float(vals[2]))
            Nconfigs_U1.append(float(vals[3]))

    with open("Ordered_2_4") as file:
        lines = file.readlines()
        #Skip the first two lines
        for j in range(2,len(lines)):
            line = lines[j]
            vals = line.split()
            E_mean_O2.append(float(vals[1]))
            M_mean_O2.append(float(vals[2]))
            Nconfigs_O2.append(float(vals[3]))

    with open("Unordered_2_4") as file:
        lines = file.readlines()
        #Skip the first two lines
        for j in range(2,len(lines)):
            line = lines[j]
            vals = line.split()
            E_mean_U2.append(float(vals[1]))
            M_mean_U2.append(float(vals[2]))
            Nconfigs_U2.append(float(vals[3]))

    plt.figure()
    plt.title("Mean energy, $\\langle E \\rangle$, with ordered spin")
    plt.plot(MCcycles, E_mean_O1)
    plt.xlabel("# of Monte Carlo cycles")
    plt.ylabel("$\\langle$ E $\\rangle$")
    plt.legend(["T=1.0"])
    plt.tight_layout()
    plt.savefig("Figures/E_Ordered_1.png")

    plt.figure()
    plt.title("Mean energy, $\\langle E \\rangle$, with ordered spin")
    plt.plot(MCcycles, E_mean_O2, "orange")
    plt.xlabel("# of Monte Carlo cycles")
    plt.ylabel("$\\langle$ E $\\rangle$")
    plt.legend(["T=2.4"])
    plt.tight_layout()
    plt.savefig("Figures/E_Ordered_2_4.png")

    plt.figure()
    plt.title("Mean energy, $\\langle E \\rangle$, with random spin")
    plt.plot(MCcycles, E_mean_U1)
    plt.plot(MCcycles, E_mean_U2)
    plt.xlabel("# of Monte Carlo cycles")
    plt.ylabel("$\\langle$ E $\\rangle$")
    plt.legend(["T=1.0", "T=2.4"])
    plt.tight_layout()
    plt.savefig("Figures/E_Random.png")

    plt.figure()
    plt.title("Mean absolute Magnetization, $\\langle |M| \\rangle$,"
              + " with ordered spin")
    plt.plot(MCcycles, M_mean_O1)
    plt.xlabel("# of Monte Carlo cycles")
    plt.ylabel("$\\langle |M| \\rangle$")
    plt.legend(["T=1.0"])
    plt.tight_layout()
    plt.savefig("Figures/M_Ordered_1.png")

    plt.figure()
    plt.title("Mean absolute Magnetization, $\\langle |M| \\rangle$,"
              + " with ordered spin")
    plt.plot(MCcycles, M_mean_O2, "orange")
    plt.xlabel("# of Monte Carlo cycles")
    plt.ylabel("$\\langle |M| \\rangle$")
    plt.legend(["T=2.4"])
    plt.tight_layout()
    plt.savefig("Figures/M_Ordered_2_4.png")

    plt.figure()
    plt.title("Mean absolute Magnetization, $\\langle |M| \\rangle$,"
              + "with random spin")
    plt.plot(MCcycles, M_mean_U1)
    plt.plot(MCcycles, M_mean_U2)
    plt.xlabel("# of Monte Carlo cycles")
    plt.ylabel("$\\langle |M| \\rangle$")
    plt.legend(["T=1.0", "T=2.4"])
    plt.tight_layout()
    plt.savefig("Figures/M_Random.png")

    plt.figure()
    plt.title("# of accepted configurations (normalized) w/ ordered spin\n"
              + "as function of MC cycles")
    plt.plot(MCcycles, Nconfigs_O1)
    plt.xlabel("# of Monte Carlo cycles")
    plt.ylabel("Accepted configurations (normalized)")
    plt.legend(["T=1.0"])
    plt.tight_layout()
    plt.savefig("Figures/Accpt_Ordered_1.png")

    plt.figure()
    plt.title("# of accepted configurations (normalized) w/ ordered spin\n"
              + "as function of MC cycles")
    plt.plot(MCcycles, Nconfigs_O2, "orange")
    plt.xlabel("# of Monte Carlo cycles")
    plt.ylabel("Accepted configurations (normalized)")
    plt.legend(["T=2.4"])
    plt.tight_layout()
    plt.savefig("Figures/Accpt_Ordered_2_4.png")

    plt.figure()
    plt.title("# of accepted configurations (normalized) w/ random spin\n\
              as function of MC cycles")
    plt.plot(MCcycles, Nconfigs_U1)
    plt.plot(MCcycles, Nconfigs_U2)
    plt.xlabel("# of Monte Carlo cycles")
    plt.ylabel("Accepted configurations (normalized)")
    plt.legend(["T=1.0", "T=2.4"])
    plt.tight_layout()
    plt.savefig("Figures/Accpt_Random.png")

    with open("Nconfig_vs_Temp") as file:
        lines = file.readlines()

        T = []
        Nconfigs = []
        for i in range(2, len(lines)):
            line = lines[i]
            vals = line.split()
            T.append(float(vals[0]))
            Nconfigs.append(float(vals[1]))

        plt.figure()
        plt.title("# of accepted configurations (normalized)\n as function of T")
        plt.plot(T, Nconfigs)
        plt.xlabel("Temperatures [kT/J]")
        plt.ylabel("Accepted configurations (normalized)")
        plt.tight_layout()
        plt.savefig("Figures/AccptVsT.png")
        plt.show()


if Task == "D":
    filenames = ["Probability_1","Probability_24"]

    for i in filenames:
        with open(i) as file:
            lines = file.readlines()
            E = []
            counts = []
            max_count = 0
            most_probable_energy = 0
            for j in range(1,len(lines)):
                line = lines[j]
                vals = line.split()
                energy = float(vals[0])
                count = float(vals[1])
                E.append((energy))
                counts.append((count))
                if count > max_count:
                    max_count = count
                    most_prob_energy = energy
            plt.figure()
            props = dict(boxstyle="square", facecolor="wheat", alpha=1)
            if i == "Probability_1":
                plt.title("Probability distribution, P(E), with L=20 and T=1.0")
                plt.bar(E, counts, width = 4)
                plt.xlim(-805,-770)
                plt.text(0.3*(plt.xlim()[1]-plt.xlim()[0])+plt.xlim()[0] ,plt.ylim()[1]*0.85,
                 "Most probable energy:\n" + str(most_prob_energy), bbox = props)
                plt.xlabel("Energy")
                plt.ylabel("Energy counts, P(E)")
                plt.tight_layout()
                plt.savefig("Figures/Prob_1.png")
            else:
                plt.title("Probability distribution, P(E), with L=20 and T=2.4")
                plt.bar(E, counts, width = 3)
                plt.xlim(-705,-305)
                plt.text(0.05*(plt.xlim()[1]-plt.xlim()[0])+plt.xlim()[0] ,plt.ylim()[1]*0.85,
                 "Most probable energy:\n" + str(most_prob_energy), bbox = props)
                plt.xlabel("Energy")
                plt.ylabel("Energy counts, P(E)")
                plt.tight_layout()
                plt.savefig("Figures/Prob_2_4.png")
    plt.show()

if Task == "E":
    filename = "Phase_transitions"

    with open(filename) as file:
        lines = file.readlines()
        T = []
        E = []
        M = []
        Cv = []
        Xi = []
        split = []
        for j in range(2, len(lines)):
            line = lines[j]
            pieces = line.split()
            split.append(len(line))
            T.append(float(pieces[0]))
            E.append(float(pieces[1]))
            M.append(float(pieces[2]))
            Cv.append(float(pieces[3]))
            Xi.append(float(pieces[4]))

        t_split = len(split)/4.
        t1 = int(t_split)
        t2 = int(2*t_split)
        t3 = int(3*t_split)
        t4 = int(4*t_split)
        temps = T[:t1]

        plt.figure()
        plt.title("Phase transition plot for $\\langle E(T) \\rangle$:\n"
                  + "L=[40,60,80,100], $\\Delta$t=0.002, MC=5e5")
        plt.plot(temps, E[:t1], label="40x40")
        plt.plot(temps, E[t1:t2], label="60x60")
        plt.plot(temps, E[t2:t3], label="80x80")
        plt.plot(temps, E[t3:t4], label="100x100")
        plt.xlabel("Temperatures [kT/J]")
        plt.ylabel("$\\langle$ E(T) $\\rangle$")
        plt.tight_layout()
        plt.legend()
        plt.savefig("Figures/Phase_E(T).png")

        plt.figure()
        plt.title("Phase transition plot for $\\langle |M(T)| \\rangle$:\n"
                  + "L=[40,60,80,100], $\\Delta$t=0.002, MC=5e5")
        plt.plot(temps, M[:t1], label="40x40")
        plt.plot(temps, M[t1:t2], label="60x60")
        plt.plot(temps, M[t2:t3], label="80x80")
        plt.plot(temps, M[t3:t4], label="100x100")
        plt.xlabel("Temperatures [kT/J]")
        plt.ylabel("$\\langle |M(T)| \\rangle$")
        plt.tight_layout()
        plt.legend()
        plt.savefig("Figures/Phase_M(T).png")

        plt.figure()
        plt.title("Phase transition plot for $C_v$:\n"
                  + "L=[40,60,80,100], $\\Delta$t=0.002, MC=5e5")
        plt.plot(temps, Cv[:t1], label="40x40")
        plt.plot(temps, Cv[t1:t2], label="60x60")
        plt.plot(temps, Cv[t2:t3], label="80x80")
        plt.plot(temps, Cv[t3:t4], label="100x100")
        plt.xlabel("Temperatures [kT/J]")
        plt.ylabel("Specific heat $C_v$")
        plt.tight_layout()
        plt.legend()
        plt.savefig("Figures/Phase_Cv(T).png")

        plt.figure()
        plt.title("Phase transition plot for $\\chi$:\n"
                  + "L=[40,60,80,100], $\\Delta$t=0.002, MC=5e5")
        plt.plot(temps, Xi[:t1], label="40x40")
        plt.plot(temps, Xi[t1:t2], label="60x60")
        plt.plot(temps, Xi[t2:t3], label="80x80")
        plt.plot(temps, Xi[t3:t4], label="100x100")
        plt.xlabel("Temperatures [kT/J]")
        plt.ylabel("Susceptibility $\\chi$")
        plt.tight_layout()
        plt.legend()
        plt.savefig("Figures/Phase_Xi(T).png")
        plt.show()

if Task == "F":
    with open("") as file:
        lines = file.readlines()
    temps = []
    energylist = []
    maglist = []
    Cvlist = []
    Suscplist = []
    indeks = 0
    for i in range(1, len(lines)):
        pieces = lines[i].split()
        temps.append(float(pieces[0]))
        energylist.append(float(pieces[1]))
        maglist.append(float(pieces[2]))
        Cvlist.append(float(pieces[3]))
        Suscplist.append(float(pieces[4]))

    firstTemp = temps[0]
    for i in range(1,len(temps)):
        if temps[i] == firstTemp:
            temps = temps[0:i]
            break

    TCCv = []
    TCX = []
    for i in range(int(len(energylist)/len(temps))):
        max_temp = 0
        sublistCv = Cvlist[i*len(temps):len(temps)*(i+1)]
        sublistSuscp = Suscplist[i*len(temps):len(temps)*(i+1)]
        maxCv = max(sublistCv)
        maxSuscp = max(sublistSuscp)
        TCCv.append(temps[sublistCv.index(maxCv)])
        TCX.append(temps[sublistSuscp.index(maxSuscp)])
        print("Tc for Cv =",temps[sublistCv.index(maxCv)])
        print("Tc for X =",temps[sublistSuscp.index(maxSuscp)])

    plt.figure()
    plt.title("Mean Energy")
    plt.xlabel("T [kT/J]")
    plt.ylabel("Energy expectation value $\\langle$E$\\rangle$ [J]")
    for i in range(int(len(energylist)/len(temps))):
        plt.plot(temps,energylist[i*len(temps):len(temps)*(i+1)],"")
    plt.legend(["L = 40","L = 60","L = 80","L = 100"])

    plt.figure()
    plt.title("Absolute mean Magnetization")
    plt.xlabel("T [kT/J]")
    plt.ylabel("Magnetization expectation value $\\langle$|M|$\\rangle$ [1]")
    for i in range(int(len(energylist)/len(temps))):
        plt.plot(temps,maglist[i*len(temps):len(temps)*(i+1)],"")
    plt.legend(["L = 40","L = 60","L = 80","L = 100"])


    plt.figure()
    plt.title("Specific heat")
    plt.xlabel("T [kT/J]")
    plt.ylabel("Specific heat $\\langle$$C_v$$\\rangle$ [$J^2/kT^2$]")
    for i in range(int(len(energylist)/len(temps))):
        plt.plot(temps,Cvlist[i*len(temps):len(temps)*(i+1)],"")
    plt.legend(["L = 40","L = 60","L = 80","L = 100"])


    plt.figure()
    plt.title("Susceptibility")
    plt.xlabel("T [kT/J]")
    plt.ylabel("Susceptibility $\\langle$$\\chi$$\\rangle$ [1/kT]")
    for i in range(int(len(energylist)/len(temps))):
        plt.plot(temps,Suscplist[i*len(temps):len(temps)*(i+1)],"")
    plt.legend(["L = 40","L = 60","L = 80","L = 100"])

    plt.show()

    """
    Task f)
    """
    #Performing a linear regression to find critical temp in thermodyn. limit
    TCCv = np.array(TCCv)
    TCX = np.array(TCX)
    Llist = np.array([40,60,80,100])
    Llist = 1.0/Llist

    linreg1 = np.polyfit(Llist,TCCv,1)
    linreg2 = np.polyfit(Llist,TCX,1)

    plt.figure()
    plt.title("Specific heat $C_V$")
    plt.xlabel("$\\frac{1}{L}$")
    plt.ylabel("$T_C$ [kT/J]")
    plt.plot(Llist,TCCv,"o")
    plt.plot(Llist,np.polyval(linreg1,Llist))
    plt.legend(["$T_C$(L) from simulations","$T_C(L)$ = a$\\cdot$ $\\frac{1}{L}$ + $T_C(L = \infty)$ $\\to$ %g$\\cdot$x + %g" % (linreg1[0],linreg1[1])])

    plt.figure()
    plt.title("Susceptibility $\\chi$")
    plt.xlabel("$\\frac{1}{L}$")
    plt.ylabel("$T_C$ [kT/J]")
    plt.plot(Llist,TCX,"o")
    plt.plot(Llist,np.polyval(linreg2,Llist))
    plt.legend(["$T_C$(L) from simulations","$T_C(L)$ = a$\\cdot$ $\\frac{1}{L}$ + $T_C(L = \infty)$ $\\to$ %g$\\cdot$x + %g" % (linreg2[0],linreg2[1])])

    print("\n")
    print("The estimated Critical Temperature from our simulations is Tc = %g " % (0.5*(linreg1[1]+linreg2[1])))

    plt.show()