#include "classes.h"

namespace StellarObjectsLibrary
{   /* Library for Stellar initial parameters from NASA.
     * Objects given as {Name of object, M, x, y, z, vx, vy, vz}.
     * M: Mass of object scaled with the Sun mass.
     * (x, y, z): Initial positions in AU from the Sun.
     * (vx, vy, vz): Initial velocities in AU/yr.
     */

    double days = 365.25;
    double SunMassFactor = 2e30;

    struct MassObject Sun = {"Sun", 1.0,
                             -3.445794798894677E-03, 7.514129097010522E-03, 1.331363550190799E-05,
                             -8.442031504341172E-06*days, -1.537062205934645E-06*days, 2.314790690956382E-07*days};

    struct MassObject Mercury = {"Mercury", 3.3e23/SunMassFactor,
                                 -6.212132461242996E-02, 3.123897172704480E-01, 3.030862955458767E-02,
                                 -3.328093197835003E-02*days, -4.292401922786234E-03*days, 2.701869400450026E-03*days};

    struct MassObject Venus = {"Venus", 4.9e24/SunMassFactor,
                               3.334775526428279E-01, -6.372726935878117E-01, -2.827692608909614E-02,
                               1.778225491438890E-02*days, 9.297275890011069E-03*days, -8.988268998481397E-04*days};

    struct MassObject Earth = {"Earth", 6e24/SunMassFactor,
                               5.335858905497771E-01, 8.370731944610594E-01, -2.812649127895525E-05,
                               -1.472817722343729E-02*days, 9.288879282185889E-03*days, -5.766212639367660E-07*days};

    struct MassObject Mars = {"Mars", 6.6e23/SunMassFactor,
                              -1.583214513581554E+00, -3.941213442522968E-01, 3.035798283835003E-02,
                              3.961264391127923E-03*days, -1.236724612249653E-02*days, -3.562885769488203E-04*days};

    struct MassObject Jupiter = {"Jupiter", 1.9e27/SunMassFactor,
                                 2.101106711987771E-01, -5.230907987162060E+00, 1.699342433382969E-02,
                                 7.447883568716252E-03*days, 6.628387409912372E-04*days, -1.693570419889271E-04*days};

    struct MassObject Saturn = {"Saturn", 5.5e26/SunMassFactor,
                                3.588562490729777E+00, -9.366213156994505E+00, 2.000047694251728E-02,
                                4.900113984196986E-03*days, 1.979143311996603E-03*days, -2.292900593455939E-04*days};

    struct MassObject Uranus = {"Uranus", 8.8e25/SunMassFactor,
                                1.631727196770728E+01, 1.125858689379779E+01, -1.695777226814316E-01,
                                -2.262627777080596E-03*days, 3.053992619552440E-03*days, 4.055554640501760E-05*days};

    struct MassObject Neptune = {"Neptune", 1.03e26/SunMassFactor,
                                 2.921158221126241E+01, -6.489302419222984E+00, -5.395752928544737E-01,
                                 6.601407945057027E-04*days, 3.082831731434076E-03*days, -7.910742860418487E-05*days};

    struct MassObject Pluto = {"Pluto", 1.31e22/SunMassFactor,
                               1.284801582259121E+01, -3.138483363860192E+01, -3.580323393621311E-01,
                               2.983920616636910E-03*days, 5.229257132807072E-04*days, -9.283492007312656E-04*days};
}
