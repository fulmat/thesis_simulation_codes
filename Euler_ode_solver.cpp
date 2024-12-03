#include <iostream>
 #include <vector>
 #include <cmath>
 #include <random>
using namespace std;

    //A modell állandó paraméterei globális változóként
double kw5p = 0.25;
double kw5s = 2;
double kw6 = 1;
double kw2p = 0.2;
double kw2s = 2;
double Jw2 = 0.2;
double kw1 = 0.4;
double Jw1 = 0.2;
double kwd = 1;
double kc3p = 0.1;
double kc3s = 1;
double Jc3 = 0.05;
double Jc4 = 0.05;
double kc4 = 0.4;
double K1p = 0.1;
double K1 = 0.6;
double J1 = 0.1;
double K2p = 0.05;
double K2 = 20;
double K2s = 1;
double kwee1p = 0.08;
double kwee1s = 10;
double kcdc25p = 0.05;
double kcdc25s = 10;
double k15 = 0.25;
double k16 = 0.25;
double J15 = 0.1;
double k17p = 0.35;
double k17 = 10;
double J17 = 0.3;
double k18 = 10;
double K9 = 2.5;
double K10 = 5;
double k24 = 1000;
double k24r = 10;
double K7p = 0;
double K7 = 0.6;
double K8p = 0.1;
double K8 = 2;
double K25 = 1000;
double K25R = 10;
double J8 = 0.1;
double YE = 1;
double YB = 0.05;
double K29 = 0.05;
double K30 = 20;
double K5 = 20;
double K6p = 10;
double K6 = 100;
double HE = 0.5;
double HB = 1;
double HA = 0.5;
double RBT = 10;
double LD = 3.3;
double LE = 5;
double LB = 5;
double LA = 3;
double K20 = 10;
double K19p = 0;
double K19 = 20;
double K21 = 1;
double PP1T = 1;
double FE = 25;
double FB = 2;
double K3p = 7.5;
double K3 = 140;
double J3 = 0.01;
double J4 = 0.01;
double K4 = 40;
double GE = 0;
double GB = 1;
double GA = 0.3;
double K33 = 0.05;
double K34 = 0.05;
double K31 = 0.7;
double K32 = 1.8;
double J31 = 0.01;
double J32 = 0.01;
double K13 = 5;
double K14 = 2.5;
double J13 = 0.005;
double J14 = 0.005;
double K11p = 0;
double K11 = 1.5;
double K12 = 1.5;
double E2FT = 5;
double K22 = 1;
double K23p = 0.005;
double K23 = 1;
double K26 = 10000;
double K26R = 200;
double K27 = 0.2;
double K28 = 0.2;
double MU = 0.033;
double Kez = 0.2;
double DEPRIV = -10;
double amp = 0.005;

// Circadian rhythm modul paraméterei
double n = 2;
double J = 0.3;
double kms = 1;
double kmd = 0.1;
double kcps = 0.5;
double kcpd = 0.525;
double ka = 100;
double kd = 0.01;
double kp1 = 10;
double kcp2d = 0.0525;
double kicd = 0.01;
double kica = 20;
double kp2 = 0.1;
double Jp = 0.05;
double TFtot = 0.5;


double eps = 1;

double timer=0.0005;

struct cell { //a modell változó paraméterei
double Wee1 = 2.3;
double Wee1P = 0.077;
double Cdc25a = 0.016;
double CYCB = 0.0025;
double CycBP = 0.0028;
double ERG = 0.01218087133020163;
double DRG = 0.9005327820777893;
double CYCD = 0.012876115970612;
double CD = 0.434299469;
double CYCE = 0.01693358868360519;
double CE = 0.4828373789787;
double CYCA = 0.001462333556264639;
double CA = 0.0403826870955527;
double P27 = 0.8491207979619503;
double Cdh1 = 0.99;
double PPX = 1;
double IEP = 0.000017;
double Cdc20 = 0.00001;
double Cdc20i = 0.00009;
double E2F = 4.957176518440247;
double GM = 0.93809711933136;
//innen vettem ki az eps-t
double MASS = 1.07598838806;
double oszt = 8;
double M = 1.4;
double CP = 0.037;
double CP2 = 0.046;
double TF = 0.13;
double IC = 0.37;
};

void osc (cell &A)
{

    double lower_bound = -1;
    double upper_bound = 1;
    uniform_real_distribution<double> unif(lower_bound,upper_bound);
    default_random_engine re;
    //egyszerű kiírás, hogy lássak számokat
    
cout << A.Wee1 << "  "
         << A.Wee1P << "  "
         << A.Cdc25a << "  "
         << A.CYCB << "  "
         << A.CycBP << "  "
         << A.ERG << "  "
         << A.DRG << "  "
         << A.CYCD << "  "
         << A.CD << "  "
         << A.CYCE << "  "
         << A.CE << "  "
         << A.CYCA << "  "
         << A.CA << "  "
         << A.P27 << "  "
         << A.Cdh1 << "  "
         << A.PPX << "  "
         << A.IEP << "  "
         << A.Cdc20 << "  "
         << A.Cdc20i << "  "
         << A.E2F << "  "
         << A.GM << endl;


//cout<<A.CYCB<<endl;

// Számolt paraméterek
double CYCET = A.CYCE + A.CE;
double CYCDT = A.CYCD + A.CD;
double CYCAT = A.CYCA + A.CA;//nincs használva
double P27T = A.P27 + A.CD + A.CE + A.CA;//nincs használva
double Cdc20T= A.Cdc20 + A.Cdc20i;//nincs használva



double V4 = K4 * (GE * A.CYCE + GA * A.CYCA + GB * A.CYCB);
double V6 = K6p + K6 * (HE * A.CYCE + HA * A.CYCA + HB * A.CYCB);
double V8 = K8p + K8 * (YE * (A.CYCE + A.CYCA) + YB * A.CYCB) / (J8 + CYCET);
double PP1A = PP1T / (K21 * (FE * (A.CYCE + A.CYCA) + FB * A.CYCB) + 1);
double RBH = RBT / (K20 * (LD * CYCDT + LE * A.CYCE + LA * A.CYCA + LB * A.CYCB) / (K19p * (PP1T - PP1A) + K19 * PP1A) + 1);
double L = (K26R + K20 * (LD * CYCDT + LE * A.CYCE + LA * A.CYCA + LB * A.CYCB)) / K26;
double E2RBC = 2 * E2FT * RBH / (E2FT + RBH + L + sqrt(pow(E2FT + RBH + L, 2) - 4 * E2FT * RBH));
double E2FA = (E2FT - E2RBC) * A.E2F / E2FT;
double V2 = K2p * (1 - A.Cdh1) + K2 * A.Cdh1 + K2s * A.Cdc20;
// A változó paraméterekhez tartozó egyenletek
// Differential equations for cell cycle module
double _CYCB = eps * (K1p + K1 * pow(A.CYCB / J1, 2) / (1 + pow(A.CYCB / J1, 2))) * A.MASS - V2 * A.CYCB - (kwee1p + kwee1s * A.Wee1) * A.CYCB + (kcdc25p + kcdc25s * A.Cdc25a) * A.CycBP + unif(re) * sqrt(2 * amp * A.CYCB);
double _CycBP = (kwee1p + kwee1s * A.Wee1) * A.CYCB - ((kcdc25p +kcdc25s * A.Cdc25a)* CYCET) - V2 * CYCET +unif(re)* sqrt(2 * amp * CYCET);
double _Wee1 = (kw5p + kw5s * A.M) - kw6 * A.Wee1 - ((kw2p + kw2s * A.CYCB) * A.Wee1) / (Jw2 + A.Wee1) + (kw1 * A.Wee1P) / (Jw1 + A.Wee1P) + unif(re)* sqrt(2 * amp * A.Wee1);
double _Wee1P = ((kw2p + kw2s * A.CYCB) * A.Wee1) / (Jw2 + A.Wee1) - (kw1 * A.Wee1P) / (Jw1 + A.Wee1P) - kwd * A.Wee1P + unif(re)* sqrt(2 * amp * A.Wee1P);
double _Cdc25a = ((kc3p + kc3s * A.CYCB) * (1 - A.Cdc25a)) / (Jc3 + (1 - A.Cdc25a)) - (kc4 * A.Cdc25a) / (Jc4 + A.Cdc25a) + unif(re)* sqrt(2 * amp * A.Cdc25a);
double _ERG = eps * k15 / (1 + pow(A.DRG / J15, 2)) - k16 * A.ERG + unif(re)* sqrt(2 * amp * A.ERG);
double _DRG = eps * (k17p * A.ERG + k17 * pow(A.DRG / J17, 2) / (1 + pow(A.DRG / J17, 2))) - k18 * A.DRG + unif(re)* sqrt(2 * amp * A.DRG);
double _CYCD = eps * K9 * A.DRG + V6 * A.CD + k24r * A.CD - k24 * A.CYCD * A.P27 - K10 * A.CYCD + unif(re)* sqrt(2 * amp * A.CYCD);
double _CD = k24 * A.CYCD * A.P27 - k24r * A.CD - V6 * A.CD - K10 * A.CD + unif(re)* sqrt(2 * amp * A.CD);
double _CYCE = eps * (K7p + K7 * E2FA) - V8 * A.CYCE - K25 * A.CYCE * A.P27 + K25R * A.CE + V6 * A.CE + unif(re)*  sqrt(2 * amp * A.CYCE);
double _CE = K25 * A.CYCE * A.P27 - K25R * A.CE - V8 * A.CE - V6 * A.CE + unif(re)* sqrt(2 * amp * A.CE);
double _CYCA = eps * K29 * E2FA * A.MASS - K30 * A.Cdc20 * A.CYCA - K25 * A.CYCA * A.P27 + K25R * A.CA + V6 * A.CA + unif(re)* sqrt(2 * amp * A.CYCA);
double _CA = K25 * A.CYCA * A.P27 - K25R * A.CA - (V6 + K30 * A.Cdc20) * A.CA + unif(re)* sqrt(2 * amp * A.CA);
double _P27 = eps * K5 + K10 * A.CD + k24r * A.CD - V6 * A.P27 - K25 * A.P27 * (A.CYCE + A.CYCA) - k24 * A.CYCD * A.P27 + K25R * (A.CE + A.CA) + V8 * A.CE + K30 * A.Cdc20 * A.CA + unif(re)* sqrt(2 * amp * A.P27);
double _Cdh1 = (K3p + K3 * A.Cdc20) * (1 - A.Cdh1) / (J3 + 1 - A.Cdh1) - V4 * A.Cdh1 / (J4 + A.Cdh1) + unif(re)* sqrt(2 * amp * A.Cdh1);

double _PPX = eps * K33 - K34 * A.PPX + unif(re)* sqrt(2 * amp * A.PPX);
double _IEP = K31 * A.CYCB * (1 - A.IEP) / (J31 + 1 - A.IEP) - K32 * A.PPX * A.IEP / (J32 + A.IEP) + unif(re)* sqrt(2 * amp * A.IEP);
double _Cdc20 = K13 * A.IEP * (A.Cdc20i) / (J13 + A.Cdc20i) - K14 * A.Cdc20 / (J14 + A.Cdc20) - K12 * A.Cdc20 + unif(re)* sqrt(2 * amp * A.Cdc20);
double _Cdc20i = eps * (K11p + K11 * A.CYCB) - K13 * A.IEP * (A.Cdc20i) / (J13 + A.Cdc20i) + K14 * A.Cdc20 / (J14 + A.Cdc20) - K12 * A.Cdc20i + unif(re)* sqrt(2 * amp * A.Cdc20i);
double _MASS = eps*MU*A.GM;
// Circadian rhythm module equations
double _M = kms * pow(A.TF, n) / (pow(J, n) + pow(A.TF, n)) - kmd * A.M;
double _CP = kcps * A.M - kcpd * A.CP - 2 * ka * pow(A.CP, 2) + 2 * kd * A.CP2 - kp1 * A.CP / (Jp + A.CP + 2 * A.CP2 + 2 * A.IC);
double _CP2 = ka * pow(A.CP, 2) - kd * A.CP2 - kcp2d * A.CP2 + kicd * A.IC - kica * A.CP2 * A.TF - kp2 * A.CP2 / (Jp + A.CP + 2 * A.CP2 + 2 * A.IC);
double _TF = kcp2d * A.IC + kicd * A.IC - kica * A.TF * A.CP2 + kp2 * A.IC / (Jp + A.CP + 2 * A.CP2 + 2 * A.IC);
double _IC = TFtot - A.TF;

// adding
A.CYCB += _CYCB* timer;
A.CycBP += _CycBP* timer;
A.Wee1 += _Wee1* timer;
A.Wee1P += _Wee1P* timer;
A.Cdc25a += _Cdc25a* timer;
A.ERG += _ERG* timer;
A.DRG += _DRG* timer;
A.CYCD += _CYCD* timer;
A.CD += _CD* timer;
A.CYCE += _CYCE* timer;
A.CE += _CE* timer;
A.CYCA += _CYCA* timer;
A.CA += _CA* timer;
A.P27 += _P27* timer;
A.Cdh1 += _Cdh1* timer;
A.PPX += _PPX* timer;
A.IEP += _IEP* timer;
A.Cdc20 += _Cdc20* timer;
A.Cdc20i += _Cdc20i* timer;
A.MASS += _MASS* timer;
//circadian
A.M += _M* timer;
A.CP += _CP* timer;
A.CP2 += _CP2* timer;
A.TF += _TF* timer;
A.IC += _IC* timer;
}



int main()
{
    
    cell A;
    //egy 100 lépéses próbakör, hogy lássak számokat
    for (int i=0; i<100; i++)
    {
        osc(A);
    }
    return 0;
}