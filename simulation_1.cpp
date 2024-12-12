#include <iostream>
#include <vector>
#include <cstring>
#include <cmath>
#include <random>
#include <iterator>
//#include <bits/stdc++.h>

#include <afxwin.h> 
#include <afxext.h> 

#include <windows.h>
#include <winuser.h>
//#include "include/CCanvasWnd.h"

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
// Circadian rhythm modul
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
// Simulation dependent constants
double eps = 1;
double timer = 0.0005;
double WR = 17;
int id = 1;
double maxWNT = 0.00032;
double pWNTDG = 0.2; //pWNT degradation factor
int gridsize = 60;
double dfrate = 0.4;
int numofit = 2400;
double pcdcs = 0.002411; //paneth cell death chance slope 
double egcdcs = 0.0868; //enterocyte and goblet cell death chance slope 
double maxdistance = sqrt(2 * pow(gridsize / 2, 2));
double WNTPC = maxWNT / sqrt(2 * pow(gridsize / 4, 2));
double BMPPC = 100 / maxdistance;
double spddd = 55; //secretory progenitor probability (division daughter determinator)

struct intra {
    string type = "";
    int ID = 0;
    int pID = 0;
    float age = 0; //age being float is important in the data analysis function
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


struct extra {
    double gradWNT;
    double pWNT = 0;
    double WNT;
    double BMP;
};

namespace Simulation {
    static void move(int a, int b, intra temp, vector<vector<intra>>& v);
    static void pmove(vector<vector<intra>>v);
    static void show(vector<vector<vector<intra>>> full);
}

void inload(vector<vector<intra>>& v)
{
    //creating the grid with empty structs
    vector<intra> vv;
    intra z;
    memset(&z, 0, sizeof(z));

    for (int i = 0; i < gridsize; i++)
    {
        vv.push_back(z);
    }
    for (int j = 0; j < gridsize; j++)
    {
        v.push_back(vv);
    }

    //stem cell colony with initial values
    intra s;
    s.type = "s";
    for (int i = 0; i < v.size(); i++)
    {
        for (int j = 0; j < v[i].size(); j++)
        {
            if (pow(i - gridsize / 2, 2) + pow(j - gridsize / 2, 2) <= 4)
            {
                v[i][j] = s;

                v[i][j].ID = id;

                id++;
            }
        }
    }
}


void exload(vector<vector<extra>>& v)
{
    vector<extra>vv;
    extra z;

    for (int i = 0; i < gridsize; i++)
    {
        vv.push_back(z);
    }
    for (int j = 0; j < gridsize; j++)
    {
        v.push_back(vv);
    }


    for (int i = 0; i < v.size(); i++)
    {
        for (int j = 0; j < v[i].size(); j++)
        {
            double distance = sqrt(pow(i - gridsize / 2, 2) + pow(j - gridsize / 2, 2));
            if (distance < maxdistance / 2) v[i][j].WNT = maxWNT - distance * WNTPC;
            v[i][j].BMP = BMPPC * distance;
            v[i][j].gradWNT = v[i][j].WNT;
        }
    }
}

void aging(vector<vector<intra>>& v)
{
    for (int i = 0; i < v.size(); i++)
    {
        for (int j = 0; j < v[i].size(); j++)
        {
            if (v[i][j].ID != 0)   v[i][j].age += 1;
        }
    }
}

void intraODEs(vector<vector<intra>>& v, vector<vector<extra>>w)
{
    double lower_bound = -1;
    double upper_bound = 1;
    uniform_real_distribution<double> unif(lower_bound, upper_bound);
    default_random_engine re;

    for (int i = 0; i < v.size(); i++)
    {
        for (int j = 0; j < v[i].size(); j++)
        {
            if (v[i][j].ID != 0 && (v[i][j].type == "s" || v[i][j].type == "t" || v[i][j].type == "ap" || v[i][j].type == "sp"))
            {
                int step = 0;
                while (step != 4)
                {
                    double CYCET = v[i][j].CYCE + v[i][j].CE;
                    double CYCDT = v[i][j].CYCD + v[i][j].CD;
                    double CYCAT = v[i][j].CYCA + v[i][j].CA;//nincs használva
                    double P27T = v[i][j].P27 + v[i][j].CD + v[i][j].CE + v[i][j].CA;//nincs használva
                    double Cdc20T = v[i][j].Cdc20 + v[i][j].Cdc20i;//nincs használva
                    double V4 = K4 * (GE * v[i][j].CYCE + GA * v[i][j].CYCA + GB * v[i][j].CYCB);
                    double V6 = K6p + K6 * (HE * v[i][j].CYCE + HA * v[i][j].CYCA + HB * v[i][j].CYCB);
                    double V8 = K8p + K8 * (YE * (v[i][j].CYCE + v[i][j].CYCA) + YB * v[i][j].CYCB) / (J8 + CYCET);
                    double PP1A = PP1T / (K21 * (FE * (v[i][j].CYCE + v[i][j].CYCA) + FB * v[i][j].CYCB) + 1);
                    double RBH = RBT / (K20 * (LD * CYCDT + LE * v[i][j].CYCE + LA * v[i][j].CYCA + LB * v[i][j].CYCB) / (K19p * (PP1T - PP1A) + K19 * PP1A) + 1);
                    double L = (K26R + K20 * (LD * CYCDT + LE * v[i][j].CYCE + LA * v[i][j].CYCA + LB * v[i][j].CYCB)) / K26;
                    double E2RBC = 2 * E2FT * RBH / (E2FT + RBH + L + sqrt(pow(E2FT + RBH + L, 2) - 4 * E2FT * RBH));
                    double E2FA = (E2FT - E2RBC) * v[i][j].E2F / E2FT;
                    double V2 = K2p * (1 - v[i][j].Cdh1) + K2 * v[i][j].Cdh1 + K2s * v[i][j].Cdc20;
                    // Differential equations for cell cycle module
                    double _CYCB = eps * (K1p + K1 * pow(v[i][j].CYCB / J1, 2) / (1 + pow(v[i][j].CYCB / J1, 2))) * v[i][j].MASS - V2 * v[i][j].CYCB - (kwee1p + kwee1s * v[i][j].Wee1) * v[i][j].CYCB + (kcdc25p + kcdc25s * v[i][j].Cdc25a) * v[i][j].CycBP - w[i][j].WNT * WR + unif(re) * sqrt(2 * amp * v[i][j].CYCB);
                    double _CycBP = (kwee1p + kwee1s * v[i][j].Wee1) * v[i][j].CYCB - ((kcdc25p + kcdc25s * v[i][j].Cdc25a) * CYCET) - V2 * CYCET + unif(re) * sqrt(2 * amp * CYCET);
                    double _Wee1 = (kw5p + kw5s * v[i][j].M) - kw6 * v[i][j].Wee1 - ((kw2p + kw2s * v[i][j].CYCB) * v[i][j].Wee1) / (Jw2 + v[i][j].Wee1) + (kw1 * v[i][j].Wee1P) / (Jw1 + v[i][j].Wee1P) + unif(re) * sqrt(2 * amp * v[i][j].Wee1);
                    double _Wee1P = ((kw2p + kw2s * v[i][j].CYCB) * v[i][j].Wee1) / (Jw2 + v[i][j].Wee1) - (kw1 * v[i][j].Wee1P) / (Jw1 + v[i][j].Wee1P) - kwd * v[i][j].Wee1P + unif(re) * sqrt(2 * amp * v[i][j].Wee1P);
                    double _Cdc25a = ((kc3p + kc3s * v[i][j].CYCB) * (1 - v[i][j].Cdc25a)) / (Jc3 + (1 - v[i][j].Cdc25a)) - (kc4 * v[i][j].Cdc25a) / (Jc4 + v[i][j].Cdc25a) + unif(re) * sqrt(2 * amp * v[i][j].Cdc25a);
                    double _ERG = eps * k15 / (1 + pow(v[i][j].DRG / J15, 2)) - k16 * v[i][j].ERG + unif(re) * sqrt(2 * amp * v[i][j].ERG);
                    double _DRG = eps * (k17p * v[i][j].ERG + k17 * pow(v[i][j].DRG / J17, 2) / (1 + pow(v[i][j].DRG / J17, 2))) - k18 * v[i][j].DRG + unif(re) * sqrt(2 * amp * v[i][j].DRG);
                    double _CYCD = eps * K9 * v[i][j].DRG + V6 * v[i][j].CD + k24r * v[i][j].CD - k24 * v[i][j].CYCD * v[i][j].P27 - K10 * v[i][j].CYCD + unif(re) * sqrt(2 * amp * v[i][j].CYCD);
                    double _CD = k24 * v[i][j].CYCD * v[i][j].P27 - k24r * v[i][j].CD - V6 * v[i][j].CD - K10 * v[i][j].CD + unif(re) * sqrt(2 * amp * v[i][j].CD);
                    double _CYCE = eps * (K7p + K7 * E2FA) - V8 * v[i][j].CYCE - K25 * v[i][j].CYCE * v[i][j].P27 + K25R * v[i][j].CE + V6 * v[i][j].CE + unif(re) * sqrt(2 * amp * v[i][j].CYCE);
                    double _CE = K25 * v[i][j].CYCE * v[i][j].P27 - K25R * v[i][j].CE - V8 * v[i][j].CE - V6 * v[i][j].CE + unif(re) * sqrt(2 * amp * v[i][j].CE);
                    double _CYCA = eps * K29 * E2FA * v[i][j].MASS - K30 * v[i][j].Cdc20 * v[i][j].CYCA - K25 * v[i][j].CYCA * v[i][j].P27 + K25R * v[i][j].CA + V6 * v[i][j].CA + unif(re) * sqrt(2 * amp * v[i][j].CYCA);
                    double _CA = K25 * v[i][j].CYCA * v[i][j].P27 - K25R * v[i][j].CA - (V6 + K30 * v[i][j].Cdc20) * v[i][j].CA + unif(re) * sqrt(2 * amp * v[i][j].CA);
                    double _P27 = eps * K5 + K10 * v[i][j].CD + k24r * v[i][j].CD - V6 * v[i][j].P27 - K25 * v[i][j].P27 * (v[i][j].CYCE + v[i][j].CYCA) - k24 * v[i][j].CYCD * v[i][j].P27 + K25R * (v[i][j].CE + v[i][j].CA) + V8 * v[i][j].CE + K30 * v[i][j].Cdc20 * v[i][j].CA + unif(re) * sqrt(2 * amp * v[i][j].P27);
                    double _Cdh1 = (K3p + K3 * v[i][j].Cdc20) * (1 - v[i][j].Cdh1) / (J3 + 1 - v[i][j].Cdh1) - V4 * v[i][j].Cdh1 / (J4 + v[i][j].Cdh1) + unif(re) * sqrt(2 * amp * v[i][j].Cdh1);
                    double _PPX = eps * K33 - K34 * v[i][j].PPX + unif(re) * sqrt(2 * amp * v[i][j].PPX);
                    double _IEP = K31 * v[i][j].CYCB * (1 - v[i][j].IEP) / (J31 + 1 - v[i][j].IEP) - K32 * v[i][j].PPX * v[i][j].IEP / (J32 + v[i][j].IEP) + unif(re) * sqrt(2 * amp * v[i][j].IEP);
                    double _Cdc20 = K13 * v[i][j].IEP * (v[i][j].Cdc20i) / (J13 + v[i][j].Cdc20i) - K14 * v[i][j].Cdc20 / (J14 + v[i][j].Cdc20) - K12 * v[i][j].Cdc20 + unif(re) * sqrt(2 * amp * v[i][j].Cdc20);
                    double _Cdc20i = eps * (K11p + K11 * v[i][j].CYCB) - K13 * v[i][j].IEP * (v[i][j].Cdc20i) / (J13 + v[i][j].Cdc20i) + K14 * v[i][j].Cdc20 / (J14 + v[i][j].Cdc20) - K12 * v[i][j].Cdc20i + unif(re) * sqrt(2 * amp * v[i][j].Cdc20i);
                    double _MASS = eps * MU * v[i][j].GM;

                    v[i][j].CYCB += _CYCB * timer;
                    v[i][j].CycBP += _CycBP * timer;
                    v[i][j].Wee1 += _Wee1 * timer;
                    v[i][j].Wee1P += _Wee1P * timer;
                    v[i][j].Cdc25a += _Cdc25a * timer;
                    v[i][j].ERG += _ERG * timer;
                    v[i][j].DRG += _DRG * timer;
                    v[i][j].CYCD += _CYCD * timer;
                    v[i][j].CD += _CD * timer;
                    v[i][j].CYCE += _CYCE * timer;
                    v[i][j].CE += _CE * timer;
                    v[i][j].CYCA += _CYCA * timer;
                    v[i][j].CA += _CA * timer;
                    v[i][j].P27 += _P27 * timer;
                    v[i][j].Cdh1 += _Cdh1 * timer;
                    v[i][j].PPX += _PPX * timer;
                    v[i][j].IEP += _IEP * timer;
                    v[i][j].Cdc20 += _Cdc20 * timer;
                    v[i][j].Cdc20i += _Cdc20i * timer;
                    v[i][j].MASS += _MASS * timer;

                    if (v[i][j].type == "s") step += 2;
                    if (v[i][j].type == "t" || v[i][j].type == "ap") step++;
                    if (v[i][j].type == "sp") step = 4;
                }
                //circadian
                double _M = kms * pow(v[i][j].TF, n) / (pow(J, n) + pow(v[i][j].TF, n)) - kmd * v[i][j].M;
                double _CP = kcps * v[i][j].M - kcpd * v[i][j].CP - 2 * ka * pow(v[i][j].CP, 2) + 2 * kd * v[i][j].CP2 - kp1 * v[i][j].CP / (Jp + v[i][j].CP + 2 * v[i][j].CP2 + 2 * v[i][j].IC);
                double _CP2 = ka * pow(v[i][j].CP, 2) - kd * v[i][j].CP2 - kcp2d * v[i][j].CP2 + kicd * v[i][j].IC - kica * v[i][j].CP2 * v[i][j].TF - kp2 * v[i][j].CP2 / (Jp + v[i][j].CP + 2 * v[i][j].CP2 + 2 * v[i][j].IC);
                double _TF = kcp2d * v[i][j].IC + kicd * v[i][j].IC - kica * v[i][j].TF * v[i][j].CP2 + kp2 * v[i][j].IC / (Jp + v[i][j].CP + 2 * v[i][j].CP2 + 2 * v[i][j].IC);
                double _IC = TFtot - v[i][j].TF;
                v[i][j].M += _M * timer;
                v[i][j].CP += _CP * timer;
                v[i][j].CP2 += _CP2 * timer;
                v[i][j].TF += _TF * timer;
                v[i][j].IC += _IC * timer;
            }
        }
    }
}

void exocytosis(vector<vector<intra>>& v, vector<vector<extra>>& w)
{
    for (int i = 0; i < gridsize; i++)
    {
        for (int j = 0; j < gridsize; j++)
        {
            if (v[i][j].type == "p" && pow(i - gridsize / 2, 2) + pow(j - gridsize / 2, 2) <= 25)
            {
                w[i][j].pWNT += 0.000009; //5% of maxWNT
            }
        }
    }
}

void diffusion(vector<vector<extra>>& v)
{
    vector<vector<extra>> w;
    exload(w);
    for (int i = 0; i < gridsize; i++)
    {
        for (int j = 0; j < gridsize; j++)
        {
            vector<pair<int, int>>pos;

            if (i % 2 == 0)
            {
                if (j + 1 < gridsize) pos.push_back(make_pair(i, j + 1));
                if (j - 1 >= 0) pos.push_back(make_pair(i, j - 1));
                if (i + 1 < gridsize) pos.push_back(make_pair(i + 1, j));
                if (i + 1 < gridsize && j - 1 >= 0) pos.push_back(make_pair(i + 1, j - 1));
                if (i - 1 >= 0) pos.push_back(make_pair(i - 1, j));
                if (i - 1 >= 0 && j - 1 >= 0) pos.push_back(make_pair(i - 1, j - 1));



            }
            else
            {
                if (j - 1 >= 0) pos.push_back(make_pair(i, j - 1));
                if (j + 1 < gridsize) pos.push_back(make_pair(i, j + 1));
                if (i + 1 < gridsize) pos.push_back(make_pair(i + 1, j));
                if (i + 1 < gridsize && j + 1 < gridsize) pos.push_back(make_pair(i + 1, j + 1));
                if (i - 1 >= 0) pos.push_back(make_pair(i - 1, j));
                if (i - 1 >= 0 && j + 1 < gridsize) pos.push_back(make_pair(i - 1, j + 1));
            }
            for (int n = 0; n < pos.size(); n++)
            {
                w[i][j].pWNT += (v[pos[n].first][pos[n].second].pWNT - v[i][j].pWNT) * dfrate;
            }

        }
    }
    for (int i = 0; i < gridsize; i++)
    {
        for (int j = 0; j < gridsize; j++)
        {
            v[i][j].pWNT += w[i][j].pWNT;
            v[i][j].pWNT = v[i][j].pWNT * pWNTDG; //pWNT degrdation 
            v[i][j].WNT = v[i][j].gradWNT + v[i][j].pWNT;
        }
    }

}
void datanal(vector<vector<intra>>& v, vector<vector<extra>>& w, vector<pair<int, int>>& dtpos, vector<pair<int, int>>& diffpos, vector<pair<int, int>>& divpos, vector<int>& divid)
{
    double lower_bound = 0;
    double upper_bound = 100;
    uniform_real_distribution<double> unif(lower_bound, upper_bound);
    default_random_engine re;
    for (int i = 0; i < gridsize; i++)
    {
        for (int j = 0; j < gridsize; j++)
        {
            pair<int, int>a;
            if (v[i][j].ID != 0)
            {

                if (v[i][j].type == "p" && v[i][j].age >= 16 * 24)
                {
                    if (unif(re) < (v[i][j].age - 16 * 24) * pcdcs) //0% os kiesési esélyel kezdi majd futásonként (óránként) 0.002411%-al növekszik. fv és analízis?
                    {
                        a.first = i;
                        a.second = j;
                        dtpos.push_back(a);
                    }
                }
                if ((v[i][j].type == "e" || v[i][j].type == "g") && v[i][j].age >= 3 * 24)
                {

                    if (unif(re) < (v[i][j].age - 3 * 24) * egcdcs) //0%-on kezd, majd 0.0868%-okkal növekszik
                    {
                        a.first = i;
                        a.second = j;
                        dtpos.push_back(a);
                    }
                }

                if (v[i][j].CYCB < 0.002 && (v[i][j].type == "s" || v[i][j].type == "ap" || v[i][j].type == "sp" || v[i][j].type == "t"))
                {
                    int r = rand() % 10;
                    if (r < 5)
                    {
                        a.first = i;
                        a.second = j;
                        divpos.push_back(a);
                        divid.push_back(v[i][j].ID);
                    }
                }
                if (v[i][j].type == "t")
                {   //5% chance on average
                    if (10 - ((w[i][j].WNT / maxWNT) * 10) > unif(re))
                    {

                        pair<int, int>b;
                        b.first = i;
                        b.second = j;
                        diffpos.push_back(b);
                        for (int k = 0; k < divpos.size(); k++)
                        {
                            if (divpos[k] == b)
                            {
                                divpos.erase(divpos.begin() + k);
                                divid.erase(divid.begin() + k);
                            }
                        }
                    }
                }
            }
        }
    }


}

void death(vector<vector<intra>>& v, vector<pair<int, int>>& dtpos)
{
    intra z = {};
    z.ID = 0;
    for (int i = 0; i < dtpos.size(); i++)
    {
        v[dtpos[i].first][dtpos[i].second] = z;
    }
    vector<pair<int, int>>c;
    dtpos = c;
}

void differentiation(vector<vector<intra>>& v, vector<pair<int, int>>& diffpos)
{

    for (int i = 0; i < diffpos.size(); i++)
    {
        int a = diffpos[i].first;
        int b = diffpos[i].second;
        if (a < gridsize - 1 && b < gridsize - 1 && a>0 && b>0)
        {
            intra nw;
            nw.ID = v[a][b].ID;
            nw.pID = v[a][b].pID;
            int notch = 0;

            if (a % 2)
            {
                if (v[a + 1][b - 1].type == "s" || v[a + 1][b - 1].type == "t") notch += 1;
                if (v[a - 1][b - 1].type == "s" || v[a - 1][b - 1].type == "t") notch += 1;
            }
            if (a % 2)
            {
                if (v[a + 1][b + 1].type == "s" || v[a + 1][b + 1].type == "t") notch += 1;
                if (v[a - 1][b + 1].type == "s" || v[a - 1][b + 1].type == "t") notch += 1;
            }

            if (v[a][b - 1].type == "s" || v[a][b - 1].type == "t") notch += 1;
            if (v[a][b + 1].type == "s" || v[a][b + 1].type == "t") notch += 1;
            if (v[a + 1][b].type == "s" || v[a + 1][b].type == "t") notch += 1;
            if (v[a - 1][b].type == "s" || v[a - 1][b].type == "t") notch += 1;

            if (notch >= 4)
            {
                nw.type = "ap";
            }
            else
            {
                nw.type = "sp";
            }
            v[a][b] = nw;
        }
    }
    vector<pair<int, int>> d;
    diffpos = d;
}

void division(vector<vector<intra>>& v, vector<vector<extra>>& w, vector<pair<int, int>>& divpos, vector<int>& divid)
{
    double lower_bound = 0;
    double upper_bound = 100;
    uniform_real_distribution<double> unif(lower_bound, upper_bound);
    default_random_engine re;

    for (int k = 0; k < divpos.size(); k++)
    {
        intra z = {};
        intra temp;
        int a = divpos[k].first;
        int b = divpos[k].second;
        bool differentiates = false;
        if (v[a][b].ID == divid[k])
        {
            z.type = v[a][b].type;
            z.ID = v[a][b].ID;
            z.pID = v[a][b].pID;
            z.age = v[a][b].age;
            if (v[a][b].type == "s")
            {
                if (a % 2 == 0 && v[a][b - 1].type == "s" && v[a][b + 1].type == "s" && v[a + 1][b - 1].type == "s" && v[a + 1][b].type == "s" && v[a - 1][b - 1].type == "s" && v[a - 1][b].type == "s")
                {
                    temp.type = "s";
                }
                else
                {
                    if (a % 2 == 1 && v[a][b - 1].type == "s" && v[a][b + 1].type == "s" && v[a + 1][b + 1].type == "s" && v[a + 1][b].type == "s" && v[a - 1][b + 1].type == "s" && v[a - 1][b].type == "s")
                    {
                        temp.type = "s";
                    }
                    else
                    {
                        if (w[a][b].WNT >= maxWNT - 3 * WNTPC)
                        {
                            temp.type = "t";
                        }
                        else
                        {
                            differentiates = true;
                        }
                    }
                }
                if (!differentiates)
                {

                    temp.pID = divid[k];
                    temp.ID = id;
                    id++;
                    Simulation::move(a, b, temp, v);
                    temp.ID = divid[k];
                    temp.pID = v[a][b].pID;
                    temp.type = "s";
                    temp.age = v[a][b].age;
                    v[a][b] = temp;
                }
            }
            if (v[a][b].type == "t")
            {
                temp.pID = divid[k];
                temp.ID = id;
                temp.type = "t";
                id++;
                Simulation::move(a, b, temp, v);
                temp.pID = v[a][b].pID;
                temp.ID = v[a][b].ID;
                temp.age = v[a][b].age;
                temp.type = v[a][b].type;
                v[a][b] = temp;
            }
            if (v[a][b].type == "ap")
            {
                if (unif(re) < 50)//fixed 50% E+E
                {

                    z.pID = divid[k];
                    z.ID = id;
                    z.type = "e";
                    z.age = 0;
                    id++;
                    Simulation::move(a, b, z, v);
                    z.pID = v[a][b].pID;
                    z.type = "e";
                    z.age = 0;
                    z.ID = v[a][b].ID;
                    v[a][b] = z;
                }
                else
                {
                    if (unif(re) < w[a][b].BMP)
                    {

                        z.pID = divid[k];
                        z.ID = id;
                        z.age = 0;
                        z.type = "e";
                        id++;
                        Simulation::move(a, b, z, v);
                        temp.ID = divid[k];
                        temp.pID = v[a][b].pID;
                        temp.age = v[a][b].age;
                        temp.type = "ap";
                        v[a][b] = temp;
                    }
                    else
                    {
                        temp.pID = divid[k];
                        temp.ID = id;
                        temp.type = "ap";
                        id++;
                        Simulation::move(a, b, temp, v);
                        temp.pID = v[a][b].pID;
                        temp.ID = divid[k];
                        temp.age = v[a][b].age;
                        temp.type = "ap";
                        v[a][b] = temp;
                    }
                }
            }
            if (v[a][b].type == "sp")
            {
                if (unif(re) < spddd)
                {
                    z.pID = divid[k];
                    z.ID = id;
                    z.type = "g";
                    z.age = 0;
                    id++;
                    Simulation::move(a, b, z, v);

                    z.pID = v[a][b].pID;
                   
                    z.ID = v[a][b].ID;
                    v[a][b] = z;
                }
                else
                {
                    if (unif(re) < spddd)
                    {

                        z.pID = divid[k];
                        z.ID = id;
                        z.age = 0;
                        z.type = "g";
                        id++;
                        Simulation::move(a, b, z, v);
                        temp.type = "sp";
                        temp.pID = v[a][b].pID;
                        temp.ID = v[a][b].ID;
                        temp.age = v[a][b].age;
                        v[a][b] = temp;
                    }
                    else
                    {
                        if (unif(re) < spddd)
                        {
                            temp.pID = divid[k];
                            temp.ID = id;
                            temp.type = "sp";
                            id++;
                            Simulation::move(a, b, temp, v);
                            temp.pID = v[a][b].pID;
                            temp.ID = v[a][b].ID;
                            temp.age = v[a][b].age;
                            v[a][b] = temp;
                        }
                        else
                        {
                            z.pID = divid[k];
                            z.ID = id;
                            z.type = "p";
                            z.age = 0;
                            id++;
                            Simulation::move(a, b, z, v);
                            temp.age = v[a][b].age;
                            temp.ID = v[a][b].ID;
                            temp.pID = v[a][b].pID;
                            temp.type = "sp";
                            v[a][b] = temp;
                        }
                    }
                }
            }
            if (differentiates)
            {
                temp.type = "t";
                temp.ID = divid[k];
                temp.pID = v[a][b].pID;
                temp.age = v[a][b].age;
                v[a][b] = temp;
            }
        }
        else //if dislocated, the cell is put back to the end of the containers with the new coordinates
        {
            for (int i = 0; i < gridsize; i++)
            {
                for (int j = 0; j < gridsize; j++)
                {
                    if (v[i][j].ID == divid[k])
                    {
                        divpos.push_back(make_pair(i, j));
                        divid.push_back(v[i][j].ID);
                    }
                }
            }

        }
    }
    vector<pair<int, int>>divpos1;
    vector<int>divid1;
    divpos = divpos1;
    divid = divid1;
}
namespace Simulation {
    static void move(int a, int b, intra temp, vector<vector<intra>>& v)
    {
        if (temp.ID != 0)
        {
            if (a % 2 == 0 && a < gridsize - 1 && b < gridsize - 1 && b>0 && a>0 && (v[a + 1][b - 1].type == "" || v[a - 1][b - 1].type == "" || v[a][b - 1].type == "" || v[a][b + 1].type == "" || v[a + 1][b].type == "" || v[a - 1][b].type == ""))
            {

                vector<pair<int, int>> frepos;
                if (v[a + 1][b - 1].type == "") frepos.push_back(make_pair(a + 1, b - 1));
                if (v[a - 1][b - 1].type == "") frepos.push_back(make_pair(a - 1, b - 1));
                if (v[a][b - 1].type == "") frepos.push_back(make_pair(a, b - 1));
                if (v[a][b + 1].type == "") frepos.push_back(make_pair(a, b + 1));
                if (v[a + 1][b].type == "") frepos.push_back(make_pair(a + 1, b));
                if (v[a - 1][b].type == "") frepos.push_back(make_pair(a - 1, b));
                int r = rand() % frepos.size();
                v[frepos[r].first][frepos[r].second] = temp;

            }
            else
            {
                if (a % 2 == 1 && a < gridsize - 1 && b < gridsize - 1 && b>0 && a>0 && (v[a + 1][b + 1].type == "" || v[a - 1][b + 1].type == "" || v[a][b - 1].type == "" || v[a][b + 1].type == "" || v[a + 1][b].type == "" || v[a - 1][b].type == ""))
                {
                    vector<pair<int, int>> frepos;
                    if (v[a + 1][b + 1].type == "") frepos.push_back(make_pair(a + 1, b + 1));
                    if (v[a - 1][b + 1].type == "") frepos.push_back(make_pair(a - 1, b + 1));
                    if (v[a][b - 1].type == "") frepos.push_back(make_pair(a, b - 1));
                    if (v[a][b + 1].type == "") frepos.push_back(make_pair(a, b + 1));
                    if (v[a + 1][b].type == "") frepos.push_back(make_pair(a + 1, b));
                    if (v[a - 1][b].type == "") frepos.push_back(make_pair(a - 1, b));
                    int r = rand() % frepos.size();
                    v[frepos[r].first][frepos[r].second] = temp;
                }
                else
                {
                    int newa;
                    int newb;
                    intra newtemp;
                    double y = 0;
                    double x = 0;
                    //direction vectors
                    if (a - gridsize / 2 != 0)
                    {
                        y = (a - gridsize / 2) / abs(a - gridsize / 2);
                    }

                    if (b - gridsize / 2 != 0)
                    {
                        x = (b - gridsize / 2) / abs(b - gridsize / 2);
                    }

                    int r1 = rand() % 3;
                    int r2 = rand() % 2;
                    int r6 = rand() % 6;
                    if (x == -1)
                    {
                        if (y == -1)
                        {
                            if (a % 2 == 0)
                            {
                                if (r1 == 0) newa = a - 1; newb = b;
                                if (r1 == 1) newa = a - 1; newb = b - 1;
                                if (r1 == 2) newa = a; newb = b - 1;
                            }
                            else //a%2==1
                            {
                                if (r1 == 0) newa = a - 1; newb = b;
                                if (r1 == 1) newa = a - 1; newb = b + 1;
                                if (r1 == 2) newa = a; newb = b - 1;
                            }

                        }
                        else
                        {
                            if (y == 0)
                            {
                                if (a % 2 == 0)
                                {
                                    if (r1 == 0 || r1 == 1) { newa = a; newb = b - 1; }
                                    else
                                    {
                                        if (r2 == 1) { newa = a + 1; newb = b - 1; }
                                        else { newa = a - 1; newb = b - 1; }
                                    }
                                }
                                else //a%2==1
                                {
                                    if (r1 == 0 || r1 == 1) { newa = a; newb = b - 1; }
                                    else
                                    {
                                        if (r2 == 1) { newa = a + 1; newb = b; }
                                        else { newa = a - 1; newb = b; }
                                    }
                                }

                            }
                            else //y=1
                            {
                                if (a % 2 == 0)
                                {
                                    if (r1 == 0) newa = a + 1; newb = b - 1;
                                    if (r1 == 1) newa = a + 1; newb = b;
                                    if (r1 == 2) newa = a; newb = b - 1;
                                }
                                else //a%2==1
                                {
                                    if (r1 == 0) newa = a + 1; newb = b;
                                    if (r1 == 1) newa = a + 1; newb = b + 1;
                                    if (r1 == 2) newa = a; newb = b - 1;
                                }

                            }

                        }
                    }
                    else
                    {
                        if (x == 0)
                        {
                            if (y == -1)
                            {
                                if (a % 2 == 0)
                                {
                                    if (r2 == 0) newa = a - 1; newb = b;
                                    if (r2 == 1) newa = a - 1; newb = b - 1;
                                }
                                else //a%2==1
                                {
                                    if (r2 == 0) newa = a - 1; newb = b;
                                    if (r2 == 1) newa = a - 1; newb = b + 1;
                                }

                            }
                            else
                            {
                                if (y == 0)
                                {
                                    if (a % 2 == 0)
                                    {
                                        if (r6 == 0) newa = a; newb = b - 1;
                                        if (r6 == 1) newa = a - 1; newb = b - 1;
                                        if (r6 == 2) newa = a - 1; newb = b;
                                        if (r6 == 3) newa = a; newb = b + 1;
                                        if (r6 == 4) newa = a + 1; newb = b;
                                        if (r6 == 5) newa = a + 1; newb = b - 1;
                                    }
                                    else //a%2=1
                                    {
                                        if (r6 == 0) newa = a; newb = b - 1;
                                        if (r6 == 1) newa = a - 1; newb = b;
                                        if (r6 == 2) newa = a - 1; newb = b + 1;
                                        if (r6 == 3) newa = a; newb = b + 1;
                                        if (r6 == 4) newa = a + 1; newb = b + 1;
                                        if (r6 == 5) newa = a + 1; newb = b;
                                    }
                                }
                                else //y=1
                                {
                                    if (a % 2 == 0)
                                    {
                                        if (r2 == 0) newa = a + 1; newb = b - 1;
                                        if (r2 == 1) newa = a + 1; newb = b;
                                    }
                                    else //a%2==1
                                    {
                                        if (r2 == 0) newa = a + 1; newb = b;
                                        if (r2 == 1) newa = a + 1; newb = b + 1;
                                    }
                                }

                            }

                        }
                        else //x=1
                        {
                            if (y == -1)
                            {
                                if (a % 2 == 0)
                                {
                                    if (r1 == 0) newa = a; newb = b + 1;
                                    if (r1 == 1) newa = a - 1; newb = b;
                                    if (r1 == 2) newa = a - 1; newb = b - 1;
                                }
                                else //a%2==1
                                {
                                    if (r1 == 0) newa = a - 1; newb = b;
                                    if (r1 == 1) newa = a - 1; newb = b + 1;
                                    if (r1 == 2) newa = a; newb = b + 1;
                                }

                            }
                            else
                            {
                                if (y == 0)
                                {
                                    if (a % 2 == 0)
                                    {
                                        if (r1 == 0 || r1 == 1) { newa = a; newb = b + 1; }
                                        else
                                        {
                                            if (r2 == 1) { newa = a + 1; newb = b; }
                                            else { newa = a - 1; newb = b; }
                                        }
                                    }
                                    else //a%2==1
                                    {
                                        if (r1 == 0 || r1 == 1) { newa = a; newb = b + 1; }
                                        else
                                        {
                                            if (r2 == 1) { newa = a + 1; newb = b + 1; }
                                            else { newa = a - 1; newb = b + 1; }
                                        }
                                    }

                                }
                                else //y=1
                                {
                                    if (a % 2 == 0)
                                    {
                                        if (r1 == 0) newa = a + 1; newb = b - 1;
                                        if (r1 == 1) newa = a + 1; newb = b;
                                        if (r1 == 2) newa = a; newb = b + 1;
                                    }
                                    else //a%2==1
                                    {
                                        if (r1 == 0) newa = a + 1; newb = b;
                                        if (r1 == 1) newa = a + 1; newb = b + 1;
                                        if (r1 == 2) newa = a; newb = b + 1;
                                    }

                                }

                            }

                        }
                    }


                    if (newa >= 0 && newb >= 0 && newa < gridsize && newb < gridsize)
                    {
                        newtemp = v[newa][newb];
                        v[newa][newb] = temp;
                        move(newa, newb, newtemp, v);
                    }
                }
            }
        }
    }
    static void pmove(vector<vector<intra>>v)
    {
        double lower_bound = 0;
        double upper_bound = 100;
        uniform_real_distribution<double> unif(lower_bound, upper_bound);
        default_random_engine re;
        vector<int>ids;
        for (int i = 0; i < gridsize; i++)
        {
            for (int j = 0; j < gridsize; j++)
            {
                auto it = find(ids.begin(), ids.end(), v[i][j].ID);
                if (v[i][j].type == "p" && it == ids.end())
                {
                    double y = i - gridsize / 2;
                    double x = j - gridsize / 2;

                    if (unif(re) < sqrt(pow(y, 2) + pow(x, 2)) / maxdistance * 10) //higher distance from origin -> higher chance of switching [0%-10%]
                    {
                        if (y != 0 && x != 0) 
                        {
                            y = y / abs(y);
                            x = x / abs(x);
                            intra temp;
                            temp = v[i - y][j - x];
                            v[i - y][j - x] = v[i][j];
                            v[i][j] = temp;
                        }                        
                    }
                }
            }
        }
    }
    static void show(vector<vector<vector<intra>>> full)
    {
        for (int k = 0; k < full.size(); k++)
        {
            for (int i = 0; i < full[k].size(); i++)
            {
                for (int j = 0; j < full[k][i].size(); j++)
                {
                    cout << full[k][i][j].type << " ";
                    // cout << full[k][i][j].ID<< " ";


                }
                cout << endl;
            }
            cout << endl;
        }
    }
    static void showex(vector<vector<vector<extra>>> fullex)
    {
        for (int k = 0; k < fullex.size(); k++)
        {
            for (int i = 0; i < fullex[k].size(); i++)
            {
                for (int j = 0; j < fullex[k][i].size(); j++)
                {
                    cout << fullex[k][i][j].WNT * 100 << " ";
                    // cout << full[k][i][j].ID<< " ";


                }
                cout << endl;
            }
            cout << endl;
        }
    }
}


class CCanvasWnd : public CWnd
{
public:
    CCanvasWnd(const std::vector<std::vector<intra>>& canvasData)
        : m_canvasData(canvasData) {}

protected:
    afx_msg void OnPaint();
    DECLARE_MESSAGE_MAP()

private:
    const std::vector<std::vector<intra>>& m_canvasData;
};

BEGIN_MESSAGE_MAP(CCanvasWnd, CWnd)
    ON_WM_PAINT()
END_MESSAGE_MAP()

class MyApp : public CWinApp
{
public:
    MyApp(const std::vector<std::vector<std::vector<intra>>>& full)
        : m_data(full) {}

    virtual BOOL InitInstance();

private:
    const std::vector<std::vector<std::vector<intra>>>& m_data;
};

BOOL MyApp::InitInstance()
{
    CWinApp::InitInstance();

    // Register the window class
    CString className = AfxRegisterWndClass(
        CS_HREDRAW | CS_VREDRAW,        // Redraw styles
        ::LoadCursor(NULL, IDC_ARROW),  // Default cursor
        (HBRUSH)::GetStockObject(WHITE_BRUSH), // White background
        NULL                           // No icon
    );
    if (className.IsEmpty())
    {
        DWORD error = ::GetLastError();
        CString errorMsg;
        errorMsg.Format(L"Failed to register window class! Error code: %lu", error);
        AfxMessageBox(errorMsg);
        return FALSE;
    }

    // Create a hidden main window to keep the message loop running
    CWnd* pMainWnd = new CWnd();
    if (!pMainWnd->CreateEx(0, className,
        _T("Main Hidden Window"),
        WS_OVERLAPPEDWINDOW,
        0, 0, 1, 1,   // Small hidden window
        NULL, NULL))
    {
        AfxMessageBox(L"Failed to create main hidden window!");
        delete pMainWnd;
        return FALSE;
    }
    m_pMainWnd = pMainWnd; // Set the hidden window as the main window

    // Create canvas windows
    for (size_t i = 0; i < m_data.size(); ++i)
    {
        CCanvasWnd* pWnd = new CCanvasWnd(m_data[i]);
        if (!pWnd->CreateEx(0, className,
            _T("Canvas Window"),
            WS_OVERLAPPEDWINDOW,
            50 + (int)i * 20, 50 + (int)i * 20, // Adjust starting positions
            560, 560,                          // Window size for 60x60 grid
            NULL, NULL))
        {
            AfxMessageBox(L"Failed to create canvas window!");
            delete pWnd;
            return FALSE;
        }
        pWnd->ShowWindow(SW_SHOW);
    }

    return TRUE; // Return TRUE to start the message loop
}


void CCanvasWnd::OnPaint()
{
    CPaintDC dc(this);

    if (m_canvasData.empty()) return;

    int pixelSize = 9;
    int rows = (int)m_canvasData.size();
    if (rows == 0) return;
    int cols = (int)m_canvasData[0].size();

    for (int row = 0; row < rows; ++row)
    {
        for (int col = 0; col < cols; ++col)
        {
            intra value = m_canvasData[row][col];

            // Get the type as a string
            std::string type = value.type;

            // Map type to color
            COLORREF color;
            if (type == "s")
                color = RGB(0, 255, 255);       // Cyan for "s"
            else if (type == "p")
                color = RGB(0, 255, 0);       // Green for "p"
            else if (type == "ap")
                color = RGB(0, 0, 255);       // Blue for "ap"
            else if (type == "sp")
                color = RGB(255, 255, 0);     // Yellow for "sp"
            else if (type == "g")
                color = RGB(255, 0, 0);     // Red for "g"
            else if (type == "e")
                color = RGB(128, 128, 128);     // Magenta for "e"
            else if (type == "t")
                color = RGB(128, 0, 128);     // Purple for "t"
            else if (type.empty())
                color = RGB(255, 255, 255);         
            else
                color = RGB(128, 0, 128);   

            // Create a brush with the selected color
            CBrush brush(color);
            CBrush* pOldBrush = dc.SelectObject(&brush);

            // Calculate the position of the circle
            int x = col * pixelSize;
            int y = row * pixelSize;

            // Shift every second row by 5 pixels
            if (row % 2 == 1)
                x += 5;

            // Draw the circle
            dc.Ellipse(x, y, x + pixelSize, y + pixelSize);

            // Restore the old brush
            dc.SelectObject(pOldBrush);
        }
    }
}

int main()
{
    vector<vector<vector<intra>>> full;
    vector<vector<vector<extra>>> fullex;

    vector<vector<intra>> ingrid;
    vector<vector<extra>> exgrid;

    inload(ingrid);
    exload(exgrid);

    vector<pair<int, int>> diffpos;
    vector<pair<int, int>> dtpos;
    vector<pair<int, int>> divpos;
    vector<int> divid;

    for (int currun = 0; currun <= 2400; currun++)
    {

        aging(ingrid);
        intraODEs(ingrid, exgrid);
        exocytosis(ingrid, exgrid);
        bool pnt = false;
        for (int i = 0; i < gridsize; i++)
        {
            for (int j = 0; j < gridsize; j++)
            {
                if (exgrid[i][j].pWNT != 0) pnt = true; i = gridsize; j = gridsize;
            }
        }
        if (pnt == true)
        {
            diffusion(exgrid);
        }
        datanal(ingrid, exgrid, dtpos, diffpos, divpos, divid);
        death(ingrid, dtpos);
        differentiation(ingrid, diffpos);
        division(ingrid, exgrid, divpos, divid);
        Simulation::pmove(ingrid);
        //if (currun%100==0)fullex.push_back(exgrid);
        if (currun % 300 == 0)full.push_back(ingrid);

    }
    MyApp theApp(full);
    theApp.InitInstance();
    return theApp.Run();


    return 0;
}