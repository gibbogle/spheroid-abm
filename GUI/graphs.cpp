#include <qstring.h>
#include "graphs.h"
#include "log.h"

LOG_USE();

Graphs::Graphs()
{
GRAPH_SET tsGraphSet[] = {

    {"nlive",
    "Live Cells",
    "No. of cells",
    "Number of live cells in the blob",
    1, true, 0, 1, 0, TS_TYPE},

    {"nanoxiadead",
    "Anoxia-killed Cells",
    "No. of cells",
     "Total number of cells that have been killed by anoxia",
    2, true, 0, 1, 0, TS_TYPE},

    {"ndrugdead",
    "Drug-killed Cells",
    "No. of cells",
     "Total number of cells that have been killed by drug",
    3, false, 0, 1, 0, TS_TYPE},

    {"nradiationdead",
    "Radiation-killed Cells",
    "No. of cells",
     "Total number of cells that have been killed by radiation",
    4, true, 0, 1, 0, TS_TYPE},

    {"nanoxiatagged",
    "Anoxia-tagged Cells",
    "No. of cells",
     "Current number of cells tagged to die by anoxia",
    5, true, 0, 1, 0, TS_TYPE},

    {"ndrugtagged",
    "Drug-tagged Cells",
    "No. of cells",
     "Current number of cells tagged to die by drug",
    6, false, 0, 1, 0, TS_TYPE},

    {"nradiationtagged",
    "Radiation-tagged Cells",
    "No. of cells",
     "Current number of cells tagged to die by radiation",
    7, false, 0, 1, 0, TS_TYPE},

    {"diameter",
    "Spheroid Diameter",
    "Diameter (um)",
     "Spheroid diameter (um)",
    8, true, 0, 1, 0, TS_TYPE},

    {"volume",
    "Spheroid Volume",
    "Volume (mm3)",
     "Spheroid volume (mm3)",
    9, true, 0, 0.001, 0, TS_TYPE},

    {"hypoxicfraction",
    "Hypoxic Fraction",
    "%",
     "Fraction of cells with oxygen level below the specified threshold for hypoxia",
    10, true, 0, 0.1, 0, TS_TYPE},

    {"growthfraction",
    "Slow growth Fraction",
    "%",
     "Percentage of cells that are growing at a rate less than the specified fraction of the mean growth rate with no nutrient limits",
    11, true, 0, 0.1, 0, TS_TYPE},

    {"necroticfraction",
    "Necrotic Fraction",
    "%",
     "Percentage of the spheroid that is necrotic = (number of vacant sites)/(number of sites taken up by the spheroid)",
    12, true, 0, 0.1, 0, TS_TYPE},

// Profiles
    {"CFSE",
    "CFSE",
    "",
    "CFSE description",
    PROFILE_CFSE, false, 0, 1, 0, PROF_TYPE},

    {"Oxygen",
    "Oxygen",
    "",
    "Oxygen description",
    PROFILE_OXYGEN, true, 0, 1, 0, PROF_TYPE},

    {"Glucose",
    "Glucose",
    "",
    "Glucose description",
    PROFILE_GLUCOSE, false, 0, 1, 0, PROF_TYPE},

    {"Tracer",
    "Tracer",
    "",
    "Tracer description",
    PROFILE_TRACER, false, 0, 1, 0, PROF_TYPE},

    {"TPZdrug",
    "TPZ drug",
    "",
    "TPZ drug description",
    PROFILE_TPZ_DRUG, false, 0, 1, 0, PROF_TYPE},

    {"TPZmetab1",
    "TPZ metabolite 1",
    "",
    "TPZmetab1 description",
    PROFILE_TPZ_METAB_1, false, 0, 1, 0, PROF_TYPE},

    {"TPZmetab2",
    "TPZ metabolite 2",
    "",
    "TPZmetab2 description",
    PROFILE_TPZ_METAB_2, false, 0, 1, 0, PROF_TYPE},

    {"DNBdrug",
    "DNB drug",
    "",
    "DNB drug description",
    PROFILE_DNB_DRUG, false, 0, 1, 0, PROF_TYPE},

    {"DNBmetab1",
    "DNB metabolite 1",
    "",
    "DNBmetab1 description",
    PROFILE_DNB_METAB_1, false, 0, 1, 0, PROF_TYPE},

    {"DNBmetab2",
    "DNB metabolite 2",
    "",
    "DNBmetab2 description",
    PROFILE_DNB_METAB_2, false, 0, 1, 0, PROF_TYPE},

    {"growthrate",
    "Growth rate",
    "",
    "Growth rate description",
    PROFILE_GROWTH_RATE, false, 0, 1, 0, PROF_TYPE},

    {"cellvolume",
    "Cell volume",
    "",
    "Cell volume description",
    PROFILE_CELL_VOLUME, false, 0, 1, 0, PROF_TYPE},
};

    n_tsGraphs = sizeof(tsGraphSet)/sizeof(GRAPH_SET);
    tsGraphs = new GRAPH_SET[n_tsGraphs];
    for (int i=0; i<n_tsGraphs; i++) {
        tsGraphs[i] = tsGraphSet[i];
    }
    graphList = new GRAPH_SET[maxGraphs];
    nGraphs = maxGraphs;
}


GRAPH_SET Graphs::get_graph(int k)
{
	return graphList[k];
}

int Graphs::get_dataIndex(int k)
{
	return graphList[k].dataIndex;
}

QString Graphs::get_tag(int k)
{
	return graphList[k].tag;
}

QString Graphs::get_title(int k)
{
	return graphList[k].title;
}

QString Graphs::get_yAxisTitle(int k)
{
	return graphList[k].yAxisTitle;
}

QString Graphs::get_description(int k)
{
    return graphList[k].description;
}

double Graphs::get_maxValue(int k) {
	return graphList[k].maxValue;
}

double Graphs::get_scaling(int k) {
	return graphList[k].scaling;
}

double Graphs::get_yscale(int k) {
    return graphList[k].yscale;
}

double Graphs::get_xscale(double xmax) {
    int n = 1;
    for (;;) {
        if (xmax <= n) break;
        n++;
    }
    return double(n);
}


bool Graphs::isActive(int k)
{
	return graphList[k].active;
}

int Graphs::get_type(int k) {
    return graphList[k].type;
}

bool Graphs::isTimeseries(int k)
{
    return (graphList[k].type == TS_TYPE);
}

bool Graphs::isProfile(int k)
{
    return (graphList[k].type == PROF_TYPE);
}

bool Graphs::isDistribution(int k)
{
    return (graphList[k].type == DIST_TYPE);
}

void Graphs::set_maxValue(int k, double v)
{
	graphList[k].maxValue = v;
}

void Graphs::makeGraphList(int non_ts)
{
    int k = maxGraphs;
    int nts = 0;
    for (int i=0; i<n_tsGraphs; i++) {
        if (tsGraphs[i].active) {
            k--;
            graphList[k] = tsGraphs[i];
            nts++;
            if (nts == maxGraphs - non_ts) break;
        }
    }
    int ndummy = maxGraphs - nts - non_ts;
    for (k=0; k<ndummy; k++) {
        graphList[k].tag = "dummy";
        graphList[k].active = false;
//        graphList[k].ts = true;
        graphList[k].type = TS_TYPE;
        graphList[k].scaling = 1;
    }
    for (k=ndummy; k<ndummy + non_ts; k++) {
        graphList[k].tag = "non_ts";
        graphList[k].active = true;
//        graphList[k].ts = false;
        graphList[k].type = DIST_TYPE;  //????
        graphList[k].scaling = 1;
    }
    nGraphs = maxGraphs;

    char msg[128];
    sprintf(msg,"nGraphs: %d  non_ts: %d  nts: %d",nGraphs,non_ts,nts);
    LOG_MSG(msg);
//    for (k=0; k<nGraphs; k++) {
//        LOG_QMSG(graphList[k].tag);
//        sprintf(msg,"k: %d scaling: %f",k,graphList[k].scaling);
//        LOG_MSG(msg);
//    }
}

