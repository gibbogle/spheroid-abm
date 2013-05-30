#include <qstring.h>
#include "graphs.h"
#include "log.h"

LOG_USE();

// summaryData(1:13) = (/ istep, Ncells, Nanoxia_dead, Ndrug_dead, Nradiation_dead, &
//  Ntagged_anoxia, Ntagged_drug, Ntagged_radiation, &
//	diam_um, vol_mm3_1000, hypoxic_percent_10, growth_percent_10, necrotic_percent_10 /)

Graphs::Graphs()
{
GRAPH_SET tsGraphSet[] = {

    {"nlive",
    "Live Cells",
    "No. of cells",
    1, true, 0, 1, true},

    {"nanoxiadead",
    "Anoxia-killed Cells",
    "No. of cells",
    2, true, 0, 1, true},

    {"ndrugdead",
    "Drug-killed Cells",
    "No. of cells",
    3, false, 0, 1, true},

    {"nradiationdead",
    "Radiation-killed Cells",
    "No. of cells",
    4, true, 0, 1, true},

    {"nanoxiatagged",
    "Anoxia-tagged Cells",
    "No. of cells",
    5, true, 0, 1, true},

    {"ndrugtagged",
    "Drug-tagged Cells",
    "No. of cells",
    6, false, 0, 1, true},

    {"nradiationtagged",
    "Radiation-tagged Cells",
    "No. of cells",
    7, false, 0, 1, true},

    {"diameter",
    "Spheroid Diameter",
    "Diameter (um)",
    8, true, 0, 1, true},

    {"volume",
    "Spheroid Volume",
    "Volume (mm3)",
    9, true, 0, 0.001, true},

    {"hypoxicfraction",
    "Hypoxic Fraction",
    "%",
    10, true, 0, 0.1, true},

    {"growthfraction",
    "Growth Fraction",
    "%",
    11, true, 0, 0.1, true},

    {"necroticfraction",
    "Necrotic Fraction",
    "%",
    12, false, 0, 0.1, true}

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

double Graphs::get_maxValue(int k) {
	return graphList[k].maxValue;
}

double Graphs::get_scaling(int k) {
	return graphList[k].scaling;
}

bool Graphs::isActive(int k)
{
	return graphList[k].active;
}

bool Graphs::isTimeseries(int k)
{
    return graphList[k].ts;
}

void Graphs::set_maxValue(int k, double v)
{
	graphList[k].maxValue = v;
}

void Graphs::makeGraphList(int non_ts)
{
//    char msg[128];
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
//    sprintf(msg,"nts: %d  ndummy: %d",nts,ndummy);
//    LOG_MSG(msg);
    for (k=0; k<ndummy; k++) {
        graphList[k].active = false;
        graphList[k].ts = true;
        graphList[k].tag = "dummy";
        graphList[k].scaling = 1;
    }
    for (k=ndummy; k<ndummy + non_ts; k++) {
        graphList[k].tag = "non_ts";
        graphList[k].active = true;
        graphList[k].ts = false;
        graphList[k].scaling = 1;
    }
    nGraphs = maxGraphs;
//    sprintf(msg,"nGraphs: %d",nGraphs);
//    LOG_MSG(msg);
//    for (k=0; k<nGraphs; k++) {
//        LOG_QMSG(graphList[k].tag);
//        sprintf(msg,"k: %d scaling: %f",k,graphList[k].scaling);
//        LOG_MSG(msg);
//    }
}
