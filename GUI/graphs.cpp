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
    1, true, 0, 1, true},

    {"nradiationdead",
    "Radiation-killed Cells",
    "No. of cells",
    2, true, 0, 1, true},

    {"ndrugdead",
    "Drug-killed Cells",
    "No. of cells",
    3, false, 0, 1, true},

    {"nanoxiadead",
    "Anoxia-killed Cells",
    "No. of cells",
    7, true, 0, 1, true},

    {"nanoxiatagged",
    "Anoxia-tagged Cells",
    "No. of cells",
    8, true, 0, 1, true},

    {"ntagged",
    "Tagged Cells",
    "No. of cells",
    4, false, 0, 1, true},

    {"diameter",
    "Spheroid Diameter",
    "Diameter (um)",
    5, true, 0, 1, true},

    {"volume",
    "Spheroid Volume",
    "Volume (mm3)",
    6, true, 0, 0.001, true},

    {"hypoxicfraction",
    "Hypoxic Fraction",
    "%",
    9, true, 0, 0.1, true},

    {"necroticfraction",
    "Necrotic Fraction",
    "%",
    10, false, 0, 0.1, true}

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
