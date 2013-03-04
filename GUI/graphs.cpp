#include <qstring.h>
#include "graphs.h"
#include "log.h"

LOG_USE();


// summaryData(1:9) = (/ int(tnow/60),istep,ntot,ncogseed,ncog,Ndead,int(InflowTotal*60/DELTA_T), int(100*vascularity), teffgen/)

// summaryData(1:4) = (/ istep, nlive, ndead, int(radius in um) /)
// summaryData(1:6) = (/ istep, Ncells, Nradiationdead, Ndrugdead, Ntagged, diam_um, vol_mm3, Nanoxiadead /)

Graphs::Graphs()
{
GRAPH_SET tsGraphSet[] = {

    {"nlive",
    "Number of Live Cells",
    "No. of cells",
    1, true, 0, 1, true},

    {"nradiationdead",
    "Number of Radiation-killed Cells",
    "No. of cells",
    2, true, 0, 1, true},

    {"ndrugdead",
    "Number of Drug-killed Cells",
    "No. of cells",
    3, true, 0, 1, true},

    {"nanoxiadead",
    "Number of Anoxia-killed Cells",
    "No. of cells",
    7, true, 0, 1, true},

    {"ntagged",
    "Number of Tagged Cells",
    "No. of cells",
    4, false, 0, 1, true},

    {"diameter",
    "Spheroid Diameter",
    "Diameter (um)",
    5, false, 0, 1, true},

    {"volume",
    "Spheroid Volume",
    "Volume (mm3)",
    6, true, 0, 0.001, true}

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
    char msg[128];
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
