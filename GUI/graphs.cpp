#include <qstring.h>
#include "graphs.h"

// summaryData(1:9) = (/ int(tnow/60),istep,ntot,ncogseed,ncog,Ndead,int(InflowTotal*60/DELTA_T), int(100*vascularity), teffgen/)

// summaryData(1:4) = (/ istep, nlive, ndead, int(radius in um) /)
// summaryData(1:6) = (/ istep, Ncells, Ndead, Ndrugdead, Ntagged, diam_um /)

Graphs::Graphs()
{
GRAPH_SET graphs[] = {

    {"dummy1",
    "",
    "",
    0, false, 0, 1},

    {"dummy2",
    "",
    "",
    0, false, 0, 1},

    {"nlive",
    "Number of Live Cells",
    "No. of cells",
    1, true, 0, 1},

    {"nradiationdead",
    "Number of Radiation-killed Cells",
    "No. of cells",
    2, true, 0, 1},

    {"ndrugdead",
    "Number of Drug-killed Cells",
    "No. of cells",
    3, true, 0, 1},

    {"ntagged",
    "Number of Tagged Cells",
    "No. of cells",
    4, true, 0, 1},

    {"diameter",
    "Spheroid Diameter",
    "Diameter (um)",
    5, true, 0, 1},

};


//    diam_number = 7;
	nGraphs = sizeof(graphs)/sizeof(GRAPH_SET);
	graphList = new GRAPH_SET[nGraphs];
	for (int i=0; i<nGraphs; i++) {
		graphList[i] = graphs[i];
	}
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

void Graphs::set_maxValue(int k, double v)
{
	graphList[k].maxValue = v;
}
