#ifndef GRAPHS_H
#define GRAPHS_H

struct graph_set {
	QString tag;
	QString title;
	QString yAxisTitle;
	int dataIndex;		// this must be consistent with the ordering of summary_data[]
	bool active;		// false for dummy graphs
	double maxValue;
	double scaling;		// multiplying factor for scaling of summary_data
};

typedef graph_set GRAPH_SET;

class Graphs
{
	GRAPH_SET *graphList;

public:

	Graphs();
	~Graphs();
	GRAPH_SET get_graph(int);
	int nGraphs;
	int get_dataIndex(int);
	QString get_tag(int);
	QString get_title(int);
	QString get_yAxisTitle(int);
	double get_maxValue(int);
	double get_scaling(int);
	bool isActive(int);
	void set_maxValue(int, double);

};

#endif
