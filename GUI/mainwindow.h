/****************************************************************************

****************************************************************************/

#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <string>
#include <fstream>
#include <QTcpServer>
#include <QTcpSocket>

using namespace std;

//#ifdef __DISPLAY768
//#include "ui_spheroid_GUI-768.h"
//#else
#include "ui_spheroid_GUI.h"
//#endif
#include <qwt_plot_curve.h>
#include "params.h"
#include "misc.h"
#include "plot.h"
#include "myvtk.h"
#include "field.h"
#include "result_set.h"
#include "log.h"
#include "SimpleView3DUI.h"
#include "SimpleView2DUI.h"

QT_BEGIN_NAMESPACE
class QAction;
class QMenu;
class QPlainTextEdit;
class QMdiArea;
class QTcpServer;
class QTcpSocket;
QT_END_NAMESPACE

class SliderPlus 
{
	QString name;
	int pindex;
	int windex;
	double vmin;
	double vmax;
	double dv;
	int n;

public:

	SliderPlus(QString, double, double, int,int, int);
	~SliderPlus();
	int val_to_int(double);
	double int_to_val(int);
	QString val_to_str(double);
	double str_to_val(QString);
	int pIndex();
	int wIndex();
	int nTicks();
};


class MainWindow : public QMainWindow, private Ui::MainWindow
{
    Q_OBJECT

public:
    MainWindow(QWidget *parent = 0);

	char msg[2048];

protected:
    void closeEvent(QCloseEvent *event);

private slots:
    void on_action_show_gradient2D_triggered();
    void on_action_show_gradient3D_triggered();
    void newFile();
    void open();
    void about();
    void documentWasModified();

    bool save();
    bool saveAs();
	void readInputFile();
	void loadResultFile();
    void goToInputs();
    void goToOutputs();
    void goToVTK();
    void goToField();
    void runServer();
    void pauseServer();
    void stopServer();
	void changeParam();
	void redrawDistPlot();
	void showMore(QString);
	void updateSliderBox();

	double getMaximum(RESULT_SET *, double *);
	void addGraph();
	void removeGraph();
	void removeAllGraphs();
	void playVTK();
	void setVTKSpeed();
	void saveSnapshot();
    void showGradient3D();
    void showGradient2D();
    void setSavePosStart();

    void on_radioButton_oxygen_clicked();
    void on_radioButton_glucose_clicked(bool checked);

    void onSelectConstituent();

public slots:
	void preConnection();
	void outputData(QString);
	void postConnection();
	void timer_update();
	void errorPopup(QString);
	void displayScene();
//    void displayFld();
    void showSummary(int);
    void startRecorder();
    void stopRecorder();
 //   void on_buttonButton_constituent_clicked();
    void buttonClick_constituent(QAbstractButton* button);
    void buttonClick_plane(QAbstractButton* button);
	void textChanged_fraction(QString text);
	void textEdited_fraction(QString text);
    void setupConc(int nc, bool *used);

    void on_cbox_use_drugA_toggled(bool checked);
//    void on_cbox_drugA_decay_toggled(bool checked);
    void on_cbox_drugA_metabolite_toggled(bool checked);
//    void on_cbox_drugA_metabolite_decay_toggled(bool checked);
    void on_cbox_use_drugB_toggled(bool checked);
//    void on_cbox_drugB_decay_toggled(bool checked);
    void on_cbox_drugB_metabolite_toggled(bool checked);
//    void on_cbox_drugB_metabolite_decay_toggled(bool checked);

private:
    void createActions();
	void createLists();
	void drawDistPlots();
	void setupParamList();
	void loadParams();
	void reloadParams();

    void trackError();

    void enableUseOxygen();
    void disableUseOxygen();
    void enableUseGlucose();
    void disableUseGlucose();
    void enableUseSN30K();
    void disableUseSN30K();
    void enableUseDrugA();
    void disableUseDrugA();
    void enableUseDrugB();
    void disableUseDrugB();
    void enableUseTreatmentFile();
    void disableUseTreatmentFile();

	void writeout();
	void execute_para();
	void init_VTK();
	void read_cell_positions();
	void close_sockets();
	void compareOutputs();
	void clearAllGraphs();
	void initializeGraphs(RESULT_SET *);
	void drawGraphs();
	QString selectResultSet();
	int selectGraphCase();
    void setupGraphSelector();
    void setGraphsActive();

    void showBool(QString, bool);

	double erf(double z);
    double pnorm(double x1, double x2, double mu, double sig);
    double plognorm(double x1, double x2, double mu, double sig);
    void create_lognorm_dist(double p1, double p2,int n, double *x, double *prob);
	int dist_limit(double *p, int n);
	QString parse_rbutton(QString wtag, int *rbutton_case);
	void setBdryRadioButton(QRadioButton *w_rb, int val);
	void setLineEditVisibility(QString wname, int val);

	PARAM_SET get_param(int);

    void createMenus();
    void createToolBars();
    void createStatusBar();
    void readSettings();
    void writeSettings();
    bool maybeSave();
    void loadFile(const QString &fileName);
//    bool saveFile(const QString &fileName);
    void setCurrentFile(const QString &fileName);
    QString strippedName(const QString &fullFileName);

    QPlainTextEdit *textEdit;
    QString curFile;
	QList<QLineEdit *> lineEdit_list;
	QList<QSpinBox *> spin_list;
	QList<QComboBox *> combo_list;
	QList<QCheckBox *> checkbox_list;
	QList<QRadioButton *> radiobutton_list;
	QList<QSlider *> slider_list;
	QList<QLabel *> label_list;
	QList<SliderPlus *> sliderplus_list;
	QList<QWidget *> sliderParam;
	QList<RESULT_SET *> result_list;

	QwtPlot *distplot_list[5];
	QwtPlotCurve *curve_list[5];

	QList<QWidget *> widget_list;

	int nDistPts;
	int nTicks;
	int nParams;
	int nSliders;
	int nWidgets;
	int nLabels;
	int *param_to_sliderIndex;
	bool paramSaved;
	bool paused;
	bool posdata;
	bool DCmotion;
	bool done;
	bool first;
	bool started;
	bool firstVTK;
	bool playingVTK;
	int tickVTK;
	int currentDescription;
	QString defaultInputFile;
	QString inputFile;
//	QString stopfile;
//	QString pausefile;
	QString cellfile;
//	QString dll_path;
	QString vtkfile;
	QTextBrowser *box_outputData;
	SocketHandler *sthread0;
	SocketHandler *sthread1;
	QTimer *timer;

	int step;
	int ntimes;
	int savepos_start;
	int ncpu;
	double hours;
	double hour;
	int progress;
	int nGraphs;		// act, ntot_LN, ncog_PER, ...
	int nGraphCases;

	RESULT_SET *newR;

//	Plot *graph_act;
//	Plot *graph_ntot_LN;
//	Plot *graph_ncog_PER;
//	Plot *graph_ncog_LN;
//	Plot *graph_ncog;
//	Plot *graph_ncogseed;
//	Plot *graph_nDC;
//	Plot *graph_teffgen;
//	Plot *graph_nbnd;
	Plot *graph_dummy;	// placeholder

    Plot *pGraph[32];
    QCheckBox *checkBox_conc;
    QCheckBox *checkBox_vol;
    QCheckBox *checkBox_oxy;
    QCheckBox **cbox_ts;

	QString graphCaseName[Plot::ncmax];
	RESULT_SET *graphResultSet[Plot::ncmax];
	static const bool show_outputdata = false;
	static const bool use_CPORT1 = false;

	static const int CPORT0 = 5000;
	static const int CPORT1 = 5001;
	static const bool USE_RANGES = false;

	MyVTK *vtk;
    Field *field;
	ExecThread *exthread;
};

class MyDoubleValidator : public QDoubleValidator
{
public:
	MyDoubleValidator( double bottom, double top, int decimals, QObject* parent = 0)
		: QDoubleValidator( bottom, top, decimals, parent)
	{}

	QValidator::State validate ( QString &input, int &pos ) const
	{
		if ( input.isEmpty() || input == "." ) {
			return Intermediate;
		}
		bool ok;
		double entered = input.toDouble(&ok);
		if (!ok) return Invalid;
		if (entered < bottom())
			return Intermediate;
		if ( QDoubleValidator::validate( input, pos ) != Acceptable ) {
			return Invalid;
		}
		return Acceptable;
	}
};

//static const double DELTA_T = 0.25;

#endif
