#include <QtGui>

#include "mainwindow.h"
#include "log.h"
#include "params.h"
//#include "graphs.h"
//#include "misc.h"
//#include "plot.h"
//#include "myvtk.h"
//#include "field.h"
//#include "transfer.h"
//#include "dialog.h"

#include "global.h"

#ifdef linux
#include <QTcpServer>
#else
#include <QTcpServer.h>
#endif

LOG_USE();

extern Params *parm;	// I don't believe this is the right way, but it works

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
void MainWindow::on_comb_TPZ_currentIndexChanged(int index)
{
    text_TPZ_DRUG_NAME->setText(comb_TPZ->currentText());
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
void MainWindow::on_comb_DNB_currentIndexChanged(int index)
{
    text_DNB_DRUG_NAME->setText(comb_DNB->currentText());
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
void MainWindow::on_cbox_USE_TPZ_DRUG_toggled(bool checked)
{
    LOG_MSG("cbox_use_TPZ_drug toggled");
    QLineEdit *leb = findChild<QLineEdit *>("line_TPZ_DRUG_BDRY_CONC");
//    QCheckBox *cbd = findChild<QCheckBox *>("cbox_TPZ_DRUG_DECAY");
    QCheckBox *cbm = findChild<QCheckBox *>("cbox_TPZ_DRUG_SIMULATE_METABOLITE");
    leb->setEnabled(checked);
    cbm->setEnabled(checked);
//    setTreatmentFileUsage();
    comb_TPZ->setEnabled(checked);
    text_TPZ_DRUG_NAME->setEnabled(checked);
    text_TPZ_DRUG_NAME->setText(comb_TPZ->currentText());
//    int indexTPZ = comb_TPZ->currentIndex();
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
void MainWindow::on_cbox_USE_DNB_DRUG_toggled(bool checked)
{
//    LOG_MSG("cbox_use_DNB_drug toggled");
    QLineEdit *leb = findChild<QLineEdit *>("line_DNB_DRUG_BDRY_CONC");
//    QCheckBox *cbd = findChild<QCheckBox *>("cbox_DNB_DRUG_DECAY");
    QCheckBox *cbm = findChild<QCheckBox *>("cbox_DNB_DRUG_SIMULATE_METABOLITE");
    leb->setEnabled(checked);
    cbm->setEnabled(checked);
//    setTreatmentFileUsage();
    comb_DNB->setEnabled(checked);
    text_DNB_DRUG_NAME->setEnabled(checked);
    text_DNB_DRUG_NAME->setText(comb_DNB->currentText());
//    int indexDNB = comb_DNB->currentIndex();
}
//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
void MainWindow::on_checkBox_CELLDISPLAY_1_toggled(bool display)
{
    vtk->display_celltype[1] = display;
//    vtk->cleanup();
    vtk->renderCells();
//    LOG_QMSG("toggled display_celltype[1]");
}

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
void MainWindow::on_checkBox_CELLDISPLAY_2_toggled(bool display)
{
    vtk->display_celltype[2] = display;
    vtk->renderCells();
}

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
void MainWindow::on_comboBox_CELLCOLOUR_1_currentIndexChanged(int index)
{
    QColor qcolor;
//    vtk->celltype_colour[1] = comboBox_CELLCOLOUR_1->currentText();
    qcolor = comboColour[index];
    vtk->celltype_colour[1] = qcolor;
    vtk->renderCells();
    sprintf(msg,"changed celltype_colour[1]: index: %d r,g,b: %d %d %d",index,qcolor.red(),qcolor.green(),qcolor.blue());
    LOG_MSG(msg);
}

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
void MainWindow::on_comboBox_CELLCOLOUR_2_currentIndexChanged(int index)
{
//    vtk->celltype_colour[2] = comboBox_CELLCOLOUR_2->currentText();
    vtk->celltype_colour[2] = comboColour[index];
    vtk->renderCells();
}

//------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------
void MainWindow::on_checkBox_FACS_PLOT_toggled(bool checked)
{
    line_FACS_INTERVAL->setEnabled(checked);
    if (!checked) {
        line_FACS_INTERVAL->setText("0");
    }
}

//-------------------------------------------------------------
// Switches to the FACS screen
//-------------------------------------------------------------
void MainWindow::on_action_FACS_triggered()
{
    stackedWidget->setCurrentIndex(4);
    action_outputs->setEnabled(true);
    action_inputs->setEnabled(true);
    action_VTK->setEnabled(true);
    action_FACS->setEnabled(false);
    Global::showingFACS = true;
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
//void MainWindow::on_cbox_USE_RADIATION_toggled(bool checked)
//{
//    setTreatmentFileUsage();
//}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
//void MainWindow::setTreatmentFileUsage()
//{
//    if (cbox_USE_TPZ_DRUG->isChecked() || cbox_USE_DNB_DRUG->isChecked() || cbox_USE_RADIATION->isChecked()) {
//        enableUseTreatmentFile();
//    } else {
//        disableUseTreatmentFile();
//    }
//}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
void MainWindow::on_line_CELLPERCENT_1_textEdited(QString pc1_str)
{
    double pc1 = pc1_str.toDouble();
    double pc2 = 100 - pc1;
    QString pc2_str = QString::number(pc2);
    line_CELLPERCENT_2->setText(pc2_str);
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
void MainWindow::on_line_CELLPERCENT_2_textEdited(QString pc2_str)
{
    double pc2 = pc2_str.toDouble();
    double pc1 = 100 - pc2;
    QString pc1_str = QString::number(pc1);
    line_CELLPERCENT_1->setText(pc1_str);
}


//------------------------------------------------------------------------------------------------------
// This should be used for any radioButtonGroups for model input parameters
//------------------------------------------------------------------------------------------------------
void MainWindow::radioButtonChanged(QAbstractButton *b)
{
    QString wtag = b->objectName();
    int rbutton_case;
    if (b->isChecked()) {
        QString ptag = parse_rbutton(wtag,&rbutton_case);
        // Now need to reflect the change in the workingParameterList
        // Need to locate ptag
        for (int k=0; k<nParams; k++) {
            PARAM_SET p = parm->get_param(k);
            if (ptag.compare(p.tag) == 0) {
                parm->set_value(k,double(rbutton_case));
                break;
            }
        }
    }
}

void MainWindow::buttonClick_constituent(QAbstractButton* button)
{
    LOG_MSG("buttonClick_constituent");
    field->setConstituent(button);
}

void MainWindow::buttonClick_plane(QAbstractButton* button)
{
    LOG_MSG("buttonClick_plane");
    field->setPlane(button);
}

void MainWindow::buttonClick_canvas(QAbstractButton* button)
{
    LOG_MSG("buttonClick_canvas");
}

void MainWindow::textChanged_fraction(QString text)
{
    LOG_MSG("textChanged_fraction");
    field->setFraction(text);
}

void MainWindow::textEdited_fraction(QString text)
{
    LOG_MSG("textEdited_fraction");
    field->setFraction(text);
}

void MainWindow::onSelectConstituent()
{
    if (exthread != NULL)
        field->selectConstituent();
}

void MainWindow::on_verticalSliderTransparency_sliderMoved(int position)
{
    vtk->setOpacity(position);
}

