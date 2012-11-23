#include "field.h"
#include "log.h"

LOG_USE();

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
Field::Field(QWidget *page2D)
{
	field_page = page2D;
    axis = X_AXIS;
    fraction = 0;
    constituent = OXYGEN;
    slice_changed = true;
}

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
Field::~Field()
{
}

//------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------
void Field::setConstituent(QAbstractButton *button)
{
//    QMessageBox msgBox;
//    msgBox.setText("setConstituent");
//    msgBox.exec();
    LOG_MSG("setConstituent");
	int prev_constituent = constituent;
    QString text = button->text();
    if (text.compare("Oxygen") == 0)
        constituent = OXYGEN;
    else if (text.compare("Glucose") == 0)
        constituent = GLUCOSE;
    else if (text.compare("Drug A") == 0)
        constituent = DRUG_A;
    else if (text.compare("Drug B") == 0)
        constituent = DRUG_B;
	if (constituent != prev_constituent) {
		constituent_changed = true;
		displayField();
	}
}

//------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------
void Field::setPlane(QAbstractButton *button)
{
//    QMessageBox msgBox;
//    msgBox.setText("setPlane");
//    msgBox.exec();
    LOG_MSG("setPlane");
    QString text = button->text();
//    if (button->isChecked()) {
//        LOG_QMSG(text);
//    } else {
//        LOG_MSG("Button not checked");
//    }
	int prev_axis = axis;
    if (text.compare("X-Y") == 0)
        axis = Z_AXIS;
    else if (text.compare("Y-Z") == 0)
        axis = X_AXIS;
    else if (text.compare("X-Z") == 0)
        axis = Y_AXIS;
	if (axis != prev_axis) {
		slice_changed = true;
		displayField();
	}
}

//------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------
void Field::setFraction(QString text)
{
	double prev_fraction = fraction;
	fraction = text.toDouble();
	if (fraction != prev_fraction) {
		slice_changed = true;
		displayField();
	}
}

//------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------
void Field::chooseParameters()
{
//    chemo_select[0] = 1;
//    chemo_select[1] = 0;
//    chemo_select[2] = 0;
//    chemo_select[3] = 0;
//    chemo_displayed[0] = false;
//    chemo_displayed[1] = false;
//    chemo_displayed[2] = false;
//    chemo_displayed[3] = false;

    bool ok;
    QString text;
    QStringList items;
//    items << tr("X-Y plane") << tr("X-Z plane") << tr("Y-Z plane");
//    items << QString("X-Y plane") << QString("X-Z plane") << QString("Y-Z plane");

//    QString item = QInputDialog::getItem(this,tr("QInputDialog::getItem()"),
//                                         tr("Slice plane:"), items, 0, false, &ok);
    QString item = "";
    if (ok && !item.isEmpty()) {
        if (item.contains("X-Y")) {
            axis = 3;
            text = "X-Y plane";
        } else if (item.contains("X-Z")) {
            axis = 2;
            text = "X-Z plane";
        } else {
            axis = 1;
            text = "Y-Z plane";
        }
 //       this->ui->label_plane->setText(text);
    }

//    double d = QInputDialog::getDouble(this, tr("QInputDialog::getDouble()"),
//                                       tr("Slice fractional position:"), 0.0, -1, 1, 2, &ok);
    double d = 0;
    fraction = d;
    sprintf(msg,"Radius fraction: %5.2f",fraction);
    text = QString(msg);
//    this->ui->label_fraction->setText(text);

//    scale = QInputDialog::getDouble(this, tr("QInputDialog::getDouble()"),
//                                       tr("Scaling factor: (0 normalizes scale)"), 0.0, 0, 100, 2, &ok);
//    if (scale == 0) {
//        text = "Normalized vectors";
//    } else {
//        sprintf(msg,"Vector scaling: %5.2f",scale);
//        text = QString(msg);
//    }
//    this->ui->label_scaling->setText(text);

//    int ret = QMessageBox::question(this, tr("Chemokine relative strength"),
//                                        tr("Do you want to multiply gradients by the chemokine strength?"),
//                                        QMessageBox::Yes | QMessageBox::No, QMessageBox::No);
//    if (ret == QMessageBox::Yes) {
//        use_strength = 1;
//        text = "Using relative strength";
//    } else {
//        use_strength = 0;
//        text = "Not using relative strength";
//    }
//    this->ui->label_strength->setText(text);
}


//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
void Field::displayField()
{
//    QGraphicsScene* scene = new QGraphicsScene(QRect(0, 0, 130, 280));
    QGraphicsScene* scene = new QGraphicsScene(QRect(0, 0, 800, 800));
    QBrush brush;
    QGraphicsTextItem *text;
    int nsites, nconst;

	if (slice_changed) {
		get_fieldinfo(&axis, &fraction, &nsites, &nconst);
		sprintf(msg,"nsites: %d",nsites);
		LOG_MSG(msg);
		this->data = (FIELD_DATA *)malloc(nsites*sizeof(FIELD_DATA));
		get_fielddata(&nsites, this->data);
		slice_changed = false;
	}

    brush.setColor(QColor(150,100,0));
    brush.setStyle(Qt::SolidPattern);
    scene->addEllipse(10,10,20,20,Qt::NoPen, brush);
    text = scene->addText("FDC");
    text->setPos(35, 10);

    brush.setColor(QColor(200,60,100));
    brush.setStyle(Qt::SolidPattern);
    scene->addEllipse(10,40,20,20,Qt::NoPen, brush);
    text = scene->addText("MRC");
    text->setPos(35, 40);

    brush.setColor(QColor(30,20,255));
    scene->addEllipse(10,70,20,20,Qt::NoPen, brush);
    text = scene->addText("Naive B cell");
    text->setPos(35, 70);

    brush.setColor(QColor(0,200,255));
    scene->addEllipse(10,100,20,20,Qt::NoPen, brush);
    text = scene->addText("CCR7 UP");
    text->setPos(35, 100);

    brush.setColor(QColor(50,255,150));
    scene->addEllipse(10,130,20,20,Qt::NoPen, brush);
    text = scene->addText("EBI2 UP");
    text->setPos(35, 130);

    brush.setColor(QColor(255,255,0));
    scene->addEllipse(10,160,20,20,Qt::NoPen, brush);
    text = scene->addText("BCL6 HI");
    text->setPos(35, 160);

    brush.setColor(QColor(0,150,0));
    scene->addEllipse(10,190,20,20,Qt::NoPen, brush);
    text = scene->addText("BCL6 LO");
    text->setPos(35, 190);

    brush.setColor(QColor(128,128,128));
    scene->addEllipse(10,220,20,20,Qt::NoPen, brush);
    text = scene->addText("Max divisions");
    text->setPos(35, 220);

    brush.setColor(QColor(255,0,0));
    scene->addEllipse(10,250,20,20,Qt::NoPen, brush);
    text = scene->addText("Plasma cell");
    text->setPos(35, 250);

    brush.setColor(QColor(255,130,0));
    scene->addEllipse(10,280,20,20,Qt::NoPen, brush);
    text = scene->addText("CD4 T cell");
    text->setPos(35, 280);

    QGraphicsView* view = new QGraphicsView(field_page);
    view->setScene(scene);
    view->setGeometry(QRect(0, 0, 810, 810));
    view->show();
}
