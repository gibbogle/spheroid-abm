#include "field.h"
#include "log.h"

LOG_USE();

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
Field::Field(QWidget *page2D)
{
	field_page = page2D;
    MAX_CHEMO = 4;
    axis = Z_AXIS;
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
    LOG_MSG("setFraction");
	fraction = text.toDouble();
	if (fraction != prev_fraction) {
		slice_changed = true;
		displayField();
	}
}

//------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------
void Field::setSliceChanged()
{
    slice_changed = true;
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
    QGraphicsScene* scene = new QGraphicsScene(QRect(0, 0, 690, 690));
    QBrush brush;
    QGraphicsTextItem *text;
    int i, xindex, yindex, ix, iy, cp, dp, c, rmax, w, xp0, rgbcol[3];
    double xp, yp, d0, d, volume, scale;
    double a, b;

    LOG_MSG("displayField");
	if (slice_changed) {
        get_fieldinfo(&NX, &axis, &fraction, &nsites, &nconst, cused);
		sprintf(msg,"nsites: %d",nsites);
		LOG_MSG(msg);
		this->data = (FIELD_DATA *)malloc(nsites*sizeof(FIELD_DATA));
        get_fielddata(&axis, &fraction, &nsites, &nconst, this->data);
		slice_changed = false;
	}
    // Get picture limits, set size of square
    // Paint squares

    if (axis == X_AXIS) {           // Y-Z plane
        xindex = 1;
        yindex = 2;
    } else if (axis == Y_AXIS) {   // X-Z plane
        xindex = 0;
        yindex = 2;
    } else if (axis == Z_AXIS) {   // X-Y plane
        xindex = 0;
        yindex = 1;
    }
    c = NX/2;     // Centre is at (c,c)
    rmax = 0;
    for (i=0; i<nsites; i++) {
        ix = this->data[i].site[xindex];
        iy = this->data[i].site[yindex];
        rmax = MAX(rmax,ix - c);
        rmax = MAX(rmax,c+1 - ix);
        rmax = MAX(rmax,iy - c);
        rmax = MAX(rmax,c+1 - iy);
    }
    // rmax is the max number of site squares on any side of the centre (c,c) where c = NX/2
    // and each site square has width=1, and the site with ix extends from ix-1 to ix (etc.)
    cp = CANVAS_WIDTH/2;
    dp = int(0.95*CANVAS_WIDTH);
    w = int(dp/(2*rmax));  // width of mapped site square on canvas in pixels
    sprintf(msg,"NX,c,rmax,cp,dp,w: %d %d %d %d %d %d",NX,c,rmax,cp,dp,w);
    LOG_MSG(msg);
    // Blob slice is assumed to extend from c-rmax to c+rmax in both directions.
    // Note that a site square has a width = 1, i.e. extends (-0.5, 0.5) about the nominal position.
    // This must fit into a canvas of width and height = CANVAS_WIDTH
//    a = (dp - cp)/(2*rmax - c);
//    b = cp - a*c;
    // A site square at (ix,iy) maps to a canvas square of size (w,w) with:
    // (xp,yp)at (xp0 + (ix-1)*w, xp0 + (iy-1)*w) (presumably this is flipped in the y, i.e. about x axis)
    // where xp0 = cp - rmax*w
    xp0 = cp - rmax*w;
    d0 = w;
    double cmax = 0;
    for (i=0; i<nsites; i++) {
        cmax = MAX(cmax,data[i].conc[constituent]);
    }

    brush.setStyle(Qt::SolidPattern);
    for (i=0; i<nsites; i++) {
        ix = this->data[i].site[xindex];
        iy = this->data[i].site[yindex];
        xp = cp + (ix-c)*w;
        yp = cp + (iy-1-c)*w;
        double f = data[i].conc[constituent]/cmax;
        chooseColor(f,rgbcol);
        brush.setColor(QColor(rgbcol[0],rgbcol[1],rgbcol[2]));
        sprintf(msg,"i,rgbcol,f: %d %d %d %d %f",i,rgbcol[0],rgbcol[1],rgbcol[2],f);
        LOG_MSG(msg);
        scene->addRect(xp,yp,w,w,Qt::NoPen, brush);
    }
    for (i=0; i<nsites; i++) {
        ix = this->data[i].site[xindex];
        iy = this->data[i].site[yindex];
        xp = cp + (ix-c)*w;
        yp = cp + (iy-1-c)*w;
        volume = this->data[i].volume;
        if (volume > 0) {
            scale = pow(volume,0.3333);
            d = scale*d0;
            brush.setColor(QColor(0,255,0));
            scene->addEllipse(xp+(w-d)/2,yp+(w-d)/2,d,d,Qt::NoPen, brush);
        }
    }

//    brush.setColor(QColor(150,100,0));
//    brush.setStyle(Qt::SolidPattern);
//    scene->addEllipse(10,10,20,20,Qt::NoPen, brush);
//    text = scene->addText("FDC");
//    text->setPos(35, 10);

//    brush.setColor(QColor(200,60,100));
//    brush.setStyle(Qt::SolidPattern);
//    scene->addEllipse(10,40,20,20,Qt::NoPen, brush);
//    text = scene->addText("MRC");
//    text->setPos(35, 40);

//    brush.setColor(QColor(30,20,255));
//    scene->addEllipse(10,70,20,20,Qt::NoPen, brush);
//    text = scene->addText("Naive B cell");
//    text->setPos(35, 70);

//    brush.setColor(QColor(0,200,255));
//    scene->addEllipse(10,100,20,20,Qt::NoPen, brush);
//    text = scene->addText("CCR7 UP");
//    text->setPos(35, 100);

//    brush.setColor(QColor(50,255,150));
//    scene->addEllipse(10,130,20,20,Qt::NoPen, brush);
//    text = scene->addText("EBI2 UP");
//    text->setPos(35, 130);

//    brush.setColor(QColor(255,255,0));
//    scene->addEllipse(10,160,20,20,Qt::NoPen, brush);
//    text = scene->addText("BCL6 HI");
//    text->setPos(35, 160);

//    brush.setColor(QColor(0,150,0));
//    scene->addEllipse(10,190,20,20,Qt::NoPen, brush);
//    text = scene->addText("BCL6 LO");
//    text->setPos(35, 190);

//    brush.setColor(QColor(128,128,128));
//    scene->addEllipse(10,220,20,20,Qt::NoPen, brush);
//    text = scene->addText("Max divisions");
//    text->setPos(35, 220);

//    brush.setColor(QColor(255,0,0));
//    scene->addEllipse(10,250,20,20,Qt::NoPen, brush);
//    text = scene->addText("Plasma cell");
//    text->setPos(35, 250);

//    brush.setColor(QColor(255,130,0));
//    scene->addEllipse(10,280,20,20,Qt::NoPen, brush);
//    text = scene->addText("CD4 T cell");
//    text->setPos(35, 280);

    QGraphicsView* view = new QGraphicsView(field_page);
    view->setScene(scene);
    view->setGeometry(QRect(0, 0, 700, 700));
    view->show();
}

void Field::chooseColor(double f, int rgbcol[])
{
    if (cused[constituent] == 1) {
        rgbcol[0] = int((1-f)*255);
        rgbcol[1] = 0;
        rgbcol[2] = int(f*255);
    } else {
        rgbcol[0] = rgbcol[1] = rgbcol[2] = 255;
    }
}

void Field::makeConcPlot(QMdiArea *mdiArea)
{
    QString tag = "conc";
    QString title = "Concentration";
    pG = new Plot(tag,tag);
//    pG->setGeometry(QRect(100, 100, 300, 200));
    pG->setTitle(title);
//    pG->setAxisTitle(QwtPlot::yLeft, title);

    mdiArea->addSubWindow(pG);
    pG->show();
}

void Field::updateConcPlot()
{
    int nc, nmu, i, ichemo;
    double dx, x[1000], y[1000], conc[4000], cmax;

    get_concdata(&nc, &dx, conc);
    nmu = int(nc*dx*1.0e4);
    sprintf(msg,"updateConcPlot: %d %f %d",nc,dx,nmu);
    LOG_MSG(msg);
    pG->setAxisScale(QwtPlot::xBottom, 0, nmu, 0);
    QPen *pen = new QPen();
    QColor pencolor[] = {Qt::black, Qt::red, Qt::blue, Qt::darkGreen, Qt::magenta, Qt::darkCyan };
    pen->setColor(pencolor[0]);
    pG->curve[0]->setPen(*pen);
    ichemo = constituent;
    cmax = 0;
    for (i=0; i<nc; i++) {
        x[i] = i*dx*1.0e4;
        y[i] = conc[i*MAX_CHEMO+ichemo];
        sprintf(msg,"%d %f %f",i,x[i],y[i]);
        LOG_MSG(msg);
        cmax = MAX(cmax,y[i]);
    }
    pG->setAxisScale(QwtPlot::yLeft, 0, cmax, 0);
    pG->curve[0]->setData(x, y, nc);
    pG->replot();
    delete pen;
}
