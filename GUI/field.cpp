#include "field.h"
#include "log.h"

LOG_USE();

double concData[4000];
int conc_nc;
double conc_dx;
int MAX_CHEMO;
int NX, NY, NZ;
double dfraction;
double volProb[100];
int vol_nv;
double vol_v0;
double vol_dv;
double oxyProb[100];
int oxy_nv;
double oxy_dv;
bool goflag;

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
Field::Field(QWidget *page2D, bool save)
{
	field_page = page2D;
    save_images = save;
    axis = Z_AXIS;
    fraction = 0;
    const_name[OXYGEN] = "Oxygen";
    const_name[GLUCOSE] = "Glucose";
    const_name[DRUG_A] = "Drug A";
    const_name[DRUG_B] = "Drug B";
    const_name[DRUG_A_METAB] = "Drug A metabolite";
    const_name[DRUG_B_METAB] = "Drug B metabolite";
    const_name[GROWTH_RATE] = "Growth rate";
    constituent = OXYGEN;
    slice_changed = true;
    setConcPlot(true);
    setVolPlot(true);
    setOxyPlot(true);
    pGconc = NULL;
    pGvol = NULL;
    pGoxy = NULL;
    ifield = 0;
    view = new QGraphicsView(field_page);
}

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
Field::~Field()
{
}

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
bool Field::isConcPlot()
{
    return useConcPlot;
}

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
void Field::setConcPlot(bool status)
{
    useConcPlot = status;
}

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
bool Field::isVolPlot()
{
    return useVolPlot;
}

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
void Field::setVolPlot(bool status)
{
    useVolPlot = status;
}

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
bool Field::isOxyPlot()
{
    return useOxyPlot;
}

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
void Field::setOxyPlot(bool status)
{
    useOxyPlot = status;
}

//------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------
void Field::selectConstituent()
{
    int iconst;
    QStringList items;

    LOG_MSG("selectConstituent");
    get_fieldinfo(&NX, &axis, &fraction, &nsites, &nconst, const_used);
    for (iconst=0; iconst<MAX_CONC+1; iconst++) {
        if (iconst == constituent) continue;
        if (const_used[iconst] == 1) {
            items << const_name[iconst];
        }
    }

    bool ok;
    QString item = QInputDialog::getItem(this, tr("QInputDialog::getItem()"),
                                          tr("Constituent:"), items, 0, false, &ok);
    if (ok && !item.isEmpty()) {
        for (iconst=0; iconst<MAX_CONC+1; iconst++) {
            if (item == const_name[iconst]) {
                constituent = iconst;
                if (useConcPlot)
                    updateConcPlot();
                break;
            }
        }
    }
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
//    constituentText = text;
    if (text.compare("Oxygen") == 0)
        constituent = OXYGEN;
    else if (text.compare("Glucose") == 0)
        constituent = GLUCOSE;
    else if (text.compare("Drug A") == 0)
        constituent = DRUG_A;
    else if (text.compare("Drug A metabolite") == 0)
        constituent = DRUG_A_METAB;
    else if (text.compare("Drug B") == 0)
        constituent = DRUG_B;
    else if (text.compare("Drug B metabolite") == 0)
        constituent = DRUG_B_METAB;
    else if (text.compare("Growth rate") == 0)
        constituent = GROWTH_RATE;
    if (constituent != prev_constituent) {
		constituent_changed = true;
        LOG_MSG("setConstituent");
        displayField(hour);
	}
//    constituentText = const_name[constituent];
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
        LOG_MSG("setPlane");
        displayField(hour);
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
//		slice_changed = true;
        LOG_MSG("setFraction");
        displayField(hour);
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
void Field::setSaveImages(bool save)
{
    save_images = save;
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

//    bool ok;
    QString text;
    QStringList items;
//    items << tr("X-Y plane") << tr("X-Z plane") << tr("Y-Z plane");
//    items << QString("X-Y plane") << QString("X-Z plane") << QString("Y-Z plane");

//    QString item = QInputDialog::getItem(this,tr("QInputDialog::getItem()"),
//                                         tr("Slice plane:"), items, 0, false, &ok);
    QString item = "";
    if (!item.isEmpty()) {
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
// New version, site/cell size is fixed, the blob grows
//-----------------------------------------------------------------------------------------
void Field::displayField(int hr)
{
    QGraphicsScene* scene = new QGraphicsScene(QRect(0, 0, CANVAS_WIDTH, CANVAS_WIDTH));
    QBrush brush;
    int i, xindex, yindex, ix, iy, w, rgbcol[3];
    double xp, yp, d0, d, volume, scale;
    double a, b, Wc;
    int Nc = 50;
    bool growthRate;

    LOG_MSG("displayField");
    hour = hr;
	if (slice_changed) {
//        LOG_MSG("get_fieldinfo");
        get_fieldinfo(&NX, &axis, &fraction, &nsites, &nconst, const_used);
        if (nconst != MAX_CONC) {
            sprintf(msg,"Error: MAX_CONC != MAX_CHEMO in field.h");
            LOG_MSG(msg);
            exit(1);
        }
        this->data = (FIELD_DATA *)malloc(nsites*sizeof(FIELD_DATA));
//        LOG_MSG("get_fielddata");
        get_fielddata(&axis, &fraction, &nsites, &nconst, this->data);
		slice_changed = false;
//        LOG_MSG("got_fielddata");
    }
    goflag = true;
    if (constituent == GROWTH_RATE)
        growthRate = true;
    else
        growthRate = false;

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

/*
 Nc = # of sites to fill the canvas from side to side (or top to bottom)
 Wc = canvas width (pixels)
 w = site width = Wc/Nc
 xp = a.ix + b
 yp = a.iy + b
 blob centre at (Nx/2,Nx/2) maps to canvas centre at (Wc/2,Wc/2)
 => Wc/2 = a.Nx/2 + b
 The width of Nc sites maps to the canvas width
 => Wc = a.Nc
 => a = Wc/Nc, b = Wc/2 - a.Nx/2
*/
    Wc = CANVAS_WIDTH;
    w = Wc/Nc;
    a = w;
    b = Wc/2 - a*NX/2;
    d0 = w*dfraction;
    double cmax = 0;
    for (i=0; i<nsites; i++) {
        if (growthRate)
            cmax = MAX(cmax,data[i].dVdt);
        else
            cmax = MAX(cmax,data[i].conc[constituent]);
    }
//    LOG_MSG("got cmax");
    brush.setStyle(Qt::SolidPattern);
    brush.setColor(QColor(0,0,0));
    scene->addRect(0,0,CANVAS_WIDTH,CANVAS_WIDTH,Qt::NoPen, brush);
    for (i=0; i<nsites; i++) {
        ix = this->data[i].site[xindex];
        iy = this->data[i].site[yindex];
        xp = int(a*ix + b - w/2);
        yp = int(a*iy + b - w/2);
        if (growthRate) {
            brush.setColor(QColor(0,0,0));
        } else {
            double f = data[i].conc[constituent]/cmax;
            chooseColor(f,rgbcol);
            brush.setColor(QColor(rgbcol[0],rgbcol[1],rgbcol[2]));
//            sprintf(msg,"i,rgbcol,f: %d %d %d %d %f %f %f",i,rgbcol[0],rgbcol[1],rgbcol[2],data[i].conc[constituent],cmax,f);
//            LOG_MSG(msg);
        }
        scene->addRect(xp,yp,w,w,Qt::NoPen, brush);
        volume = this->data[i].volume;      // = 0 if there is no cell
        if (volume > 0) {
            scale = pow(volume,0.3333);
            d = scale*d0;   // fix this - need to change d0
            if (growthRate) {
//                double f = data[i].conc[constituent]/cmax;
                double f = data[i].dVdt/cmax;
                chooseRateColor(f,rgbcol);
                brush.setColor(QColor(rgbcol[0],rgbcol[1],rgbcol[2]));
//                sprintf(msg,"i,rgbcol,f: %d %d %d %d %f %f %f",i,rgbcol[0],rgbcol[1],rgbcol[2],data[i].dVdt,cmax,f);
//                LOG_MSG(msg);
            } else {
                brush.setColor(QColor(0,255,0));
            }
            scene->addEllipse(xp+(w-d)/2,yp+(w-d)/2,d,d,Qt::NoPen, brush);
        }
    }
    view->setScene(scene);
    view->setGeometry(QRect(0, 0, 700, 700));
//    view->setGeometry(QRect(0, 0, CANVAS_WIDTH+4, CANVAS_WIDTH+4));
    view->show();
    if (save_images) {
        scene->clearSelection();                                                  // Selections would also render to the file
        scene->setSceneRect(scene->itemsBoundingRect());                          // Re-shrink the scene to it's bounding contents
        QImage image(scene->sceneRect().size().toSize(), QImage::Format_ARGB32);  // Create the image with the exact size of the shrunk scene
        image.fill(Qt::transparent);                                              // Start all pixels transparent

        QPainter painter(&image);
        scene->render(&painter);
        ifield++;
        char filename[] = "image/field0000.png";
        char numstr[5];
        sprintf(numstr,"%04d",hour);
        for (int i=0; i<4; i++)
            filename[11+i] = numstr[i];
        image.save(filename);
    }
}

//-----------------------------------------------------------------------------------------
// Old version, blob scaled to fit the window
//-----------------------------------------------------------------------------------------
void Field::displayField1()
{
//    QGraphicsScene* scene = new QGraphicsScene(QRect(0, 0, 130, 280));
    QGraphicsScene* scene = new QGraphicsScene(QRect(0, 0, 690, 690));
    QBrush brush;
//    QGraphicsTextItem *text;
    int i, xindex, yindex, ix, iy, cp, dp, c, rmax, w, xp0, rgbcol[3];
    double xp, yp, d0, d, volume, scale;
//    double a, b;
    bool growthRate;

    LOG_MSG("displayField");
    if (slice_changed) {
        get_fieldinfo(&NX, &axis, &fraction, &nsites, &nconst, const_used);
        sprintf(msg,"nsites: %d",nsites);
        LOG_MSG(msg);
        if (nconst != MAX_CONC) {
            sprintf(msg,"Error: MAX_CONC != MAX_CHEMO in field.h");
            LOG_MSG(msg);
            exit(1);
        }
        this->data = (FIELD_DATA *)malloc(nsites*sizeof(FIELD_DATA));
        get_fielddata(&axis, &fraction, &nsites, &nconst, this->data);
        slice_changed = false;
    }
    LOG_MSG("got field data");
    if (constituent == GROWTH_RATE)
        growthRate = true;
    else
        growthRate = false;

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
//        sprintf(msg,"%d %d %d %d %d %d %d",i,ix,iy,ix-c,c+1-ix,iy-c,c+1-iy);
//        LOG_MSG(msg);
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
    sprintf(msg,"constituent, NX,c,rmax,cp,dp,w: %d %d %d %d %d %d %d",constituent,NX,c,rmax,cp,dp,w);
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
        if (growthRate)
            cmax = MAX(cmax,data[i].dVdt);
        else
            cmax = MAX(cmax,data[i].conc[constituent]);
    }

    brush.setStyle(Qt::SolidPattern);
    for (i=0; i<nsites; i++) {
        ix = this->data[i].site[xindex];
        iy = this->data[i].site[yindex];
        xp = cp + (ix-c)*w;
        yp = cp + (iy-1-c)*w;
        if (growthRate) {
            brush.setColor(QColor(0,0,0));
        } else {
            double f = data[i].conc[constituent]/cmax;
            chooseColor(f,rgbcol);
            brush.setColor(QColor(rgbcol[0],rgbcol[1],rgbcol[2]));
            sprintf(msg,"i,rgbcol,f: %d %d %d %d %f %f %f",i,rgbcol[0],rgbcol[1],rgbcol[2],data[i].conc[constituent],cmax,f);
            LOG_MSG(msg);
        }
        scene->addRect(xp,yp,w,w,Qt::NoPen, brush);
    }
    for (i=0; i<nsites; i++) {
        ix = this->data[i].site[xindex];
        iy = this->data[i].site[yindex];
        xp = cp + (ix-c)*w;
        yp = cp + (iy-1-c)*w;
        volume = this->data[i].volume;      // = 0 if there is no cell
        if (volume > 0) {
            scale = pow(volume,0.3333);
            d = scale*d0;
            if (growthRate) {
//                double f = data[i].conc[constituent]/cmax;
                double f = data[i].dVdt/cmax;
                chooseRateColor(f,rgbcol);
                brush.setColor(QColor(rgbcol[0],rgbcol[1],rgbcol[2]));
                sprintf(msg,"i,rgbcol,f: %d %d %d %d %f %f %f",i,rgbcol[0],rgbcol[1],rgbcol[2],data[i].dVdt,cmax,f);
                LOG_MSG(msg);
            } else {
                brush.setColor(QColor(0,255,0));
            }
            scene->addEllipse(xp+(w-d)/2,yp+(w-d)/2,d,d,Qt::NoPen, brush);
        }
    }
    LOG_MSG("make view");
    QGraphicsView* view = new QGraphicsView(field_page);
    view->setScene(scene);
    view->setGeometry(QRect(0, 0, 700, 700));
    view->show();
    LOG_MSG("showed view");

    scene->clearSelection();                                                  // Selections would also render to the file
    scene->setSceneRect(scene->itemsBoundingRect());                          // Re-shrink the scene to it's bounding contents
    QImage image(scene->sceneRect().size().toSize(), QImage::Format_ARGB32);  // Create the image with the exact size of the shrunk scene
    image.fill(Qt::transparent);                                              // Start all pixels transparent

    LOG_MSG("painter");
    QPainter painter(&image);
    LOG_MSG("render scene");
    scene->render(&painter);
    ifield++;
    char filename[] = "field0000.png";
    char numstr[5];
    sprintf(numstr,"%04d",ifield);
    for (int i=0; i<4; i++)
        filename[5+i] = numstr[i];
    image.save(filename);

}
//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
void Field::chooseColor(double f, int rgbcol[])
{
    if (const_used[constituent] == 1) {
        rgbcol[2] = int((1-f)*255);
        rgbcol[1] = 0;
        rgbcol[0] = int(f*255);
    } else {
        rgbcol[0] = rgbcol[1] = rgbcol[2] = 255;
    }
}

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
void Field::chooseRateColor(double f, int rgbcol[])
{
    if (const_used[constituent] == 1) {
        rgbcol[2] = int((1-f)*255);
        rgbcol[1] = 0;
        rgbcol[0] = int(f*255);
    } else {
        rgbcol[0] = rgbcol[1] = rgbcol[2] = 255;
    }
}

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
void Field::getTitle(QString *title)
{
    *title = const_name[constituent] + " Concentration";
}

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
void Field::makeConcPlot(QMdiArea *mdiArea)
{
    QString tag = "conc";
    QString title;
    if (pGconc != NULL) {
        LOG_MSG("pGconc not NULL");
        delete pGconc;
    }
    pGconc = new Plot(tag,tag);
    getTitle(&title);
    LOG_QMSG(title);
    pGconc->setTitle(title);

    mdiArea->addSubWindow(pGconc);
    pGconc->show();
}

//-----------------------------------------------------------------------------------------
// Note that growth_rate is constituent MAX_CHEMO
//-----------------------------------------------------------------------------------------
void Field::updateConcPlot()
{
    int nc, nmu, i, ichemo;
    double dx, x[1000], y[1000], *conc, cmax;
    QString title = "Concentration";

    LOG_MSG("UpdateConcPlot");
//    get_concdata(&nc, &dx, conc);
    dx = conc_dx;
    nc = conc_nc;
    conc = concData;
    if (nc == 0) return;
    nmu = int(nc*dx*1.0e4);
    sprintf(msg,"updateConcPlot: %d %f %d %d",nc,dx,nmu,MAX_CHEMO);
    LOG_MSG(msg);
    if (constituent < MAX_CHEMO) {
        getTitle(&title);
    } else {
        title = "Growth rate";
    }
    pGconc->setTitle(title);
    pGconc->setAxisScale(QwtPlot::xBottom, 0, nmu, 0);
    QPen *pen = new QPen();
    QColor pencolor[] = {Qt::black, Qt::red, Qt::blue, Qt::darkGreen, Qt::magenta, Qt::darkCyan };
    pen->setColor(pencolor[0]);
    pGconc->curve[0]->setPen(*pen);
    ichemo = constituent;
    cmax = 0;
    for (i=0; i<nc; i++) {
        x[i] = i*dx*1.0e4;
        y[i] = conc[i*(MAX_CHEMO+1)+ichemo];
        cmax = MAX(cmax,y[i]);
    }
    pGconc->setAxisScale(QwtPlot::yLeft, 0, cmax, 0);
    pGconc->curve[0]->setData(x, y, nc);
    pGconc->replot();
    delete pen;
}

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
void Field::makeVolPlot(QMdiArea *mdiArea)
{
    vol_nv = 20;
    QString tag = "vol";
    QString title = "Volume Distribution";
    if (pGvol != NULL) {
        LOG_MSG("pGvol not NULL");
        delete pGvol;
    }
    pGvol = new Plot(tag,tag);
    pGvol->setTitle(title);

    mdiArea->addSubWindow(pGvol);
    pGvol->show();
}

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
void Field::updateVolPlot()
{
    int i;
    double x[100], y[100], *prob, pmax, v1, v2;

    LOG_MSG("UpdateVolPlot");
    prob = volProb;
    v1 = vol_v0 - vol_dv/2;
    v2 = vol_v0 + (vol_nv-0.5)*vol_dv;
    sprintf(msg,"updateVolPlot: %d %f %f %f %f",vol_nv,vol_dv,vol_v0, v1, v2);
    LOG_MSG(msg);
    pGvol->setAxisScale(QwtPlot::xBottom, v1, v2, 0);
    QPen *pen = new QPen();
    QColor pencolor[] = {Qt::black, Qt::red, Qt::blue, Qt::darkGreen, Qt::magenta, Qt::darkCyan };
    pen->setColor(pencolor[0]);
    pGvol->curve[0]->setPen(*pen);

    pmax = 0;
    for (i=0; i<vol_nv; i++) {
        x[i] = vol_v0 + i*vol_dv;
        y[i] = prob[i];
//        sprintf(msg,"%d %f %f",i,x[i],y[i]);
//        LOG_MSG(msg);
        pmax = MAX(pmax,y[i]);
    }
    pmax = 4.0/vol_nv;  // try this
    pGvol->setAxisScale(QwtPlot::yLeft, 0, pmax, 0);
    pGvol->curve[0]->setData(x, y, vol_nv);

    pGvol->replot();
    delete pen;
}

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
void Field::makeOxyPlot(QMdiArea *mdiArea)
{
    LOG_MSG("makeOxyPlot");
    oxy_nv = 20;
    QString tag = "oxy";
    QString title = "Cell O2 Distribution";
    if (pGoxy != NULL) {
        LOG_MSG("pGoxy not NULL");
        delete pGoxy;
    }
    pGoxy = new Plot(tag,tag);
    pGoxy->setTitle(title);

    mdiArea->addSubWindow(pGoxy);
    pGoxy->show();
}

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
void Field::updateOxyPlot()
{
    int i;
    double x[100], y[100], *prob, pmax, v1, v2;

//    LOG_MSG("UpdateOxyPlot");
    prob = oxyProb;
    v1 = 0;
    v2 = oxy_nv*oxy_dv;
//    sprintf(msg,"updateOxyPlot: %d %f %f %f", oxy_nv, oxy_dv, v1, v2);
//    LOG_MSG(msg);
    pGoxy->setAxisScale(QwtPlot::xBottom, v1, v2, 0);
    QPen *pen = new QPen();
    QColor pencolor[] = {Qt::black, Qt::red, Qt::blue, Qt::darkGreen, Qt::magenta, Qt::darkCyan };
    pen->setColor(pencolor[0]);
    pGoxy->curve[0]->setPen(*pen);

    pmax = 0;
    for (i=0; i<oxy_nv; i++) {
        x[i] = (i+0.5)*oxy_dv;
        y[i] = prob[i];
//        sprintf(msg,"%d %f %f",i,x[i],y[i]);
//        LOG_MSG(msg);
        pmax = MAX(pmax,y[i]);
    }
//    pmax = 6.0/oxy_nv;  // try this
    i = pmax/0.1;
    pmax = (i+1)*0.1;
    pGoxy->setAxisScale(QwtPlot::yLeft, 0, pmax, 0);
    pGoxy->curve[0]->setData(x, y, oxy_nv);

    pGoxy->replot();
    delete pen;
}

