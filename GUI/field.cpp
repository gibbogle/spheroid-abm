#include "field.h"
#include "log.h"

#include "global.h"

LOG_USE();

//double concData[4000];
//int conc_nc;
//double conc_dx;
//int MAX_CHEMO;
//int NX, NY, NZ;
//double dfraction;
//double volProb[100];
//int vol_nv;
//double vol_v0;
//double vol_dv;
//double oxyProb[100];
//int oxy_nv;
//double oxy_dv;

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
//Field::Field(QWidget *page2D)
Field::Field(QWidget *aParent) : QWidget(aParent)
{
//	field_page = page2D;
    field_page = aParent;
    axis = Z_AXIS;
    fraction = 0;
//    const_name[CFSE] = "CFSE";
//    const_name[OXYGEN] = "Oxygen";
//    const_name[GLUCOSE] = "Glucose";
//    const_name[TRACER] = "Tracer";
//    const_name[DRUG_A] = "Drug A";
//    const_name[DRUG_A_METAB_1] = "Drug A metabolite 1";
//    const_name[DRUG_A_METAB_2] = "Drug A metabolite 2";
//    const_name[DRUG_B] = "Drug B";
//    const_name[DRUG_B_METAB_1] = "Drug B metabolite 1";
//    const_name[DRUG_B_METAB_2] = "Drug B metabolite 2";
//    const_name[GROWTH_RATE] = "Growth rate";
    constituent = OXYGEN;
    slice_changed = true;
    setConcPlot(true);
    setVolPlot(true);
    setOxyPlot(true);
    pGconc = NULL;
    pGvol = NULL;
    pGoxy = NULL;
    ifield = 0;
    view = new MyQGraphicsView(field_page);
    constituent_rb_list = NULL;
    vbox_constituent = NULL;
    buttonGroup_constituent = new QButtonGroup;
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
// To create the group of radiobuttons for constituent selection.
// This uses information about active constituents fetched from the DLL.
//------------------------------------------------------------------------------------------------
void Field::setConstituentButtons(QGroupBox *gbox, QButtonGroup *bg, QVBoxLayout **vbox, QRadioButton ***rb_list, QString tag)
{
    int ivar;
    QString name, str;
    int **ip;
    QRadioButton **p;
    QRadioButton *rb;

    p = *rb_list;
    LOG_QMSG("setConstituentButtons: " + tag);
    if (p) {
        LOG_MSG("rb_list not NULL, delete it");
        for (ivar=0; ivar<Global::nvars_used; ivar++) {
            rb = p[ivar];
            bg->removeButton(rb);
            delete rb;
        }
        delete p;
    }
    if (!*vbox) {
        LOG_MSG("vbox = NULL, create it");
        *vbox = new QVBoxLayout;
        gbox->setLayout(*vbox);
    }
    name = "rb_constituent_"+tag;
    LOG_QMSG(name);
    *rb_list = new QRadioButton*[Global::nvars_used];
    p = *rb_list;
//    sprintf(msg,"rb_list: %p vbox: %p bg: %p nvars_used: %d",p,*vbox,bg,Global::nvars_used);
//    LOG_MSG(msg);
    for (ivar=0; ivar<Global::nvars_used; ivar++) {
        str = Global::var_string[ivar];
        p[ivar] = new QRadioButton;
        p[ivar]->setText(str);
        p[ivar]->setObjectName(name+ivar);
        (*vbox)->addWidget(p[ivar]);
        p[ivar]->setEnabled(true);
        bg->addButton(p[ivar],ivar);
    }
    p[1]->setChecked(true);   // Oxygen
    QRect rect = gbox->geometry();
    rect.setHeight(25*Global::nvars_used);
    gbox->setGeometry(rect);
    gbox->show();
}

//------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------
void Field::selectConstituent()
{
    int iconst, res;
    QStringList items;

    LOG_MSG("selectConstituent");
    /*
    get_fieldinfo(&NX, &axis, &fraction, &nsites, &nconst, const_used, &res);
    if (nconst != MAX_CONC) {
        sprintf(msg,"Error: selectConstituent: get_fieldinfo: MAX_CONC != MAX_CHEMO: %d %d",MAX_CONC,Global::MAX_CHEMO);
        LOG_MSG(msg);
        exit(1);
    }
    for (iconst=0; iconst<MAX_CONC+2; iconst++) {
        if (iconst == constituent) continue;
        if (const_used[iconst] == 1) {
            items << const_name[iconst];
        }
    }
    */
    for (iconst=0; iconst<Global::nvars_used; iconst++) {
        if (iconst == constituent) continue;
        items << Global::var_string[iconst];
    }
    bool ok;
    QString item = QInputDialog::getItem(this, tr("QInputDialog::getItem()"),
                                          tr("Constituent:"), items, 0, false, &ok);
    if (ok && !item.isEmpty()) {
//        for (iconst=0; iconst<MAX_CONC+2; iconst++) {
//            if (item == const_name[iconst]) {
        for (iconst=0; iconst<Global::nvars_used; iconst++) {
            if (item == Global::var_string[iconst]) {
                constituent = iconst;
                if (useConcPlot)
                    updateConcPlot();
                break;
            }
        }
    }
}

//------------------------------------------------------------------------------------------------
// This is to pass the info from MainWindow to Field
// The constituent with index ic in the list of used constituents has index cvar_index[ic] in the DLL.
//------------------------------------------------------------------------------------------------
void Field::setConstUsage(int nv, int *cv_index)
{
    nvars_used = nv;
    for (int i=0; i<nv; i++) {
        cvar_index[i] = cv_index[i];
    }
}


//------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------
void Field::setConstituent(QAbstractButton *button)
{
    int res;
//    QMessageBox msgBox;
//    msgBox.setText("setConstituent");
//    msgBox.exec();
    LOG_MSG("setConstituent");
	int prev_constituent = constituent;
    QString text = button->text();
    for (int ivar=0; ivar<Global::nvars_used; ivar++) {
        if (text == Global::var_string[ivar]) {
            constituent = ivar;
        }
    }
/*
//    constituentText = text;
    if (text.compare("CFSE") == 0)
        constituent = CFSE;
    else if (text.compare("Oxygen") == 0)
        constituent = OXYGEN;
    else if (text.compare("Glucose") == 0)
        constituent = GLUCOSE;
    else if (text.compare("Drug A") == 0)
        constituent = DRUG_A;
    else if (text.compare("Drug A metabolite 1") == 0)
        constituent = DRUG_A_METAB_1;
    else if (text.compare("Drug A metabolite 2") == 0)
        constituent = DRUG_A_METAB_2;
    else if (text.compare("Drug B") == 0)
        constituent = DRUG_B;
    else if (text.compare("Drug B metabolite 1") == 0)
        constituent = DRUG_B_METAB_1;
    else if (text.compare("Drug B metabolite 2") == 0)
        constituent = DRUG_B_METAB_2;
    else if (text.compare("Growth rate") == 0)
        constituent = GROWTH_RATE;
 */
    if (constituent != prev_constituent) {
		constituent_changed = true;
        LOG_MSG("setConstituent");
        displayField(hour,&res);
	}
}

//------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------
void Field::setPlane(QAbstractButton *button)
{
    int res;
//    QMessageBox msgBox;
//    msgBox.setText("setPlane");
//    msgBox.exec();
    LOG_MSG("setPlane");
    QString text = button->text();
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
        displayField(hour,&res);
	}
}

//------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------
void Field::setFraction(QString text)
{
    int res;
	double prev_fraction = fraction;
	fraction = text.toDouble();
	if (fraction != prev_fraction) {
        displayField(hour,&res);
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
void Field::setUseLogScale(bool use_logscale)
{
    use_log = use_logscale;
}

//-----------------------------------------------------------------------------------------
// New version, site/cell size is fixed, the blob grows
//-----------------------------------------------------------------------------------------
void Field::displayField(int hr, int *res)
{
    QGraphicsScene* scene = new QGraphicsScene(QRect(0, 0, CANVAS_WIDTH, CANVAS_WIDTH));
    QBrush brush;
    int i, xindex, yindex, ix, iy, w, rgbcol[3];
    double xp, yp, d0, d, volume, scale, cmin, cmax, rmax;
    double a, b, Wc;
    int Nc;

    LOG_MSG("displayField");
    use_log = false;    // temporary
    *res = 0;
    hour = hr;
	if (slice_changed) {
        LOG_MSG("get_fieldinfo:");
        get_fieldinfo(&NX, &axis, &fraction, &nsites, &nconst, const_used, res);
        if (*res != 0) {
            printf("Error: get_fieldinfo: FAILED\n");
            LOG_MSG("Error: get_fieldinfo: FAILED");
            return;
        }
        if (nconst != MAX_CONC) {
            sprintf(msg,"Error: displayField: get_fieldinfo: MAX_CONC != MAX_CHEMO: %d %d",MAX_CONC,Global::MAX_CHEMO);
            LOG_MSG(msg);
            exit(1);
        }
        this->data = (FIELD_DATA *)malloc(nsites*sizeof(FIELD_DATA));
        LOG_MSG("get_fielddata:");
        get_fielddata(&axis, &fraction, &nsites, &nconst, this->data, res);
        if (*res != 0) {
            LOG_MSG("Error: get_fielddata: FAILED");
            return;
        }
        LOG_MSG("get_fielddata: OK");
        slice_changed = false;
    }

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
 NX = size of lattice
 Nc = # of sites to fill the canvas from side to side (or top to bottom) = (2/3)NX
 Wc = canvas width (pixels)
 w = site width = Wc/Nc
 xp = a.ix + b
 yp = a.iy + b
 blob centre at (NX/2,NX/2) maps to canvas centre at (Wc/2,Wc/2)
 => Wc/2 = a.NX/2 + b
 The width of Nc sites maps to the canvas width
 => Wc = a.Nc
 => a = Wc/Nc, b = Wc/2 - a.NX/2
*/
    Nc = (2*NX)/3;
    Wc = CANVAS_WIDTH;
    w = Wc/Nc;
    a = w;
    b = Wc/2 - a*NX/2;
    d0 = w*Global::dfraction;
    cmin = 1.0e10;
    cmax = 0;
    rmax = 0;
    for (i=0; i<nsites; i++) {
//        rmax = MAX(rmax,data[i].dVdt);    // use conc[MAX_CHEMO+1] for dVdt
        rmax = MAX(rmax,data[i].conc[Global::MAX_CHEMO+1]);
        cmin = MIN(MAX(cmin,0),data[i].conc[constituent]);
        cmax = MAX(cmax,data[i].conc[constituent]);
        // Flip it
//        data[i].site[yindex] = NX - data[i].site[yindex];
    }
    brush.setStyle(Qt::SolidPattern);
    brush.setColor(QColor(0,0,0));
    scene->addRect(0,0,CANVAS_WIDTH,CANVAS_WIDTH,Qt::NoPen, brush);
    view->setScene(scene);
    view->setGeometry(QRect(0, 0, 700, 700));
    if (cmax == 0) {
        view->show();
        return;
    }
    for (i=0; i<nsites; i++) {
        ix = this->data[i].site[xindex];
        iy = this->data[i].site[yindex];
        xp = int(a*ix + b - w);
        yp = int(a*iy + b - w);
        chooseFieldColor(data[i].conc[constituent],cmin,cmax,use_log,rgbcol);
//        sprintf(msg,"c: %f %f %f rgbcol: %d %d %d",data[i].conc[constituent],cmin,cmax,rgbcol[0],rgbcol[1],rgbcol[2]);
//        LOG_MSG(msg);
        brush.setColor(QColor(rgbcol[0],rgbcol[1],rgbcol[2]));
        scene->addRect(xp,yp,w,w,Qt::NoPen, brush);
    }
    for (i=0; i<nsites; i++) {
        ix = this->data[i].site[xindex];
        iy = this->data[i].site[yindex];
        xp = int(a*ix + b - w);
        yp = int(a*iy + b - w);
        volume = this->data[i].volume;      // = 0 if there is no cell
        if (volume > 0) {
            scale = pow(volume,0.3333);
            d = scale*d0;   // fix this - need to change d0
//            double f = data[i].dVdt/rmax;
            double f = data[i].conc[GROWTH_RATE]/rmax;
            chooseRateColor(f,rgbcol);
            brush.setColor(QColor(rgbcol[0],rgbcol[1],rgbcol[2]));
            scene->addEllipse(xp+(w-d)/2,yp+(w-d)/2,d,d,Qt::NoPen, brush);
        }
    }
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
//-----------------------------------------------------------------------------------------
void Field::chooseFieldColor(double c, double cmin, double cmax, bool use_logscale, int rgbcol[])
{
    double f, denom, logcmin, logc;
    int rgb_lo[3], rgb_hi[3], i;

    if (use_logscale) {
        if (cmin == cmax) {
            f = 1;
        } else {
//            if (cmin > 0.0001)
//                logcmin = log(cmin);
//            else
//                logcmin = 1.0e6;
//            if (c > 0.0001)
//                logc = log(c);
//            else
//                logc = 1.0e6;
            cmin = max(cmin, 0.00001);
            c = max(c, 0.00001);
            denom = (log(cmax) - log(cmin));
            if (denom < 0.001)
                f = 1;
            else
                f = (log(c) - log(cmin))/denom;
        }
    } else {
        f = c/cmax;
    }
    if (constituent == OXYGEN) {
        rgb_hi[0] =   0; rgb_hi[1] =   0; rgb_hi[2] = 0;
        rgb_lo[0] = 255; rgb_lo[1] =   0; rgb_lo[2] = 0;
        for (i=0; i<3; i++) {
            rgbcol[i] = int((1-f)*rgb_lo[i] + f*rgb_hi[i]);
            if (rgbcol[i] < 0 || rgbcol[i] > 255) {
                sprintf(msg,"chooseFieldColor: %f %f %f %f %d %d",c,cmin,cmax,f,i,rgbcol[i]);
                LOG_MSG(msg);
                exit(1);
            }
        }
    } else {
        rgb_hi[0] =   0; rgb_hi[1] =   255; rgb_hi[2] = 255;
        rgb_lo[0] =   0; rgb_lo[1] =   0; rgb_lo[2] = 0;
        for (i=0; i<3; i++) {
            rgbcol[i] = int((1-f)*rgb_lo[i] + f*rgb_hi[i]);
            if (rgbcol[i] < 0 || rgbcol[i] > 255) {
                sprintf(msg,"chooseFieldColor: %f %f %f %f %d %d",c,cmin,cmax,f,i,rgbcol[i]);
                LOG_MSG(msg);
                exit(1);
            }
        }
    }
}

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
void Field::chooseRateColor(double f, int rgbcol[])
{
    int rgb_lo[3], rgb_hi[3], i;

    rgb_hi[0] = 0; rgb_hi[1] = 255; rgb_hi[2] = 0;
    rgb_lo[0] = 0; rgb_lo[1] =  64; rgb_lo[2] = 0;
    for (i=0; i<3; i++) {
        rgbcol[i] = int((1-f)*rgb_lo[i] + f*rgb_hi[i]);
    }
}

//-----------------------------------------------------------------------------------------
// Now 'constituent' is an index of the active constituents: 0 - nvars_used
//-----------------------------------------------------------------------------------------
void Field::getTitle(int iconst, QString *title)
{
//    if (constituent < Global::MAX_CHEMO+1) {
//        *title = const_name[constituent] + " Concentration";
//    } else if (constituent == GROWTH_RATE) {
//        *title = "Growth rate";
//    }
    QString name = Global::var_string[iconst];
    if (Global::GUI_to_DLL_index[iconst] <= Global::MAX_CHEMO) {
        *title = name + " Concentration";
    } else {
        *title = name;
    }
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
    getTitle(constituent,&title);
    LOG_QMSG(title);
    pGconc->setTitle(title);

    mdiArea->addSubWindow(pGconc);
    pGconc->show();
}

//-----------------------------------------------------------------------------------------
// Note that growth_rate is now constituent MAX_CHEMO+1
//-----------------------------------------------------------------------------------------
void Field::updateConcPlot()
{
    int nc, nmu, i, ichemo;
    double dx, x[1000], y[1000], *conc, cmax, cmin;
    QString title = "Concentration";

//    LOG_MSG("UpdateConcPlot");
    dx = Global::conc_dx;
    nc = Global::conc_nc;
    conc = Global::concData;
    if (nc == 0) return;
    nmu = int(nc*dx*1.0e4);
    getTitle(constituent,&title);
    pGconc->setTitle(title);
    pGconc->setAxisScale(QwtPlot::xBottom, 0, nmu, 0);
    QPen *pen = new QPen();
    QColor pencolor[] = {Qt::black, Qt::red, Qt::blue, Qt::darkGreen, Qt::magenta, Qt::darkCyan };
    pen->setColor(pencolor[0]);
    pGconc->curve[0]->setPen(*pen);
//    ichemo = constituent;
    ichemo = Global::GUI_to_DLL_index[constituent];
    cmin = 1.0e10;
    cmax = 0;
    for (i=0; i<nc; i++) {
        x[i] = i*dx*1.0e4;
        y[i] = conc[i*(Global::MAX_CHEMO+2)+ichemo];
        cmin = MIN(cmin,y[i]);
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
    Global::vol_nv = 20;
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

//    LOG_MSG("UpdateVolPlot");
    prob = Global::volProb;
    v1 = Global::vol_v0 - Global::vol_dv/2;
    v2 = Global::vol_v0 + (Global::vol_nv-0.5)*Global::vol_dv;
    pGvol->setAxisScale(QwtPlot::xBottom, v1, v2, 0);
    QPen *pen = new QPen();
    QColor pencolor[] = {Qt::black, Qt::red, Qt::blue, Qt::darkGreen, Qt::magenta, Qt::darkCyan };
    pen->setColor(pencolor[0]);
    pGvol->curve[0]->setPen(*pen);

    pmax = 0;
    for (i=0; i<Global::vol_nv; i++) {
        x[i] = Global::vol_v0 + i*Global::vol_dv;
        y[i] = prob[i];
        pmax = MAX(pmax,y[i]);
    }
    pmax = 4.0/Global::vol_nv;  // try this
    pGvol->setAxisScale(QwtPlot::yLeft, 0, pmax, 0);
    pGvol->curve[0]->setData(x, y, Global::vol_nv);

    pGvol->replot();
    delete pen;
}

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
void Field::makeOxyPlot(QMdiArea *mdiArea)
{
    LOG_MSG("makeOxyPlot");
    Global::oxy_nv = 20;
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
    prob = Global::oxyProb;
    v1 = 0;
    v2 = Global::oxy_nv*Global::oxy_dv;
//    sprintf(msg,"updateOxyPlot: %d %f %f %f", oxy_nv, oxy_dv, v1, v2);
//    LOG_MSG(msg);
    pGoxy->setAxisScale(QwtPlot::xBottom, v1, v2, 0);
    QPen *pen = new QPen();
    QColor pencolor[] = {Qt::black, Qt::red, Qt::blue, Qt::darkGreen, Qt::magenta, Qt::darkCyan };
    pen->setColor(pencolor[0]);
    pGoxy->curve[0]->setPen(*pen);

    pmax = 0;
    for (i=0; i<Global::oxy_nv; i++) {
        x[i] = (i+0.5)*Global::oxy_dv;
        y[i] = prob[i];
        pmax = MAX(pmax,y[i]);
    }
    i = pmax/0.1;
    pmax = (i+1)*0.1;
    pGoxy->setAxisScale(QwtPlot::yLeft, 0, pmax, 0);
    pGoxy->curve[0]->setData(x, y, Global::oxy_nv);

    pGoxy->replot();
    delete pen;
}

