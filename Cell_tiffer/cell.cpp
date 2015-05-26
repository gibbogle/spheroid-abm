#include "mainwindow.h"
#include "ui_mainwindow.h"

#include <cstdio>
#include <vector>

#include <algorithm>
#include <math.h>
#include <string.h>
#include <string>
#include <sstream>
#include <assert.h>
#include <ctime>


#define STR_LEN 64
#define BIG 1.0e6


//-----------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------
int MainWindow::readCellData(const char *cellFile)
{
    char line[STR_LEN];
    int icell, ix, iy, iz, itype;
    double r;

    fprintf(fpout,"voxelsize: %f\n",voxelsize);
    FILE *fpcell = fopen(cellFile,"r");
    fgets(line, STR_LEN, fpcell);
    fprintf(fpout,"%s",line);
    sscanf(line,"%lf",&delta_x);
    fprintf(fpout,"delta_x: %lf\n",delta_x);
    fgets(line, STR_LEN, fpcell);
    fprintf(fpout,"%s",line);
    sscanf(line,"%d",&ncells);
    fprintf(fpout,"ncells: %d\n",ncells);
    location = (LOCATION *)malloc(ncells*sizeof(LOCATION));
    for (icell=0; icell<ncells; icell++) {
        if (fgets(line, STR_LEN, fpcell) == NULL) {
            printf("Error reading cell location file\n\n");
            fclose(fpcell);
            return 1;
        }
        fprintf(fpout,"%s",line);
        sscanf(line,"%d %d %d %d %lf",&itype,&ix,&iy,&iz,&r);
        fprintf(fpout,"itype,ix,iy,iz,r: %d %d %d %d %lf\n",itype,ix,iy,iz,r);
        location[icell].type = itype;
        location[icell].x = ix*delta_x;
        location[icell].y = iy*delta_x;
        location[icell].z = iz*delta_x;
        location[icell].r = r;
    }
    return 0;
}

//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------
int MainWindow::createTiffData()  //(LOCATION *locations, int ncells)
{
    int icell, itype, selectedtype;
    int nr, ix0, iy0, iz0, idx, idy, idz, id2;
    int iz1, iz2, ix, iy, iz, datadepth;
    double xmin, ymin, zmin, xmax, ymax, zmax, zmid, rng, rngz, radiusz, r2;
    LOCATION loc;

    if (ui->radioButton_type1->isChecked()) {
        selectedtype = 1;
    } else if (ui->radioButton_type2->isChecked()) {
        selectedtype = 2;
    } else {
        selectedtype = 3;
    }

    xmin = ymin = zmin = 1.0e10;
    xmax = ymax = zmax = 0;

    for (icell=0; icell<ncells; icell++) {
        loc = location[icell];
        xmin = MIN(xmin,loc.x-loc.r);
        ymin = MIN(ymin,loc.y-loc.r);
        zmin = MIN(zmin,loc.z-loc.r);
        xmax = MAX(xmax,loc.x+loc.r);
        ymax = MAX(ymax,loc.y+loc.r);
        zmax = MAX(zmax,loc.z+loc.r);
    }
    zmid = (zmin + zmax)/2;
    fprintf(fpout,"xmin,xmax: %f %f\n",xmin,xmax);
    fprintf(fpout,"ymin,ymax: %f %f\n",ymin,ymax);
    fprintf(fpout,"zmin,zmax: %f %f\n",zmin,zmax);
    margin = 10;

    rng = MAX(xmax-xmin,ymax-ymin);
    rngz = zmax-zmin;
    radiusz = rngz/2;
    rng += 2*margin;
    rngz += 2*margin;
    width = rng/voxelsize + 1;
    height = rng/voxelsize + 1;
    datadepth = rngz/voxelsize + 1;
    incell = (unsigned char *)malloc(width*height*datadepth*sizeof(unsigned char));
    memset(incell,0,width*height*datadepth);
    xysize = width*height;
    printf("width, height, datadepth: %d %d %d\n",width, height, datadepth);
    fprintf(fpout,"width, height, datadepth: %d %d %d\n",width, height, datadepth);
    fflush(fpout);
    depth = 1;
    p_im = (unsigned char *)malloc(width*height*depth*sizeof(unsigned char));
    memset(p_im,0,width*height*depth);

    for (icell=0; icell<ncells; icell++) {
        loc = location[icell];
        if (loc.type != selectedtype && selectedtype != 3) continue;
        r2 = loc.r*loc.r;
        ix0 = (int)((loc.x-xmin)/voxelsize) + margin;
        iy0 = (int)((loc.y-ymin)/voxelsize) + margin;
        iz0 = (int)((loc.z-zmin)/voxelsize) + margin;
        nr = (int)(loc.r/voxelsize);
        for (idx=-nr; idx<=nr; idx++) {
            for (idy=-nr; idy<=nr; idy++) {
                for (idz=-nr; idz<=nr; idz++) {
                    id2 = idx*idx + idy*idy + idz*idz;
                    if (id2*voxelsize*voxelsize < r2) {
                        Vin(ix0+idx,iy0+idy,iz0+idz) = 1;
                    }
                }
            }
        }
    }

    // zfraction = 0 corresponds to z = zmid = (zmin + zmax)/2
    // zfraction corresponds to z = zmid + zfraction*radiusz
    iz1 = (int)((zmid + zfraction*radiusz - zthickness/2 - zmin)/voxelsize) + margin;
    iz1 = MAX(iz1,0);
    iz1 = MIN(iz1,datadepth-1);
    iz2 = (int)((zmid + zfraction*radiusz + zthickness/2 - zmin)/voxelsize) + margin;
    iz2 = MAX(iz2,0);
    iz2 = MIN(iz2,datadepth-1);
    // Detect any incell(ix,iy,iz) for iz1 <= iz <= iz2
    for (ix=0; ix<width; ix++) {
        for (iy=0; iy<height; iy++) {
            for (iz=iz1; iz<=iz2; iz++) {
                if (Vin(ix,iy,iz) == 1) {
                    V(ix,iy,0) = 255;
                    break;
                }
            }
        }
    }

//    EDGE edge;
//    APOINT p;
    /*
    printf("ne: %d\n",net->ne);
    fprintf(fpout,"ne: %d\n",net->ne);
    // First determine the required buffer size to hold the voxels
    printf("voxelsize: %f %f %f\n",voxelsize[0],voxelsize[1],voxelsize[2]);
    wx = 0;
    wy = 0;
    wz = 0;
    for (ie=0; ie<net->ne; ie++) {
        edge = net->edgeList[ie];
        for (ip=0; ip<edge.npts; ip++) {
//            printf("%d %d %d %d\n",ie,edge.npts,ip,edge.pt[ip]);
//            fprintf(fpout,"%d %d %d %d\n",ie,edge.npts,ip,edge.pt[ip]);
//            fflush(fpout);
            p = net->point[edge.pt[ip]];
//            printf("%6.1f %6.1f %6.1f  %6.2f\n",p.x,p.y,p.z,p.d);
//            fprintf(fpout,"%6.1f %6.1f %6.1f  %6.2f\n",p.x,p.y,p.z,p.d);
//            fflush(fpout);
            wx = MAX(wx,(p.x + p.d/2.));
            wy = MAX(wy,(p.y + p.d/2.));
            wz = MAX(wz,(p.z + p.d/2.));
        }
    }
    printf("wx,wy,wz: %f %f %f\n",wx,wy,wz);
    width = (int)((wx+margin)/voxelsize[0]+10.);
    height = (int)((wy+margin)/voxelsize[1]+10.);
    depth = (int)((wz+margin)/voxelsize[2]+10.);
    xysize = width*height;
    printf("width, height, depth: %d %d %d\n",width, height, depth);
    fprintf(fpout,"width, height, depth: %d %d %d\n",width, height, depth);
    fflush(fpout);
    p_im = (unsigned char *)malloc(width*height*depth*sizeof(unsigned char));
    memset(p_im,0,width*height*depth);
    for (ie=0; ie<net->ne; ie++) {
        edge = net->edgeList[ie];
        for (ip=0; ip<edge.npts; ip++) {
//            fprintf(fpout,"ie, npts,ip: %d %d %d\n",ie,edge.npts,ip);
//            fflush(fpout);
            p = net->point[edge.pt[ip]];
            x0 = p.x/voxelsize[0];   // voxel nearest to the point
            y0 = p.y/voxelsize[1];
            z0 = p.z/voxelsize[2];
            r = p.d/2 + margin;
//            fprintf(fpout,"point: %d  %6.1f %6.1f %6.1f  %d %d %d  %6.2f\n",ip,p.x,p.y,p.z,x0,y0,z0,r);
//            fflush(fpout);
            r2 = r*r;
            nx = r/voxelsize[0] + 1;
            ny = r/voxelsize[1] + 1;
            nz = r/voxelsize[2] + 1;
//            printf("nx,ny,nz,x1,y1,z1: %d %d %d  %d %d %d\n",nx,ny,nz,x1,y1,z1);
            for (ix = -nx; ix <= nx; ix++) {
                dx = ix*voxelsize[0];
                x = x0+ix;
                if (x < 0 || x >= width) continue;
                for (iy = -ny; iy<=ny; iy++) {
                    dy = iy*voxelsize[1];
                    y = y0+iy;
                    if (y < 0 || y >= height) continue;
                    for (iz = -nz; iz<=nz; iz++) {
                        dz = iz*voxelsize[2];
                        z = z0+iz;
                        if (z < 0 || z >= depth) continue;
                        d2 = dx*dx+dy*dy+dz*dz;
                        if (d2 < r2) {
                            V(x0+ix,y0+iy,z0+iz) = 255;
                        }
                    }
                }
            }
        }
    }
    */
    return 0;
}

/*
bool MainWindow::inSphere(LOCATION p)
{
    float dx, dy, dz;

    dx = p.x - sphereCentre[0];
    dy = p.y - sphereCentre[1];
    dz = p.z - sphereCentre[2];
    if (dx*dx+dy*dy+dz*dz < sphereRadius*sphereRadius)
        return true;
    else
        return false;
}
*/

