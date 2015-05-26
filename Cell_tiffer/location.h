#ifndef LOCATION_H
#define LOCATION_H

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

struct location_str
{
    int type;
    double x,y,z;
    double r;
};
typedef location_str LOCATION;
// A cell LOCATION has position (x,y,z), radius (r)

#endif
