#include <stdio.h>
#include <stddef.h>
#include <math.h>
#include <time.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

enum { JAN, FEB, MAR, APR, MAY, JUN, JUL, AUG, SEP, OCT, NOV, DEC };

enum { RAD, DEG };

typedef struct { float m[2][2];     } matrix_2x2f_t;
typedef struct { float m[3][3];     } matrix_3x3f_t;
typedef struct { double m[2][2];    } matrix_2x2d_t;
typedef struct { double m[3][3];    } matrix_3x3d_t;

typedef struct { float x, y;        } vec_2df_t;
typedef struct { float x, y, z;     } vec_3df_t;
typedef struct { double x, y;       } vec_2dd_t;
typedef struct { double x, y, z;    } vec_3dd_t;

typedef struct { float lat, lon;    } spherical_f_t;
typedef struct { double lat, lon;   } spherical_d_t;
typedef struct { double alt, az;    } spherical_alt_az_d_t;

double vec_3dd_sq_mag(vec_3dd_t * v) { return v->x*v->x + v->y*v->y + v->z*v->z; }
double vec_3dd_mag(vec_3dd_t * v) { return sqrt(vec_3dd_sq_mag(v)); }

double rad2deg(double rad) { return rad*180.0/M_PI; }
double deg2rad(double deg) { return deg*M_PI/180.0; }

void vec_3dd_to_str(vec_3dd_t * v, char * dst, size_t n) { snprintf(dst, n, "< %f, %f, %f >", v->x, v->y, v->z); }

void spherical_d_to_str(spherical_d_t * s, int fmt, char * dst, size_t n) {
    double lat = s->lat, lon = s->lon;
    if (fmt == DEG) {
        lat = rad2deg(lat);
        lon = rad2deg(lon);
    }
    snprintf(dst, n, "(latitude: %f, longitude: %f)", lat, lon);
}

void spherical_alt_az_d_to_str(spherical_alt_az_d_t * s, int fmt, char * dst, size_t n) {
    double alt = s->alt, az = s->az;
    if (fmt == DEG) {
        alt = rad2deg(alt);
        az = rad2deg(az);
    }
    snprintf(dst, n, "(altitude: %f, azimuth: %f)", alt, az);
}

vec_3dd_t
spherical_to_cartesian(spherical_d_t * s)
{
    vec_3dd_t r;
#if 0
    // where lat is measured from equator
    r.x = cos(s->lat)*cos(s->lon);
    r.y = cos(s->lat)*sin(s->lon);
    r.z = sin(s->lat);
#else
    // where lat is measured from north pole
    r.x = sin(s->lat)*cos(s->lon);
    r.y = sin(s->lat)*sin(s->lon);
    r.z = cos(s->lat);
#endif
    return r;
}

// assumes unit sphere
spherical_d_t
cartesian_to_spherical(vec_3dd_t * v)
{
    spherical_d_t r;
#if 0
    // where lat is measured from equator
    r.lat = asin(v->z);
    r.lon = atan2(v->y, v->x);
#else
    // where lat is measured from north pole
    r.lat = acos(v->z);
    r.lon = atan2(v->y, v->x);
#endif
    return r;
}

vec_3dd_t
matmult_vec_3dd(matrix_3x3d_t * mat, vec_3dd_t * vec)
{
    vec_3dd_t r;
    r.x = mat->m[0][0]*vec->x + mat->m[0][1]*vec->y + mat->m[0][2]*vec->z;
    r.y = mat->m[1][0]*vec->x + mat->m[1][1]*vec->y + mat->m[1][2]*vec->z;
    r.z = mat->m[2][0]*vec->x + mat->m[2][1]*vec->y + mat->m[2][2]*vec->z;
    return r;
}

void
rotate_2df(vec_2df_t * v, float theta)
{
}

void
rotate_3df_x(vec_3df_t * v, float theta)
{
}

void
rotate_3df_y(vec_3df_t * v, float theta)
{
}

void
rotate_3df_z(vec_3df_t * v, float theta)
{
}

void
rotate_2dd(vec_2dd_t * v, double theta)
{
}

vec_3dd_t
rotate_3dd_x(vec_3dd_t * v, double theta)
{
    matrix_3x3d_t rot_mat = {
        {
            { 1.0, 0.0,        0.0         },
            { 0.0, cos(theta), -sin(theta) },
            { 0.0, sin(theta), cos(theta)  }
        }
    };

    return matmult_vec_3dd(&rot_mat, v);
}

vec_3dd_t
rotate_3dd_y(vec_3dd_t * v, double theta)
{
    matrix_3x3d_t rot_mat = {
        {
            { cos(theta),  0.0, sin(theta) },
            { 0.0,         1.0, 0.0        },
            { -sin(theta), 0.0, cos(theta) }
        }
    };

    return matmult_vec_3dd(&rot_mat, v);
}

vec_3dd_t
rotate_3dd_z(vec_3dd_t * v, double theta)
{
    matrix_3x3d_t rot_mat = {
        {
            { cos(theta), -sin(theta), 0.0 },
            { sin(theta), cos(theta),  0.0 },
            { 0.0,         0.0,        1.0 }
        }
    };

    return matmult_vec_3dd(&rot_mat, v);
}

// TODO: should probably return in radians to be consistent
spherical_alt_az_d_t
alt_az_of_one_point_from_another(spherical_d_t p0, spherical_d_t p1)
{
    spherical_alt_az_d_t r;
    // north pole is 0
    p0.lat = M_PI/2.0 - p0.lat;
    p1.lat = M_PI/2.0 - p1.lat;

    vec_3dd_t v = spherical_to_cartesian(&p1);

    // apply transformations which would be necessary to move p0
    // to the north pole with 0 longitude pointing the same way
    v = rotate_3dd_z(&v, M_PI - p0.lon); // 0 is toward north pole
    v = rotate_3dd_y(&v, p0.lat);

    p1 = cartesian_to_spherical(&v);
    p1.lat = M_PI/2.0 - p1.lat; // horizon is 0
    //p1.lat = rad2deg(p1.lat);
    //p1.lon = 360.0 - rad2deg(p1.lon); // east is 90
    p1.lon = 2*M_PI - p1.lon;
    //if (p1.lon >= 360.0) // 359.999999999999900000
    //    p1.lon -= 360.0;
    if (p1.lon >= 2*M_PI) // 359.999999999999900000
        p1.lon -= 2*M_PI;
    r.alt = p1.lat;
    r.az  = p1.lon;
    return r;
}

#if 0
static struct tm vernal_equinox_2021 = {
    .tm_year = 2021 - 1900,
    .tm_mon = MAR,
    .tm_mday = 20,
    .tm_hour = 9,
    .tm_min = 37,
    .tm_sec = 0,
    .tm_isdst = -1
};
#endif

time_t vernal_equinox_2021_time_t = 1616251020; // mktime(&vernal_equinox_2021);

// returns coordinates in celestial sphere
spherical_d_t
sun_celestial_coordinate_at_time(time_t t)
{
    double dt = difftime(t, vernal_equinox_2021_time_t);
    double dt_years = dt/(60*60*24*365);
    // TODO: the equations below assume a circular orbit
    spherical_d_t r = {
        .lat = deg2rad(23.4) * sin(dt_years*2*M_PI),
        .lon = deg2rad(dt_years*360.0)
    };
    return r;
}

// returns coordinates in geographic coordinate system
spherical_d_t
sun_geographic_coordinate_at_time(time_t t)
{
    double dt = difftime(t, vernal_equinox_2021_time_t);

    double init_lon = deg2rad((360*(12.0 - (9+37.0/60.0))/24.0));
    double dt_days = dt/(60*60*24);
    double dt_years = dt_days/365;
    // TODO: the equations below assume a circular orbit
    spherical_d_t r = {
        .lat = deg2rad(23.4) * sin(dt_years*2*M_PI),
        //.lon = init_lon - deg2rad(dt_days*(360.0 - 360.0/365.0))
        .lon = init_lon - deg2rad(dt_days*360.0)
    };
    while (r.lon < 0)
        r.lon += 2*M_PI;
    while (r.lon > 2*M_PI)
        r.lon -= 2*M_PI;
    return r;
}

#if 0
time_t      time(time_t *);
char *      ctime(const time_t *);
char *      asctime(const struct tm *);
double      difftime(time_t, time_t);
time_t      mktime(struct tm *);
struct tm * localtime(const time_t *);
struct tm * gmtime(const time_t *);
#endif

int main()
{
    char buffer[1024];

    // resource: https://sunrise-sunset.org/us/austin-tx
    // resource: https://sun-direction.com/

    time_t curr_time;
    //struct tm * timeinfo;

    time(&curr_time);
    //timeinfo = localtime(&curr_time);
    printf("Current local time and date: %s", ctime(&curr_time));
    //curr_time += 5*60*60; // TODO: depends on time zone
    printf("Current UTC time and date: %s", asctime(gmtime(&curr_time)));

    //for (int i = -3; i < 10 ; i++) {
    for (int i = 0; i < 14*4 ; i++) {
        // TODO: convert date and time to seconds since time 0
        spherical_d_t sun = sun_geographic_coordinate_at_time(curr_time + i*15*60);
        //spherical_d_t sun = sun_geographic_coordinate_at_time((197*24*60*60) + (((11-9)*60)+i*15-37)*60);
        //spherical_d_t sun = sun_geographic_coordinate_at_time(197*24*60*60);
        //spherical_d_t sun = sun_geographic_coordinate_at_time(i*60*60);
        spherical_d_to_str(&sun, DEG, buffer, sizeof(buffer));
        int hr = 6 + i*15/60;
        if (hr > 12)
            hr -= 12;
        int min = (i*15)%60;
        printf("time: %d:%02d sun position:\t%s\t", hr, min, buffer);

        spherical_d_t me  = { .lat = deg2rad(30.0), .lon = deg2rad(263) }; // Austin, TX
        spherical_alt_az_d_t x = alt_az_of_one_point_from_another(me, sun);
        spherical_alt_az_d_to_str(&x, DEG, buffer, sizeof(buffer));
        printf("%s\n", buffer);
    }

    return 0;
}
