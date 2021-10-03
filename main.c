#include <stdio.h>
#include <stddef.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

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

void spherical_d_to_str(spherical_d_t * s, char * dst, size_t n) { snprintf(dst, n, "(latitude: %f, longitude: %f)", s->lat, s->lon); }
void spherical_alt_az_d_to_str(spherical_alt_az_d_t * s, char * dst, size_t n) { snprintf(dst, n, "(altitude: %f, azimuth: %f)", s->alt, s->az); }

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
    p1.lat = rad2deg(p1.lat);
    p1.lon = 360.0 - rad2deg(p1.lon); // east is 90
    if (p1.lon >= 360.0) // 359.999999999999900000
        p1.lon -= 360.0;
    r.alt = p1.lat;
    r.az  = p1.lon;
    return r;
}

// time: seconds since 9:37am UT on 20 March 2021
// returns coordinates in celestial sphere
spherical_d_t
sun_celestial_coordinate_at_time(int time)
{
    double time_years = ((double)time)/(60*60*24*365);
    // TODO: the equations below assume a circular orbit
    spherical_d_t r = {
        .lat = 23.4 * sin(time_years*2*M_PI),
        .lon = deg2rad(time_years*360.0)
    };
    return r;
}

// time: seconds since 9:37am UT on 20 March 2021
// returns coordinates in geographic coordinate system
spherical_d_t
sun_geographic_coordinate_at_time(int time)
{
    double init_lon = deg2rad((360*(12.0 - (9+37.0/60.0))/24.0));
    double days = ((double)time)/(60*60*24);
    double years = ((double)days)/(365);
    // TODO: the equations below assume a circular orbit
    spherical_d_t r = {
        .lat = deg2rad(23.4) * sin(years*2*M_PI),
        //.lon = init_lon - deg2rad(days*(360.0 - 360.0/365.0))
        .lon = init_lon - deg2rad(days*360.0)
    };
    while (r.lon < 0)
        r.lon += 2*M_PI;
    while (r.lon > 2*M_PI)
        r.lon -= 2*M_PI;
    return r;
}

int main()
{
    char buffer[1024];

#if 0
    spherical_d_t s;
    vec_3dd_t init_v = {0.5, 0.5, sqrt(2.0)/2.0};
    vec_3dd_t init_v2 = {sqrt(2.0)/2.0, sqrt(2.0)/2.0, 0.0};

    // print original
    vec_3dd_t v = init_v;
    // print vector
    vec_3dd_to_str(&v, buffer, sizeof(buffer));
    printf("original:                 %s, ", buffer);
    // print spherical
    s = cartesian_to_spherical(&v);
    s.lat = rad2deg(s.lat);
    s.lon = rad2deg(s.lon);
    spherical_d_to_str(&s, buffer, sizeof(buffer));
    printf("%s\n", buffer);

    // print original after being rotated 90 deg on Z-axis
    v = rotate_3dd_z(&init_v, M_PI/2); // rotate 90 deg on Z-axis
    // print vector
    vec_3dd_to_str(&v, buffer, sizeof(buffer));
    printf("rotated 90 deg on Z-axis: %s, ", buffer);
    // print spherical
    s = cartesian_to_spherical(&v);
    s.lat = rad2deg(s.lat);
    s.lon = rad2deg(s.lon);
    spherical_d_to_str(&s, buffer, sizeof(buffer));
    printf("%s\n", buffer);

    // print original v2
    v = init_v2;
    // print vector
    vec_3dd_to_str(&v, buffer, sizeof(buffer));
    printf("original v2:              %s, ", buffer);
    // print spherical
    s = cartesian_to_spherical(&v);
    s.lat = rad2deg(s.lat);
    s.lon = rad2deg(s.lon);
    spherical_d_to_str(&s, buffer, sizeof(buffer));
    printf("%s\n", buffer);

    // print original v2 after being rotated 90 deg on Y-axis
    v = rotate_3dd_y(&init_v2, deg2rad(90.0)); // rotate 90 deg on Y-axis
    // print vector
    vec_3dd_to_str(&v, buffer, sizeof(buffer));
    printf("rotated 90 deg on Y-axis: %s, ", buffer);
    // print spherical
    s = cartesian_to_spherical(&v);
    s.lat = rad2deg(s.lat);
    s.lon = rad2deg(s.lon);
    spherical_d_to_str(&s, buffer, sizeof(buffer));
    printf("%s\n", buffer);

    // print original v2 after being rotated 90 deg on Y-axis
    v = rotate_3dd_x(&init_v2, M_PI/2); // rotate 90 deg on X-axis
    // print vector
    vec_3dd_to_str(&v, buffer, sizeof(buffer));
    printf("rotated 90 deg on X-axis: %s, ", buffer);
    // print spherical
    s = cartesian_to_spherical(&v);
    s.lat = rad2deg(s.lat);
    s.lon = rad2deg(s.lon);
    spherical_d_to_str(&s, buffer, sizeof(buffer));
    printf("%s\n", buffer);
#endif

    // resource: https://sunrise-sunset.org/us/austin-tx

    // TODO: get sun position from time
    // TODO: get me position from time

    //for (int i = -3; i < 10 ; i++) {
    for (int i = 0; i < 14*4 ; i++) {
        // TODO: convert date and time to seconds since time 0
        spherical_d_t sun = sun_geographic_coordinate_at_time((197*24*60*60) + (((11-9)*60)+i*15-37)*60);
        //spherical_d_t sun = sun_geographic_coordinate_at_time(197*24*60*60);
        //spherical_d_t sun = sun_geographic_coordinate_at_time(i*60*60);
        sun.lat = rad2deg(sun.lat);
        sun.lon = rad2deg(sun.lon);
        spherical_d_to_str(&sun, buffer, sizeof(buffer));
        int hr = 6 + i*15/60;
        if (hr > 12)
            hr -= 12;
        int min = (i*15)%60;
        printf("time: %d:%02d sun position:\t%s\t", hr, min, buffer);
        sun.lat = deg2rad(sun.lat);
        sun.lon = deg2rad(sun.lon);

        spherical_d_t me  = { .lat = deg2rad(30.0), .lon = deg2rad(263) };
        //spherical_d_t me  = { .lat = deg2rad(0.0), .lon = deg2rad(0.0) };
        spherical_alt_az_d_t x = alt_az_of_one_point_from_another(me, sun);
        spherical_alt_az_d_to_str(&x, buffer, sizeof(buffer));
        printf("%s\n", buffer);
    }

    return 0;
}
