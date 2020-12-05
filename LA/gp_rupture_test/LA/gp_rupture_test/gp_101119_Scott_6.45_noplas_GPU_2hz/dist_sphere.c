/* function to calculate distance between two points with
   latitude and longitude (lat1, lon1) and (lat2, lon2) on
   a sphere   

   source: http://www.fcaglp.unlp.edu.ar/~esuarez/gmt/1997/0148.html
  
   Daniel Roten, 2005-01-27
*/

#include <math.h>
float dist_sphere(float lat1, float lon1, float lat2, float lon2, 
                   float R){
   float dlon, dlat, a, c, d;
  
   /* convert to radians */
   lat1=lat1*M_PI/180;
   lat2=lat2*M_PI/180;
   lon1=lon1*M_PI/180;
   lon2=lon2*M_PI/180;

   dlon = lon2-lon1;
   dlat = lat2-lat1;
   a=pow(sin(dlat/2),2) + cos(lat1) * cos(lat2) * pow(sin(dlon/2),2);
   if (sqrt(a) < 1)
      c = 2 * asin(sqrt(a));
   else
      c = 2 * asin(1); /* protects against possible roundoff errors 
                          if points are nearly antipodal) */
   d = R*c;
   return d;
}
