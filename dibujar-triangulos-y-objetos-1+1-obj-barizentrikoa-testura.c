//	Program developed by
//	
//	Informatika Fakultatea
//	Euskal Herriko Unibertsitatea
//	http://www.ehu.eus/if
//
// to compile it: gcc dibujar-triangulos-y-objetos.c -lGL -lGLU -lglut
//
// 
//


#include <GL/glut.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
//#include "cargar-triangulo.h"
#include "obj.h"


// testuraren informazioa
// información de textura

extern int load_ppm(char *file, unsigned char **bufferptr, int *dimxptr, int * dimyptr);
unsigned char *bufferra;
int dimx,dimy,dimentsioa;

int indexx;
//hiruki *triangulosptr;
object3d *foptr;
object3d *sel_ptr;
int denak;
int lineak;
int objektuak;
int kamera;
char aldaketa;
int ald_lokala;
int objektuaren_ikuspegia;
int paralelo;
char fitxiz[100];
int atzearpegiakmarraztu;

typedef struct  
{
    int n; //number of points in the row 
    int x[10]; //points
    float u[10]; //corresponding u coordinates
    float v[10]; //corresponding v coordinates
            
} ScanLine; //cut points
 

ScanLine *intersectionTable = NULL; // Dynamic Table [dimentsioa]



// TODO  m = m1 * m2 
// m1 bider m2 matrizeen biderketa egin eta m matrizean jaso. 
void mxm(double *m, double*m1, double *m2)
{
}


// TODO
void kamerari_aldaketa_sartu_ezk(double *m)
{
}

// TODO
void kamerari_aldaketa_sartu_esk(double m[16])
{
}

// TODO
void objektuari_aldaketa_sartu_ezk(double m[16])
{
}

// TODO
void objektuari_aldaketa_sartu_esk(double m[16])
{
}



// TODO given u,v get the color pointer
unsigned char *color_textura(float u, float v)
{
int indx,indy;
char * lag;
//printf("texturan...%x\n",bufferra);
//TODO get the desplacement for indx and indy

// negative values?
if (u < 0) {u = 0;}
if (v < 0) {v = 0;}

// values greater than 1?
if (u > 1) {u = 1;}
if (u > 1) {u = 1;}

// Convert to texel coordinates
indx= (int) (u * (dimx-1) + 0.5f);
indy= (int) (v * (dimy-1) + 0.5f);

// Invert V (PPM)
indy = (dimy-1) - indy;

// Points to the first texel (R)
lag = (unsigned char *)bufferra;
//printf("irtetera %x\n",lag[indy*dimx+indx]);
return(lag+3*(indy*dimx+indx));

}



void print_matrizea16(double *m)
{
int i;

for (i = 0;i<4;i++)
   printf("%lf, %lf, %lf, %lf\n",m[i*4],m[i*4+1],m[i*4+2], m[i*4+3]);
}


// TODO  res = m * v 
// v bektoreari m matrizea bidertu eta res erakusleak adierazten duen bektorean jaso. 
// v bektorearen eta res emaitzaren laugarren osagaia 0 dela suposatzen du.
void mxv(double *res, double *m, double *v)
{
    res[0] = v[0];
    res[1] = v[1];
    res[2] = v[2];
}



// TODO  pptr = m * p 
// ppuntuari m matrizea bidertu eta pptr erakulseak adierazten duen puntuan jaso. 
// p puntuaren laugarren osagaia 1 dela suposatzen du.
// matrizearen laugarren lerroaren arabera emaitzaren laugarren osagaia, w, ez bada 1, orduan bere baliokidea itzuli behar du: x/w, y/w eta z/w
void mxp(point3 *pptr, double m[16], point3 p)
{
    pptr->x = p.x;
    pptr->y = p.y;  
    pptr->z = p.z;
}


// TODO objektuaren erpinek eta bektore normalek kameraren erreferentzi-sisteman dituzten koordenatuak lortu
void kam_ikuspegia_lortu(object3d *optr)
{
int i;

// TODO  get point in the viewer coordenate-system and project it
// TODO get the vectors in camera system-
for (i = 0; i<optr->num_vertices; i++) 
    { 
    //TODO aldatu
    // get viewer coordinates
    optr->vertex_table[i].camcoord.x= optr->vertex_table[i].coord.x;
    optr->vertex_table[i].camcoord.y= optr->vertex_table[i].coord.y;
    optr->vertex_table[i].camcoord.z= optr->vertex_table[i].coord.z;
    // TODO aldatu
    // Get projected coordinates
    optr->vertex_table[i].proedcoord.x = optr->vertex_table[i].camcoord.x;
    optr->vertex_table[i].proedcoord.y = optr->vertex_table[i].camcoord.y;
    optr->vertex_table[i].proedcoord.z = optr->vertex_table[i].camcoord.z;
    // TODO aldatu
    // get normal vector in camera coordinates
    optr->vertex_table[i].Ncam[0] = optr->vertex_table[i].N[0];
    optr->vertex_table[i].Ncam[1] = optr->vertex_table[i].N[1];
    optr->vertex_table[i].Ncam[2] = optr->vertex_table[i].N[2];
    }
      
// TODO get the face normal in the viewer coordenate-system 
for (i = 0; i<optr->num_faces; i++) 
    {
    // TODO aldatu
    optr->face_table[i].Ncam[0] = optr->face_table[i].N[0];
    optr->face_table[i].Ncam[1] = optr->face_table[i].N[1];
    optr->face_table[i].Ncam[2] = optr->face_table[i].N[2];
    }
}



void modelview_lortu(double *m1, double *m2)
{
}

void mesa_lortu(double* M) 
{
}


void argien_kalkulua_egin(object3d *optr, int ti)
{

}

// Draws a pixel on the screen at the given x and y coordinates
static void putpixel (float x,float y){

    glBegin(GL_POINTS);
    glVertex2f(x,y);
    glEnd();
}

// Initialize the entire table to 0's
static void initTable () {
    for (int i = 0; i < dimentsioa; i++){
        intersectionTable[i].n = 0;
        for (int j = 0; j < 10; j++ ){
            intersectionTable[i].x[j] = 0;
            intersectionTable[i].u[j] = 0.0f;
            intersectionTable[i].v[j] = 0.0f;
        }
    }
}   

// Sorts the array of x values in each table position (y) from lowest to highest
// Needed because the edge algorithm stores them in the order they are found, which can be chaotic
static void sortTable (){

    /* Depending on the value of n (number of x's for the current y):
        if n == 0 → continue (no intersection with the triangle)
        if n == 1 → continue (vertical line case, already handled during drawing)
        if n == 2:
            - If both x values are equal → continue (single point)
            - Otherwise → sort
        if n > 2 → sort
    */
   
    for (int y = 0; y < dimentsioa; y++){

        if (intersectionTable[y].n < 2) continue;

        // n >= 2
        if (intersectionTable[y].n == 2 && (intersectionTable[y].x[0] == intersectionTable[y].x[1])) continue;
        
        for (int i = 0; i < intersectionTable[y].n - 1; i++) {
            for (int j = i + 1; j < intersectionTable[y].n; j++) {
                if (intersectionTable[y].x[i] > intersectionTable[y].x[j]) {

                    // Swap X
                    int tempX = intersectionTable[y].x[i];
                    intersectionTable[y].x[i] = intersectionTable[y].x[j];
                    intersectionTable[y].x[j] = tempX;

                    // Swap U
                    float tempU = intersectionTable[y].u[i];
                    intersectionTable[y].u[i] = intersectionTable[y].u[j];
                    intersectionTable[y].u[j] = tempU;

                     // Swap V
                    float tempV = intersectionTable[y].v[i];
                    intersectionTable[y].v[i] = intersectionTable[y].v[j];
                    intersectionTable[y].v[j] = tempV;

                }
            }
        }
    }   
}

// Store an intersection point in our table
static void storePointT (int X, int Y,float u,float v ){

    int idx = intersectionTable[Y].n;

    //Avoid exceeding maximum number of x's values (10)
    if (idx >= 10) return; 

    intersectionTable[Y].x[idx] = X;
    intersectionTable[Y].u[idx] = u;
    intersectionTable[Y].v[idx] = v;
    intersectionTable[Y].n++;    
}

// 1-Simplex (Line 1D)
void linearInterpolationUV ( 
    float x0, float y0, float u0, float v0, 
    float x1, float y1, float u1, float v1,   
    float x, float y,                         
    float *u, float *v) 
{
    float t;

    float dx = fabsf(x1 - x0); // Δx
    float dy = fabsf(y1 - y0); // Δy
    float length = fmaxf(dx, dy); 

    // dx == dy
    if (length == 0.0f) {
        *u = u0;
        *v = v0;
        return;
    }

    if (dx >= dy) 
        t = (x - x0) / (x1 - x0);
    else
        t = (y - y0) / (y1 - y0);

    // Clamp t to [0, 1]
    if (t < 0.0f) t = 0.0f;
    if (t > 1.0f) t = 1.0f;

    // Texture coordinates
    *u = u0 + t * (u1 - u0);
    *v = v0 + t * (v1 - v0);
}

// float coordinate [-1,1] to integer coordinate [0, dimmentsioa-1] 
static int toInt (float c){

    /* For any number in a range [a,b] that you want to map [A,B]

                            B - A   
            c' = A + (x-a) -------
                            b - a
    */
    int c_p;

    c_p = round((c + 1.0f) * (dimentsioa - 1) / (2.0f));

    return c_p;
}

// integer coordinate [0, dimentsioa-1] to float coordinate [-1,1]
static float toFloat (int c){
    float c_p;

    c_p = c * 2.0f / (dimentsioa - 1) - 1.0f;

    return roundf(c_p * 1000000.0f) / 1000000.0f; //rounded to 6 decimals
}

// Bresenham's Algorithm
void drawLine (float x0,float y0,float x1,float y1){
    int x_0,y_0,x_1,y_1,dx,dy;
    int offset; //direction
    float P; //Parameter decision

    //Convert to intengers 
    x_0 = toInt(x0);
    y_0 = toInt(y0);
    x_1 = toInt(x1);
    y_1 = toInt(y1);
    
    //Change in x and y
    dx = abs(x_1 - x_0); //Δx
    dy = abs(y_1 - y_0); //Δy

    offset = x1 - x0 > 0 ? 1 : -1; // R = 1, L = -1
    
    if (dx == 0){ //Special Case (Vertical Line)

        for (int i = 0; i < dy; i++){
            putpixel(toFloat(x_0), toFloat(y_0 + (-1) * i));
        }
    } 
    
    else if (dx >= dy) { //Case 1: 0 <= m <= 1 (gentle slope,line closer to X-axis)

        //Initializing
        P = 2*dy - dx;
        for (int i = 0; i < dx; i++){  
            putpixel(toFloat(x_0 + offset*i), toFloat(y_0));
 
            if (P > 0){ // Threshold reached
                y_0-=1; // Move one pixel down
                P = P - 2*dx; //Reduce
            }

            P = P + 2*dy; //Increment
        }

    } else { //Case 2: dx < dy | m > 1 (steep slope, line closer to Y-axis)

        P = 2*dx - dy;
        for (int i = 0; i < dy; i++){
            putpixel(toFloat(x_0), toFloat(y_0 + (-1)*i));

            if (P > 0){
                x_0+=offset;
                P = P - 2*dy;
            }

            P = P + 2*dx;       
        }
    }
}

// drawLine function, but With Texture mapping and Scanline point storage
void drawLineWT (float x0,float y0,float u0,float v0,float x1,float y1,float u1, float v1,object3d *optr){
    int x_0,y_0,x_1,y_1,dx,dy;
    int offset; //direction
    float P; //Parameter decision
    float u,v; 
    float r,g,b; //Components Red - Green - Blue
    unsigned char *texColor;

    //Convert to intengers  
    x_0 = toInt(x0);
    y_0 = toInt(y0);
    x_1 = toInt(x1);
    y_1 = toInt(y1);
    
    //Change in x and y
    dx = abs(x_1 - x_0); //Δx
    dy = abs(y_1 - y_0); //Δy

    offset = x1 - x0 > 0 ? 1 : -1; // R = 1, L = -1
    
    if (dx == 0){ //Special Case (Vertical Line)
         
    /*  Do not skip vertical edges: must be stored, 
        as they can serve as the left or right boundary 
        of the triangle
    */
        for (int i = 0; i < dy; i++){
            
            if (optr->texturaduna){ //If the object has a texture, compute texture coords for each point

            linearInterpolationUV(x0,y0,u0,v0,x1,y1,u1,v1,toFloat(x_0),toFloat(y_0),&u,&v);
            texColor =  color_textura(u,v);

            // Normalize RGB values to [0,1]
            r = texColor[0]/255.0f;
            g = texColor[1]/255.0f;
            b = texColor[2]/255.0f;

            // Change the color
            glColor3f(r,g,b);   
        }
            putpixel(toFloat(x_0),toFloat(y_0));
            storePointT(x_0,y_0,u,v);
            y_0-= 1; //Decrement y along the Y-axis

        }
    } 
    
    else if (dx >= dy) { //Case 1: 0 <= m <= 1 (gentle slope,line closer to X-axis)

        //Initializing
        P = 2*dy - dx;
        int prevY = -1;// keep track of last Y stored 
        for (int i = 0; i < dx; i++){  

            if (optr->texturaduna){ //If the object has a texture, compute texture coords for each point

            linearInterpolationUV(x0,y0,u0,v0,x1,y1,u1,v1,toFloat(x_0),toFloat(y_0),&u,&v);
            texColor =  color_textura(u,v);
            
            // Normalize RGB values to [0,1]
            r = texColor[0]/255.0f;
            g = texColor[1]/255.0f;
            b = texColor[2]/255.0f;

            // Change the color 
            glColor3f(r,g,b);
        }
            putpixel(toFloat(x_0),toFloat(y_0));

            // If the line is horizontal or low-slope (dy == 0 or very small) 
            // avoid storing the same point multiple times
            if (dy != 0 && y_0 != prevY) {//
                storePointT(x_0,y_0,u,v);
                prevY = y_0;
            }

            if (P > 0){ // Threshold reached
                y_0-=1; // Move one pixel down
                P = P - 2*dx; //Reduce
            }

            P = P + 2*dy; //Increment the ERROR committed in each iteration
            x_0+=offset;
        }

    } else { //Case 2: dx < dy | m > 1 (steep slope, line closer to Y-axis)

        P = 2*dx - dy;
        for (int i = 0; i < dy; i++){

            if (optr->texturaduna){ //If the object has a texture, compute texture coords for each point

            linearInterpolationUV(x0,y0,u0,v0,x1,y1,u1,v1,toFloat(x_0),toFloat(y_0),&u,&v);
            texColor =  color_textura(u,v);

            // Normalize RGB values to [0,1]
            r = texColor[0]/255.0f;
            g = texColor[1]/255.0f;
            b = texColor[2]/255.0f;

            // Change the color
            glColor3f(r,g,b);
        }
            putpixel(toFloat(x_0),toFloat(y_0));
            storePointT(x_0,y_0,u,v);
            
            if (P > 0){
                x_0+=offset;
                P = P - 2*dy;
            }

            P = P + 2*dx; 
            y_0-= 1;
        }
    }
}

// Scan Line Algorithm
static void drawInternalPoints (float x0, float y0, float x1, float y1, object3d *optr){

    int x_0,y_0,x_1,y_1;
    int i_x0,i_x1; //x of intersection points
    float i_u0,i_v0,i_u1,i_v1; //intersection texture coordinates
    float u,v;
    float r,g,b; //Components Red - Green - Blue
    int lastIndex; 
    unsigned char *texColor;

    //Convert to intengers  
    x_0 = toInt(x0);
    y_0 = toInt(y0); 
    x_1 = toInt(x1);
    y_1 = toInt(y1);

    // Horizontal line or identical point: no need to draw since it's already rendered
    if (y_0 == y_1) return;

    for (int i = y_0-1;i >= y_1;i--){

        if (intersectionTable[i].n < 2) continue;

        // Identical points (probably one of the vertices)
        if (intersectionTable[i].n == 2 && intersectionTable[i].x[0] == intersectionTable[i].x[1]) continue;

        lastIndex = intersectionTable[i].n-1; 

        // Prepare intersection points 
        i_x0 = intersectionTable[i].x[0];
        i_u0 = intersectionTable[i].u[0];
        i_v0 = intersectionTable[i].v[0];
        i_x1 = intersectionTable[i].x[lastIndex];
        i_u1 = intersectionTable[i].u[lastIndex];
        i_v1 = intersectionTable[i].v[lastIndex];

        // Draw the horizontal span between both intersection points
        for (int j = i_x0+1;j < i_x1;j++){

            if (optr->texturaduna){ //If the object has a texture, compute texture coords for each point

            linearInterpolationUV(toFloat(i_x0),i,i_u0,i_v0,toFloat(i_x1),i,i_u1,i_v1,toFloat(j),toFloat(i),&u,&v);
            texColor =  color_textura(u,v);

            // Normalize RGB values to [0,1]
            r = texColor[0]/255.0f;
            g = texColor[1]/255.0f;
            b = texColor[2]/255.0f;

            // Change the color
            glColor3f(r,g,b);   
    }
            putpixel(toFloat(j),toFloat(i));
        }
        
    }
} 


/* ti is the face index
** i1 if the face has more than 3 vertices there will be more than 1 triangle, 
** i1 indicates the ith triangle of the face:
**  a.- The first vertex and the next two vertices form the first triangle, so
        0, 1, 2 are the indices of the vertices and i1 will be 1
    b.- the first vertex and the third and fourth vertices form the next triangle. So
        0, 2, 3 are the indices and consequently 1i will be 2 (second triangle of the 
        face)
    c.- 0, 3, 4 form the next triangle and so i1= 3...
** atzeaurpegiada indicates that the face is a backface. Depending of the state of the aplication tha face will be drawn in red or it will not be drawn
*/
void dibujar_triangulo(object3d *optr, int ti, int i1, int atzeaurpegiada)
{
point3 *pgoiptr, *pbeheptr, *perdiptr;
float x, y, z, u, v, nx, ny;
float luz, l12, l23, l13, d1v23, d2v13, d3v12;
float lerrotartea,pixeldist;
float alfa, beta, gamma;
float goieragina, beheeragina, erdieragina;
float *luzeptr,*erdiptr, *motzptr;
float  luzeluz, paraleloarenluzera;
int luzeanparalelokop, barnekop;
int lerrokop,i,j;
int ind0, ind1, ind2, indg, inde, indb;
point3 *p1ptr, *p2ptr, *p3ptr;
double *Nkam;
double aldaketa,baldaketa;
float x1,x2, y1;
unsigned char r,g,b;
unsigned char *colorv;

// triangeluaren hiru erpinak hartu
ind0 = optr->face_table[ti].vertex_ind_table[0];
ind1 = optr->face_table[ti].vertex_ind_table[i1];
ind2 = optr->face_table[ti].vertex_ind_table[i1+1];


p1ptr = &(optr->vertex_table[ind0].proedcoord);
p2ptr = &(optr->vertex_table[ind1].proedcoord);
p3ptr = &(optr->vertex_table[ind2].proedcoord);
       
// Puntuz puntu marraztuko dut dena!!     
// hasteko bi pixelen arteko distantzia gure munduan (-1 eta 1 arteko munduan) ze distantzia den kalkulatuko dut       
pixeldist = 2.0/(float)dimentsioa;


 // lehenengo hiru erpinak ordenatu behar ditut 
//TODO erpinak ordenatu!!! (segun la coordena Y)

typedef struct {
    point3 *p;  
    int index;    
} Vertex;

Vertex vx[3];
Vertex aux;

// Add point and index of each vertex 
vx[0].p = p1ptr; vx[0].index = ind0;
vx[1].p = p2ptr; vx[1].index = ind1;
vx[2].p = p3ptr; vx[2].index = ind2;    

// Find the points v[0] > v[1] > v[2]
if (vx[0].p->y < vx[1].p->y){
    aux = vx[0]; vx[0] = vx[1]; vx[1] = aux;
}
if (vx[0].p->y < vx[2].p->y){
    aux = vx[0]; vx[0] = vx[2]; vx[2] = aux;
}
if (vx[1].p->y < vx[2].p->y){
    aux = vx[1]; vx[1] = vx[2]; vx[2] = aux;
}

// Assign ordered points
pgoiptr = vx[0].p;   perdiptr = vx[1].p;  pbeheptr = vx[2].p; 
indg = vx[0].index;  inde = vx[1].index;  indb = vx[2].index;


r =optr->face_table[ti].rgb[0]; //r:255
g= optr->face_table[ti].rgb[1]; //g:255
b= optr->face_table[ti].rgb[2]; //b:255
glColor3ub(r,g,b); //White

 //TODO draw the three vertices. 3 erpinak marraztu
    glBegin( GL_POINTS ); 

    //glColor3ub(r,g,b);

    /*In C:
        (value ≠ 0) → True
        (value == 0)→ False

        Examples:
        - !0 → 1 (True)
        - !1 - 0 (False) 
        - !25 → 0 (False) 
        - !(-4) → 0 (False)

        If atzeaurpegiada == 0 → !atzeaurpegiada (1) → True  (not back face | front face)
        If atzeaurpegiada == 1 → !atzeaurpegiada (0) → False (back face)
    */

    if(!atzeaurpegiada)
           {
            if (optr->texturaduna) //Object has texture → texturaduna == 1
                {
                colorv = color_textura(optr->vertex_table[ind0].u, optr->vertex_table[ind0].v); //Vertex color
                glColor3ub(colorv[0],colorv[1],colorv[2]); 
                }
              else //No texture
                {            
                glColor3ub(r,g,b); //White color (default) 
           }
        }
    glVertex3f(p1ptr->x, p1ptr->y, p1ptr->z );
    if(!atzeaurpegiada)
           {
            if (optr->texturaduna) 
                {
                colorv = color_textura(optr->vertex_table[ind1].u, optr->vertex_table[ind1].v);
                glColor3ub(colorv[0],colorv[1],colorv[2]);
                }
              else
                {            
                glColor3ub(r,g,b);
                }
          /* glColor3ub(optr->vertex_table[ind1].rgb[0],
                      optr->vertex_table[ind1].rgb[1],
                      optr->vertex_table[ind1].rgb[2]);*/
           }
    glVertex3f(p2ptr->x, p2ptr->y, p2ptr->z );
    if(!atzeaurpegiada)
           {
            if (optr->texturaduna)
                {
                colorv = color_textura(optr->vertex_table[ind2].u,optr->vertex_table[ind2].v);
                glColor3ub(colorv[0],colorv[1],colorv[2]);
                }
              else
                {            
                glColor3ub(r,g,b);
                }
           /*glColor3ub(optr->vertex_table[ind2].rgb[0],
                      optr->vertex_table[ind2].rgb[1],
                      optr->vertex_table[ind2].rgb[2]);*/
           }
    glVertex3f(p3ptr->x, p3ptr->y, p3ptr->z );
    glEnd(); 


if (lineak == 0) { //Empty triangles

    //TODO draw the lines of the polygon. Ertzak marraztu 
    
    // 1-2 ertza
    drawLine(pgoiptr->x,pgoiptr->y,perdiptr->x,perdiptr->y);

    // 1-3 ertza
    drawLine(pgoiptr->x,pgoiptr->y,pbeheptr->x,pbeheptr->y);

    // 2-3 ertza
    drawLine(perdiptr->x,perdiptr->y,pbeheptr->x,pbeheptr->y);

    return;
}

// Segmentuz-segmentu marratzuko dut: 
    
    // Reset the intersection table to remove all data from the previous triangle
    initTable(); 

// TODO draw the segments of the triangle.
    
    // 1-2 ertza
    drawLineWT(pgoiptr->x,pgoiptr->y,optr->vertex_table[indg].u,optr->vertex_table[indg].v,perdiptr->x,perdiptr->y,optr->vertex_table[inde].u,optr->vertex_table[inde].v,optr);
           
    // 1-3 ertza
    drawLineWT(pgoiptr->x,pgoiptr->y,optr->vertex_table[indg].u,optr->vertex_table[indg].v,pbeheptr->x,pbeheptr->y,optr->vertex_table[indb].u,optr->vertex_table[indb].v,optr);
    
    // 2-3 ertza
    drawLineWT(perdiptr->x,perdiptr->y,optr->vertex_table[inde].u,optr->vertex_table[inde].v,pbeheptr->x,pbeheptr->y,optr->vertex_table[indb].u,optr->vertex_table[indb].v,optr);

    sortTable(); // for ScanLine

// TODO draw segments from upper vertex until midle vertex- goiko eta erdikoaren arteko segmentuak
    drawInternalPoints(pgoiptr->x,pgoiptr->y,perdiptr->x,perdiptr->y,optr);

// TODO draw segments from the midle to the lower vertex.
    drawInternalPoints(perdiptr->x,perdiptr->y,pbeheptr->x,pbeheptr->y,optr);   
}



void dibujar_poligono(object3d *optr, int ti)   
{
int i, ind0, ind1, ind2;
int atzeaurpegiada;
double * Nkam;

if (ti >= optr->num_faces) return;
// lehenengo hiru erpinekin kalkulatuko dut ikusgaitasuna.
ind0 = optr->face_table[ti].vertex_ind_table[0];
ind1 = optr->face_table[ti].vertex_ind_table[1];
ind2 = optr->face_table[ti].vertex_ind_table[2];

// TODO erabaki marraztu behar den ala ez: ikuste bolumenetik kanpokoak ez marraztu

// TODO poligonoaren kolorea zein da? oraingoz zuria.
optr->face_table[ti].rgb[0] = 255;
optr->face_table[ti].rgb[1] = 255;
optr->face_table[ti].rgb[2] = 255;
atzeaurpegiada = 0; //Initialize as not back-facing (default)
// TODO atze-aurpegia?
// atzeaurpegiada = ...
if ((!atzearpegiakmarraztu)&&atzeaurpegiada)
    {
    // Back culling...
    return;
    }
if (optr->texturaduna == 0) 
        {
        // TODO erpin bakoitzaren kolorea kalkulatu: argien, kameraren eta objektuaren orientazioaren arabera.
        argien_kalkulua_egin(optr, ti);
        }
// honaino iritsi bada bere triangelu guztiak marraztu behar ditut
for (i = 1; i<(optr->face_table[ti].num_vertices-1); i++)  // triangeluka marraztu: 4 erpinekin bi triangelu, bostekin 3...
    dibujar_triangulo(optr, ti, i, atzeaurpegiada);
}




static void marraztu(void)
{
float u,v;
int i,j;
object3d *auxptr;
double Fokudir[3];


  // marrazteko objektuak behar dira
  // no se puede dibujar sin objetos
if (foptr ==0) return;

// clear viewport...
if (objektuak == 1) glClear( GL_DEPTH_BUFFER_BIT|GL_COLOR_BUFFER_BIT );
    else 
      {
      if (denak == 0) glClear( GL_DEPTH_BUFFER_BIT|GL_COLOR_BUFFER_BIT );
      }
//TODO Ikuslearen edo kameraren erreferentzia-sistemara pasatzen duen matrizea lortu

//TODO argiek kameraren ikuspegian duten informazioa eguneratu (bai kokapenak, bai direkzioak)

if (objektuak == 1)
    {
    if (denak == 1)
        {
        //printf("objektuak marraztera\n");
        for (auxptr =foptr; auxptr != 0; auxptr = auxptr->hptr)
            {
            // TODO objektua kamerak nola ikusten duen adierazi objektuaren egituan bertan.
            kam_ikuspegia_lortu(auxptr);
            //printf("objektua kameraren ikuspegian daukat\n");
            for (i =0; i < auxptr->num_faces; i++)
                {
                dibujar_poligono(auxptr,i);
                }
            }
        }
      else
        {
        // TODO objektua kamerak nola ikusten duen adierazi objektuaren egituan bertan.
        kam_ikuspegia_lortu(sel_ptr);
        for (i =0; i < sel_ptr->num_faces; i++)
            {
            dibujar_poligono(sel_ptr,i);
            }
        }
    }
  else
    {
    // TODO objektua kamerak nola ikusten duen adierazi objektuaren egituan bertan.
    kam_ikuspegia_lortu(sel_ptr);
    dibujar_poligono(sel_ptr,indexx);
    }
glFlush();
}

void obj_normalak_kalkulatu(object3d *optr)
{
}


void read_from_file(char *fitx, object3d **fptrptr)
{
int i,retval;
object3d *optr;

    printf("%s fitxategitik datuak hartzera\n",fitx);
    optr = (object3d *)malloc(sizeof(object3d));
    optr->rgb.r = 0;
    //int read_wavefront(char * file_name, object3d * object_ptr) {
    retval = read_wavefront(fitx, optr);
    if (retval != 0) 
         {
         printf("%s fitxategitik datuak hartzerakoan arazoak izan ditut\n",fitxiz);
         free(optr);
         }
       else
         {
        
         //printf("objektuaren matrizea...\n");
         //Inicializar la Matrix de Transformacion de cada objeto
         optr->mptr = (mlist *)malloc(sizeof(mlist));
         for (i=0; i<16; i++) optr->mptr->m[i] =0;
         optr->mptr->m[0] = 1.0;
         optr->mptr->m[5] = 1.0;
         optr->mptr->m[10] = 1.0;
         optr->mptr->m[15] = 1.0;
         optr->mptr->hptr = 0;
        
         //printf("objektu edo kamera zerrendara doa informazioa...\n");
         optr->hptr = *fptrptr;
         *fptrptr = optr;
         //printf("normalak kalkulatzera\n");
         //for (i = 0; i< optr->num_vertices; i++) printf("%lf %lf %lf\n", optr->vertex_table[i].coord.x,optr->vertex_table[i].coord.y,optr->vertex_table[i].coord.z);
         obj_normalak_kalkulatu(optr);
         sel_ptr = optr;
         if (optr->texturaduna && (bufferra ==0))
            {
            // we put the information of the texture in the buffer pointed by bufferra. The dimensions of the texture are loaded into dimx and dimy
            retval = load_ppm("testura.ppm", &bufferra, &dimx, &dimy);
            if (retval !=1) 
                {
                printf("Ez dago testuraren fitxategia (testura.ppm)\n");
                optr->texturaduna=0;
                //exit(-1);
                }
            }
         }
     //printf("datuak irakurrita\n");
}




void x_aldaketa(int dir)
{

if (kamera == 0)// objektua aldatzen ari naiz
    {
    if (aldaketa =='r')
        {
        // rotate cos(5) = 0.99619469809174;  // cos(5)
        }
      else
        {
        // traslate: 0.02?
        }
    }
  else 
    if (kamera == 1) // kamera aldatzen ari naiz
      {
      if (ald_lokala == 1)// hegaldi moduan ezker/eskuin begiratu (y-rekiko biraketa)
        {
        // TODO hegan egin
        }
      else // analisi moduan
        {
        if (sel_ptr !=0)
            {
             // aukeratutako objektua analizatu eskuinerago edo ezkerreragotik (biratuz!)
             }
        }
      }
    else   // argiak aldatzen
      {
      // eguzkia biratu edo bonbila mugitu
      }
}



void y_aldaketa(int dir)
{

if (kamera == 0)// objektua aldatzen ari naiz
    {
    if (aldaketa =='r')
        {
        // rotate cos(5) = 0.99619469809174;  // cos(5)
        }
      else
        {
        // traslate: 0.02?
        }
    }
  else 
    if (kamera == 1) // kamera aldatzen ari naiz
      {
      if (ald_lokala == 1)// hegaldi moduan ezker/eskuin begiratu (y-rekiko biraketa)
        {
        // TODO hegan egin
        }
      else // analisi moduan
        {
        if (sel_ptr !=0)
            {
             // aukeratutako objektua analizatu goragotik edo beheragotik (biratuz!)
             }
        }
      }
    else   // argiak aldatzen
      {
      // eguzkia biratu edo bonbila mugitu
      }
}




void z_aldaketa(int dir)
{
if (kamera == 0) // objektuari aldaketa
    {
    if (aldaketa =='r')
        {
        // rotate cos(5) = 0.99619469809174;  // cos(5) 
        }
      else
        {
        // translate
        } 
    }
  else   // kamerari aldaketa
    if (kamera == 1) 
      {
      // hegaldi moduan beti aurrera edo atzera mugitu kamera.
      // analisi moduan traslazioa egin nahi bada objektura gerturatu (pasa gabe!! distantzia kontrolatu) edo urrutiratu
      // analisi moduan biraketa (roll)
      }
    else   // argiak aldatzen
      {
      // bonbilla mugitu munduan edo eguzkia biratu?
      }
}


void undo()
{
}

void kamera_objektuari_begira()
{
}


void print_egoera()
{
if (kamera == 0)
    {
    if (ald_lokala == 1) printf("\nobjektua aldatzen ari zara, (aldaketa lokala)\n");
        else printf("\nobjektua aldatzen ari zara, (aldaketa globala)\n");
    }
if (kamera == 1) 
    {
    if (ald_lokala == 1) printf("\nkamera aldatzen ari zara hegaldi moduan\n");
        else printf("\nkamera aldatzen ari zara analisi-moduan\n");
    }
if (kamera == 2)
    {
    printf("\nargiak aldatzen ari zara\n");
    }
if (aldaketa=='t') printf("Traslazioa dago aktibatuta\n");
    else printf("Biraketak daude aktibatuta\n");
if (objektuaren_ikuspegia) printf("objektuaren ikuspuntua erakusten ari zara (`C` sakatu kamerarenera pasatzeko)\n");
}

// This function will be called whenever the user pushes one key
static void teklatua (unsigned char key, int x, int y)
{
int retval;
int i;
FILE *obj_file;

switch(key)
	{
	case 13: 
	        if (foptr != 0)  // objekturik ez badago ezer ez du egin behar
	                         // si no hay objeto que no haga nada
	            {
	            indexx ++;  // azkena bada lehenengoa bihurtu
		                // pero si es el último? hay que controlarlo!
		    if (indexx == sel_ptr->num_faces) 
		        {
		        indexx = 0;
		        if ((denak == 1) && (objektuak == 0))
		            {
		            glClear( GL_DEPTH_BUFFER_BIT|GL_COLOR_BUFFER_BIT );
		            glFlush();
		            }
		        }
		    }
		break;
	case 'd':
		if (denak == 1) denak = 0;
		    else denak = 1;
		break;
	case 'o':
		if (objektuak == 1) objektuak = 0;
		    else objektuak = 1;
		break;
	case 'c':
	        if (kamera ==0) // objektua aldatzen ari naiz.
	            {
	            if (objektuaren_ikuspegia == 1) // objektuaren ikuspegian banago argiak aldatzera pasa naiteke, ez kamera aldatzera.
	                {
	                kamera =2;
	                printf("argiak aldatzera zoaz (objektuaren ikuspegian zaude)\n");
	                }
		        else 
		            {
		            kamera = 1;
		            ald_lokala = 1;  // hegaldi moduan jarriko naiz
	                 printf("kamera aldatzera zoaz (hegaldi moduan zaude)\n");
		            }
	            }
	          else
		    {
		    if (kamera == 1) 
		            {
		            kamera = 2;
	                printf("argia aldatzera zoaz \n");
	                }
		        else 
		            {
		            kamera = 0;
		            ald_lokala = 1;  // hegaldi moduan jarriko naiz
	                    printf("objektua aldatzera zoaz (aldaketa lokala daukazu) \n");
		            }
		    } 
		break;
	case 'C':
		if (objektuaren_ikuspegia == 1) objektuaren_ikuspegia = 0;
		    else 
		        {
		        objektuaren_ikuspegia = 1;
		        if (kamera == 1) 
		            {
		            kamera=0;
	                }
		        ald_lokala = 1;
		        printf("objektuaren ikuspegian sartu zara eta aldaketak era lokalean eragingo dituzu \n");
		        }
		break;
	case 'l':
		if (lineak == 1) lineak = 0;
		    else lineak = 1;
		break;
	case 't':
	        aldaketa = 't';
	        printf("traslazioak\n");
		break;
	case 'r':
		aldaketa = 'r';
		printf("biraketak\n");
		break;
	case 'n':
		if (atzearpegiakmarraztu == 1) 
		    {
		    atzearpegiakmarraztu = 0;
		    printf("atze-aurpegiak ez dira marraztuko\n");
		    }
		  else 
		    {
		    atzearpegiakmarraztu = 1;
		    printf("atze aurpegiak gorriz marraztuko dira\n");
		    }
		break;
	case 'g':
	        if (objektuaren_ikuspegia == 0) // objektuaren ikuspegian aldaketa beti lokala izango da.
	            {
	            if (ald_lokala == 1)
	                    {
	                    ald_lokala = 0;
	                    if ((kamera==1) &&(sel_ptr!=0)) 
	                        {
	                        kamera_objektuari_begira();
	                        //print_matrizea16(M_kam);
	                        }
	                    }
	                else ald_lokala = 1;
	            }
	          else printf("objektuaren ikuspegian aldaketak beti lokalak dira\n");
		break;
    case 'x':
                x_aldaketa(1);
                break;
    case 'y':
                y_aldaketa(1);
                break;
    case 'z':
                z_aldaketa(1);
                break;
    case 'X':
                x_aldaketa(0);
                break;
    case 'Y':
                y_aldaketa(0);
                break;
    case 'Z':
                z_aldaketa(0);
                break;
    case 'u':
                undo();
                break;
	case 'p':
		if (paralelo == 1) 
		    {
		    printf("perspektiban agertu behar du\n");
		    paralelo = 0;
		    }
		    else
		    {
		    printf("paraleloan agertu behar du\n");
		     paralelo = 1;
		    }
		break;
	case '1':
		// Eguzkia piztu/itzali
		break;
	case '2':
		// Bonbilla piztu/itzali
		break;
	case '3':
		// objektuaren fokua piztu/itzali
		break;
	case '4':
		// kameraren fokua piztu/itzali 
		break;
	case '+':
		if (kamera == 2) 
		    {// Fokuen irekiera handitu
		    }
		  else 
		    {
		    if (kamera == 1)
		        {// kameraren ikuste-bolumena handitu
		        }
		    }
		break;
	case '-':
		if (kamera == 2) 
		    {// Fokuen irekiera txikitu
		    }
		  else 
		    {
		    if (kamera == 1)
		        {// kameraren ikuste-bolumena txikitu
		        }
		    }
		break;
	case 'f':
	        /*Ask for file*/
	        printf("idatzi fitxategi izena\n"); 
	        scanf("%s", &(fitxiz[0]));
	        read_from_file(fitxiz,&foptr);
	        indexx = 0;
	        if ((ald_lokala==0)&&(kamera==1) &&(sel_ptr!=0)) kamera_objektuari_begira();
                break;
    case 9: /* <TAB> */
            if (foptr != 0) // objekturik gabe ez du ezer egin behar
                            // si no hay objeto no hace nada
                {
                sel_ptr = sel_ptr->hptr;
                /*The selection is circular, thus if we move out of the list we go back to the first element*/
                if (sel_ptr == 0) sel_ptr = foptr;
                indexx =0; // the selected polygon is the first one
                if ((ald_lokala == 0)&&(kamera==1))
                    {
                    //kamera objektuari begira jarri behar da!!
                    kamera_objektuari_begira();
                    }
                }
            break;
	case 27:  // <ESC>
		exit( 0 );
		break;
	default:
		printf("%d %c\n", key, key );
	}
print_egoera();
// The screen must be drawn to show the new triangle
glutPostRedisplay();
}


void viewportberria (int zabal, int garai)
{
    if (zabal < garai)  dimentsioa = zabal;
    else  dimentsioa = garai;

    // Free the previous table if it existed    
    if (intersectionTable != NULL) free(intersectionTable);

    // Resize
    intersectionTable = (ScanLine *) malloc(sizeof(ScanLine) * dimentsioa);

    glViewport(0,0,dimentsioa,dimentsioa);
    printf("linea kopuru berria = %d\n",dimentsioa);
}

int main(int argc, char** argv) //argc: ARGument Count - argv: ARGument Vector 
{
int retval,i;

	printf(" Triangeluak: barneko puntuak eta testura\n Triángulos con puntos internos y textura \n");
	printf("Press <ESC> to finish\n");
	glutInit(&argc,argv);
	glutInitDisplayMode ( GLUT_RGB|GLUT_DEPTH );
	dimentsioa = 500;
	glutInitWindowSize ( dimentsioa, dimentsioa );
	glutInitWindowPosition ( 100, 100 );
	glutCreateWindow( "KBG/GC praktika" );

	glutDisplayFunc( marraztu );
	glutKeyboardFunc( teklatua );
	glutReshapeFunc( viewportberria);
	//  we put the information of the texture in the buffer pointed by bufferra. The dimensions of the texture are loaded into dimx and dimy
        retval = load_ppm("testura.ppm", &bufferra, &dimx, &dimy);
        //printf("Textura cargada: %d x %d\n", dimx, dimy);
        if (retval !=1) 
            {
            printf("Ez dago testuraren fitxategia (testura.ppm)\n");
            exit(-1);
            }
    //placeholder (Yo)         
    //bufferra =0; dimx=0;dimy=0;
	glClearColor( 0.0f, 0.0f, 0.7f, 1.0f );
	glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
        glEnable(GL_DEPTH_TEST); // activar el test de profundidad (Z-buffer)
        glDepthFunc(GL_GREATER);  glClearDepth(0.0);// Handiena marraztu
        //glMatrixMode(GL_PROJECTION);
        //glLoadIdentity();
        //glOrtho(-1.0, 1.0, -1.0, 1.0, 1.0, -1.0);
        //glMatrixMode(GL_MODELVIEW);
        denak = 0;
        lineak =0;
        objektuak = 1;
        kamera = 0;
        foptr = 0;
        sel_ptr = 0;
        aldaketa = 'r';
        ald_lokala =1;
        objektuaren_ikuspegia =0;
        paralelo = 1;
        atzearpegiakmarraztu = 1;
        // TODO kamera hasieratu kamera (0,0,2.5) kokapenean dago hasieran
        
        // TODO Argiak hasieratu
        
        if (argc>1) read_from_file(argv[1],&foptr);
        else 
        {
            /*
            read_from_file("abioia-1+1.obj",&foptr);
            if (sel_ptr != 0) 
                { //sel_ptr->mptr->m[3] = -1.0;
                }   
            read_from_file("abioia-1+1.obj",&foptr);
            if (sel_ptr != 0) 
                { //sel_ptr->mptr->m[3] = 1.0;
                }   
            read_from_file("abioia-1+1.obj",&foptr);
            if (sel_ptr != 0) 
                { 
                //sel_ptr->mptr->m[7] = -0.4;
                //sel_ptr->mptr->m[11] = 0.5;
                }   
            read_from_file("abioia-1+1.obj",&foptr);
            if (sel_ptr != 0) 
                { //sel_ptr->mptr->m[11] = -0.7;
                }        
            read_from_file("abioia-1+1.obj",&foptr);
            if (sel_ptr != 0) 
                { //sel_ptr->mptr->m[7] = 0.7;
                }
            */

            read_from_file("z-1+1.obj",&foptr);
            read_from_file("abioia-1+1.obj",&foptr);
            read_from_file("triangles-1+1.obj",&foptr);
        }
        print_egoera();
	glutMainLoop();

	return 0;   
}