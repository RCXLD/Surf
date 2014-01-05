#include "surf.h"
#include <cmath>
#include <bits/stl_algobase.h>  //min
SURF::SURF()
{
}

SURF::~SURF()
{
}

#include <iostream>
void SURF::eval(unsigned char *in, float *out, int xres, int yres,float scale,float orientation)
{
    img = in;
    float co,si;

    this->xres = xres;
    this->yres = yres;
    co = cos(orientation);
    si = sin(orientation);


    #pragma omp parallel for
    for(int y = 0; y < yres; y ++)
    {
        int sample_x, sample_y, count;
        int i, ix, j, jx, xs, ys;
        #ifdef SURF_128
        float dx1, dx2, dx3, dx4, dy1, dy2, dy3, dy4;
        #else
        float dx, dy, mdx, mdy;
        #endif
        float gauss_s1, gauss_s2;
        float rx, ry, rrx, rry, len;
        float cx,cy; //Subregion centers for the 4x4 gaussian weighting
        for(int x = 0; x < xres; x ++)
        {
            #ifdef SURF_128
            float *desc = out+(x+y*xres)*128;
            #else
            float *desc = out+(x+y*xres)*64;
            #endif
            count = 0;
            len = 0;
            cx = -0.5f;
            i = -8;
            //Calculate descriptor
            while(i < 12)
            {
                j = -8;
                i = i-4;
                cx += 1.f;
                cy = -0.5f;
                while(j < 12)
                {
                    #ifdef SURF_128
                    dx1 = dx2 = dx3 = dx4 = dy1 = dy2 = dy3 = dy4 = 0.f;
                    #else
                    dx=dy=mdx=mdy=0.f;
                    #endif
                    cy += 1.f;
                    j = j - 4;
                    ix = i + 5;
                    jx = j + 5;
                    xs = fRound(x + ( -jx*scale*si + ix*scale*co));
                    ys = fRound(y + ( jx*scale*co + ix*scale*si));

                    for (int k = i; k < i + 9; ++k)
                    {
                        for (int l = j; l < j + 9; ++l)
                        {
                            //Get coords of sample point on the rotated axis
                            sample_x = fRound(x + (-l*scale*si + k*scale*co));
                            sample_y = fRound(y + ( l*scale*co + k*scale*si));

                            //Get the gaussian weighted x and y responses
                            gauss_s1 = gaussian(xs-sample_x,ys-sample_y,2.5f*scale);
                            rx = haarX(sample_y, sample_x, 2*fRound(scale));
                            ry = haarY(sample_y, sample_x, 2*fRound(scale));

                            //Get the gaussian weighted x and y responses on rotated axis
                            rrx = gauss_s1*(-rx*si + ry*co);
                            rry = gauss_s1*(rx*co + ry*si);

                            #ifdef SURF_128
                            if(rry > 0)
                                if(rrx > 0)
                                {
                                    dx1 += rrx;
                                    dy1 += rry;
                                }
                                else
                                {
                                    dx2 -= rrx;
                                    dy2 += rry;
                                }
                            else
                                if(rrx < 0)
                                {
                                    dx3 -= rrx;
                                    dy3 -= rry;
                                }
                                else
                                {
                                    dx4 += rrx;
                                    dy4 -= rry;
                                }
                            #else
                            dx += rrx;
                            dy += rry;
                            mdx += fabs(rrx);
                            mdy += fabs(rry);
                            #endif
                        }
                    }

                    //Add the values to the descriptor vector
                    gauss_s2 = gaussian(cx-2.0f,cy-2.0f,1.5f);

                    #ifdef SURF_128
                    desc[count++] = dx1*gauss_s2;
                    desc[count++] = dy1*gauss_s2;
                    desc[count++] = dx2*gauss_s2;
                    desc[count++] = dy2*gauss_s2;
                    desc[count++] = dx3*gauss_s2;
                    desc[count++] = dy3*gauss_s2;
                    desc[count++] = dx4*gauss_s2;
                    desc[count++] = dy4*gauss_s2;
                    len += (dx1*dx1 +
                            dy1*dy1 +
                            dx2*dx2 +
                            dy2*dy2 +
                            dx3*dx3 +
                            dy3*dy3 +
                            dx4*dx4 +
                            dy4*dy4 ) * gauss_s2*gauss_s2;
                    #else
                    desc[count++] = dx*gauss_s2;
                    desc[count++] = dy*gauss_s2;
                    desc[count++] = mdx*gauss_s2;
                    desc[count++] = mdy*gauss_s2;
                    len += (dx*dx + dy*dy + mdx*mdx + mdy*mdy) * gauss_s2*gauss_s2;
                    #endif

                    j += 9;
                }
                i += 9;
            }
            //Convert to Unit Vector
            len = sqrt(len);
            #ifdef SURF_128
            for(int i = 0; i < 128; ++i)
            #else
            for(int i = 0; i < 64; ++i)
            #endif
                desc[i] /= len;
        }
    }
}

int SURF::indexof(int x,int y)
{
    if(x < 0)
        x = 0;
    else if(x >= xres)
        x = xres-1;
    if(y < 0)
        y = 0;
    else if(y >= yres)
        y = yres-1;
    return x+y*xres;
}

inline int SURF::fRound(float flt)
{
  return (int) floor(flt+0.5f);
}

inline float SURF::gaussian(float x, float y, float sig)
{
  return 1.0f/(2.0f*M_PI*sig*sig) * exp( -(x*x+y*y)/(2.0f*sig*sig));
}


inline float SURF::BoxIntegral(unsigned char *data, int row, int col, int rows, int cols)
{

  // The subtraction by one for row/col is because row/col is inclusive.
  int r1 = std::min(row,yres) - 1;
  int c1 = std::min(col,xres) - 1;
  int r2 = std::min(row + rows, yres) - 1;
  int c2 = std::min(col + cols, xres)  - 1;

  float A(0.0f), B(0.0f), C(0.0f), D(0.0f);
  if (r1 >= 0 && c1 >= 0) A = data[r1 * xres + c1];
  if (r1 >= 0 && c2 >= 0) B = data[r1 * xres + c2];
  if (r2 >= 0 && c1 >= 0) C = data[r2 * xres + c1];
  if (r2 >= 0 && c2 >= 0) D = data[r2 * xres + c2];

  return std::max(0.f, (float)A - B - C + D);
}

inline float SURF::haarX(int row, int column, int s)
{
  return BoxIntegral(img, row-s/2, column, s, s/2)
    -1 * BoxIntegral(img, row-s/2, column-s/2, s, s/2);
}

inline float SURF::haarY(int row, int column, int s)
{
  return BoxIntegral(img, row, column-s/2, s/2, s)
    -1 * BoxIntegral(img, row-s/2, column-s/2, s/2, s);
}
