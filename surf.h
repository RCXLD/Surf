#ifndef SURF_H
#define SURF_H


class SURF
{
    public:
        SURF();
        virtual ~SURF();
        void eval(unsigned char *in, float *out, int xres, int yres,float scale,float orientation);
    private:
        int xres, yres;
        unsigned char *img;

        int indexof(int x,int y);
        inline int fRound(float flt);
        inline float gaussian(float x, float y, float sig);
        inline float haarX(int row, int column, int s);
        inline float haarY(int row, int column, int s);
        inline float BoxIntegral(unsigned char *data, int row, int col, int rows, int cols);
};

#endif // SURF_H
