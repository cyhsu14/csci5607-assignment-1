#include "image.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <random>
#include <iostream>
#include <algorithm>
#include <list>
using namespace std;
/**
 * Image
 **/
Image::Image (int width_, int height_){

    assert(width_ > 0);
    assert(height_ > 0);

    width           = width_;
    height          = height_;
    num_pixels      = width * height;
    sampling_method = IMAGE_SAMPLING_POINT;

    data.raw = new uint8_t[num_pixels*4];
	int b = 0; //which byte to write to
	for (int j = 0; j < height; j++){
		for (int i = 0; i < width; i++){
			data.raw[b++] = 0;
			data.raw[b++] = 0;
			data.raw[b++] = 0;
			data.raw[b++] = 0;
		}
	}

    assert(data.raw != NULL);
}

Image::Image (const Image& src){

	width           = src.width;
    height          = src.height;
    num_pixels      = width * height;
    sampling_method = IMAGE_SAMPLING_POINT;

    data.raw = new uint8_t[num_pixels*4];

    //memcpy(data.raw, src.data.raw, num_pixels);
    *data.raw = *src.data.raw;
}

Image::Image (char* fname){

	int numComponents; //(e.g., Y, YA, RGB, or RGBA)
	data.raw = stbi_load(fname, &width, &height, &numComponents, 4);

	if (data.raw == NULL){
		printf("Error loading image: %s", fname);
		exit(-1);
	}


	num_pixels = width * height;
	sampling_method = IMAGE_SAMPLING_POINT;

}

Image::~Image (){
    delete data.raw;
    data.raw = NULL;
}

void Image::Write(char* fname){

	int lastc = strlen(fname);
	switch (fname[lastc-1]){
	   case 'g': //jpeg (or jpg) or png
	     if (fname[lastc-2] == 'p' || fname[lastc-2] == 'e') //jpeg or jpg
	        stbi_write_jpg(fname, width, height, 4, data.raw, 95);  //95% jpeg quality
	     else //png
	        stbi_write_png(fname, width, height, 4, data.raw, width*4);
	     break;
	   case 'a': //tga (targa)
	     stbi_write_tga(fname, width, height, 4, data.raw);
	     break;
	   case 'p': //bmp
	   default:
	     stbi_write_bmp(fname, width, height, 4, data.raw);
	}
}

void Image::AddNoise (double factor)
{
    int x,y;
    default_random_engine gen;

    normal_distribution<double> d(0,15);
	for (x = 0 ; x < Width() ; x++)
	{
		for (y = 0 ; y < Height() ; y++)
		{
			Pixel p = GetPixel(x, y);
            double ran = d(gen);
            // printf("%f %f\n",g, factor*(g));
            double r,g,b;
            Pixel noised_p;
            r = (int)p.r + factor*ran;
            g = (int)p.g + factor*ran;
            b = (int)p.b + factor*ran;
            if(r<0) r=0;
            if(g<0) g=0;
            if(b<0) b=0;
            if(r>255)r=255;
            if(g>255)g=255;
            if(b>255)b=255;
            // noised_p = p*factor*g;
            noised_p.r = (int)r;
            noised_p.g = (int)g;
            noised_p.b = (int)b;
            // printf("%d %d %d\n",noised_p.r, noised_p.g,noised_p.b);
            // noised_p = p*factor*g;
            noised_p.SetClamp(noised_p.r,noised_p.g,noised_p.b);

			GetPixel(x,y) = noised_p;
		}
	}

}

void Image::Brighten (double factor)
{
	int x,y;
	for (x = 0 ; x < Width() ; x++)
	{
		for (y = 0 ; y < Height() ; y++)
		{
			Pixel p = GetPixel(x, y);
            Pixel scaled_p;
			scaled_p.r = factor*p.r;
			scaled_p.g = factor*p.g;
			scaled_p.b = factor*p.b;
			GetPixel(x,y) = scaled_p;
		}
	}
}


void Image::ChangeContrast (double factor)
{
    int x,y;
    unsigned int avg_gray = 0;
	for (x = 0 ; x < Width() ; x++)
	{
		for (y = 0 ; y < Height() ; y++)
		{
			Pixel p = GetPixel(x, y);
			avg_gray += p.Luminance();

		}
	}
    avg_gray = avg_gray/(Width()*Height());
    for (x = 0 ; x < Width() ; x++)
	{
		for (y = 0 ; y < Height() ; y++)
		{
			Pixel p = GetPixel(x, y);


            p.r = ComponentLerp(avg_gray, p.r,factor);
            p.g = ComponentLerp(avg_gray, p.g,factor);
            p.b = ComponentLerp(avg_gray, p.b,factor);

            GetPixel(x,y) = p;
		}
	}


}


void Image::ChangeSaturation(double factor)
{
	/* WORK HERE */
    int x,y;
    for (x = 0 ; x < Width() ; x++)
	{
		for (y = 0 ; y < Height() ; y++)
		{
			Pixel p = GetPixel(x, y);
            // if factor > 1
            int lu = p.Luminance();
			p.r = ComponentLerp(lu, p.r, factor);
            p.g = ComponentLerp(lu, p.g, factor);
            p.b = ComponentLerp(lu, p.b, factor);

            GetPixel(x,y) = p;
		}
	}
}


Image* Image::Crop(int x, int y, int w, int h)
{
	/* WORK HERE */
    Image *result = new Image(w,h);

    int i,j;
    int a=0, b=0;
	for (i = x-w/2 ; i < x+w/2 ; i++)
	{
		for (j = y-h/2 ; j < y+h/2 ; j++)
		{
            result->data.pixels[a+b*h] = GetPixel(i,j);
            b++;
		}
        a++;
        b=0;
	}
	return result;
}


void Image::ExtractChannel(int channel)
{
	/* WORK HERE */
    // 0-> r, 1->g, 2->b
    if (channel>2) {
        printf("Wrong channel!\n");
        return;
    }
    int x,y;
	for (x = 0 ; x < Width() ; x++)
	{
		for (y = 0 ; y < Height() ; y++)
		{
			Pixel p = GetPixel(x, y);
            if(channel == 0) {
                p.g = 0;
                p.b = 0;
            }
            else if(channel == 1) {
                p.r = 0;
                p.b = 0;
            }
            else if(channel == 2) {
                p.r = 0;
                p.g = 0;
            }

			GetPixel(x,y) = p;
		}
	}


}

void Image::Quantize (int nbits)
{
	/* WORK HERE */
    if (nbits >8 || nbits<=0) {
        printf("Invalid bits!\n");
        return;
    }
    int x,y;
    Pixel p;
    int margin = 256 / pow(2,nbits);

    for (x = 0 ; x < Width() ; x++)
	{
		for (y = 0 ; y < Height() ; y++)
		{
			p = GetPixel(x, y);
            int tmp_r, tmp_g, tmp_b;
            for(int i = 0; i < pow(2,nbits); i++){
                if((int)p.r > i*margin) {
                    if(i == pow(2,nbits)-1) tmp_r = 255;
                    else tmp_r = i*margin;
                }
                if((int)p.g > i*margin) {
                    if(i == pow(2,nbits)-1) tmp_g = 255;
                    else tmp_g = i*margin;
                }
                if((int)p.b > i*margin) {
                    if(i == pow(2,nbits)-1) tmp_b = 255;
                    else tmp_b = i*margin;
                }

            }
            // if(x==250 && y==500)printf("%d %d %d \n", p.r,p.g,p.b);
            p.r = (unsigned char)tmp_r;
            p.g = (unsigned char)tmp_g;
            p.b = (unsigned char)tmp_b;
            p.SetClamp(p.r,p.g,p.b);

            // if(x==250 && y==500)printf("%d %d %d \n", p.r,p.g,p.b);
            GetPixel(x,y) = p;

		}
	}


}

void Image::RandomDither (int nbits)
{
	/* WORK HERE */
    int n = pow(2,nbits);
    int x, y, k;
    int ran = 2;
    AddNoise(1.5);
    Quantize(nbits);

}


static int Bayer4[4][4] =
{
    {15,  7, 13,  5},
    { 3, 11,  1,  9},
    {12,  4, 14,  6},
    { 0,  8,  2, 10}
};


void Image::OrderedDither(int nbits)
{
	/* WORK HERE */
}

/* Error-diffusion parameters */
const double
    ALPHA = 7.0 / 16.0,
    BETA  = 3.0 / 16.0,
    GAMMA = 5.0 / 16.0,
    DELTA = 1.0 / 16.0;

void Image::FloydSteinbergDither(int nbits)
{
	/* WORK HERE */
    int n = pow(2,nbits);
    int x, y;
    int margin = 256 / pow(2,nbits);
    for (y = 0 ; y < Height() ; y++)
	{
		for (x = 0 ; x < Width() ; x++)
		{
			Pixel p = GetPixel(x, y);
            int tmp_r, tmp_g, tmp_b;
            for(int i = 0; i < pow(2,nbits); i++){
                if((int)p.r > i*margin) {
                    if(i == pow(2,nbits)-1) tmp_r = 255;
                    else tmp_r = i*margin;
                }
                if((int)p.g > i*margin) {
                    if(i == pow(2,nbits)-1) tmp_g = 255;
                    else tmp_g = i*margin;
                }
                if((int)p.b > i*margin) {
                    if(i == pow(2,nbits)-1) tmp_b = 255;
                    else tmp_b = i*margin;
                }

            }

            // pass difference to neighbors
            double diff_r, diff_g, diff_b;
            diff_r = (int)p.r - tmp_r;
            diff_g = (int)p.g - tmp_g;
            diff_b = (int)p.b - tmp_b;
            Pixel z;
            double r=0,g=0,b=0;
            if(x!=Width()-1){
                z = GetPixel(x+1,y);
                // r = z.r+diff_r*ALPHA;
                // g = z.g+diff_g*ALPHA;
                // b = z.b+diff_b*ALPHA;
                // z.SetClamp(r,g,b);
                z = Pixel((int)(z.r+diff_r*ALPHA), (int)(z.g+diff_g*ALPHA),(int)(z.b+diff_b*ALPHA),z.a);
                z.SetClamp(z.r,z.g,z.b);
                GetPixel(x+1,y) = z;
            }
            if(y!=Height()-1){
                if(x!=0){
                    z = GetPixel(x-1,y+1);
                    z = Pixel((int)(z.r+diff_r*BETA),(int)(z.g+diff_g*BETA),(int)(z.b+diff_b*BETA),z.a);
                    z.SetClamp(z.r,z.g,z.b);
                    GetPixel(x-1,y+1) = z;
                }
                z = GetPixel(x,y+1);
                z = Pixel((int)(z.r+diff_r*GAMMA),(int)(z.g+diff_g*GAMMA),(int)(z.b+diff_b*GAMMA),z.a);
                z.SetClamp(z.r,z.g,z.b);
                GetPixel(x,y+1) = z;
                if(x!=Width()-1){
                    z = GetPixel(x+1,y+1);
                    z = Pixel((int)(z.r+diff_r*DELTA),(int)(z.g+diff_g*DELTA),(int)(z.b+diff_b*DELTA),z.a);
                    z.SetClamp(z.r,z.g,z.b);
                    GetPixel(x+1,y+1) = z;
                }
            }
            p.SetClamp(tmp_r,tmp_g,tmp_b);
            GetPixel(x,y) = p;

		}
	}

}

void Image::Blur(int n)
{
	/* WORK HERE */
    double sigma = 5;
    double sum = 0.0;

    int x,y;
    double** kernel = new double*[n];
    for(x = 0;x < n; x++) kernel[x] = new double[n];
    // generate kernel
    for (x = 0; x < n; x++) {
        for (y = 0; y < n; y++) {
            kernel[x][y] = exp(-(x*x+y*y)/(2*sigma*sigma))/(2*M_PI*sigma*sigma);
            sum += kernel[x][y];
        }
    }
    // normalising
    for (x = 0; x < n; x++) {
        for (y = 0; y < n; y++) {
            kernel[x][y] /= sum;
        }
    }

    // convolute
    int i,j;

    for(x = 0;x < Width()-n+1; x++){
        for(y = 0;y < Height()-n+1; y++){
            // Pixel p = new_image->GetPixel(x, y);
            Pixel p = Pixel(0,0,0,0);
            for(i = x; i < x+n;i++){
                for(j = y; j < y+n; j++){
                    Pixel tmp = this->GetPixel(i, j);
                    p = p + tmp*kernel[(i-x)][(j-y)];
                }
            }

            GetPixel(x, y) = p;
        }
    }

}

void Image::Sharpen(int n)
{
	/* WORK HERE */
    int x,y;
    Image *before_sharp = new Image(*this);
    for(x = 0;x<Width();x++){
        for(y=0;y<Height();y++)
            before_sharp->GetPixel(x,y) = GetPixel(x,y);
    }
    this->Blur(n);

    double factor = 2;
    for(x = 0;x < Width(); x++){
        for(y = 0;y < Height(); y++){
            Pixel a = before_sharp->GetPixel(x,y);
            // printf("%d %d %d\n", a.r,a.g,a.b);
            Pixel b = this->GetPixel(x,y);
            b = PixelLerp(b,a,factor);
            b.SetClamp(b.r,b.g,b.b);
            // b = 2*a - b;

            this->GetPixel(x,y) = b;
        }
    }
}

void Image::EdgeDetect()
{
	/* WORK HERE */
    float kernel[3][3] = {-1,-1,-1,
                          -1, 8,-1,
                          -1,-1,-1};
    int x,y;
    int i,j;
    int n=3;
    // this->Quantize(4);
    for(x = 0;x < Width()-n+1; x++){
      for(y = 0;y < Height()-n+1; y++){
          Pixel p = Pixel(0,0,0,0);
          int r=0,g=0,b=0;
          for(i = x; i < x+n;i++){
              for(j = y; j < y+n; j++){
                  Pixel t = this->GetPixel(i, j);
                  r = r + (int)t.r*kernel[(i-x)][(j-y)];
                  g = g + (int)t.g*kernel[(i-x)][(j-y)];
                  b = b + (int)t.b*kernel[(i-x)][(j-y)];
              }
          }
          p.r = r;
          p.g = g;
          p.b = b;

          p.SetClamp(r,g,b);

          GetPixel(x, y) = p;
      }
    }


}

Image* Image::Scale(double sx, double sy)
{
	/* WORK HERE */
    int x,y;
    int sw = Width()*sx;
    int sh = Height()*sy;
    Image *scale = new Image(sw,sh);
    // IMAGE_SAMPLING_POINT
    Pixel p;
    for(x = 0; x < sw; x++){
        for(y = 0; y<sh; y++){
            p = Sample(x/sx,y/sy);
            scale->GetPixel(x,y) = p;
        }
    }
	return scale;
}

Image* Image::Rotate(double angle)
{
	/* WORK HERE */
    int x,y;

    // convert angle in radius
    double rad = angle*3.14/180;
    int rw = Height()*abs(sin(rad))+Width()*abs(cos(rad));
    int rh = Height()*abs(cos(rad))+Width()*abs(sin(rad));
    // int rw = Height()+ Width();
    // int rh = Height()+ Width();
    Image *rotate = new Image(rw,rh);
    // printf("%f %f\n",u,v );
    for(x = 0; x < rw; x++){
        for(y = 0; y < rh; y++){
            double tmp = x;
            double u = cos(rad)*tmp + sin(rad)*y;
            double v = -sin(rad)*tmp + cos(rad)*y;
            if(angle <= 90){ //v
                u = u-Height()*cos(rad)*sin(rad);
                v = v+Height()*sin(rad)*sin(rad);
            }
            else if(angle <=180){
                double ux = Width()*sin(rad-1.57)*sin(rad-1.57);
                double uv = Width()*cos(rad-1.57)*sin(rad-1.57)+Height();
                u = u+ux;
                v = v+uv;
            }
            else if(angle <=270){
                double ux = Height()*cos(rad-3.14)*sin(rad-3.14)+Width();
                double uv = Height()*cos(rad-3.14)*cos(rad-3.14);
                u = u+ux;
                v = v+uv;
            }
            else { //v
                u = u+Width()*cos(rad-4.71)*cos(rad-4.71);
                v = v-Width()*cos(rad-4.71)*sin(rad-4.71);
            }

            rotate->GetPixel(x,y) = (u<0 || v<0 || u>Width() || v>Height()) ? Pixel(255,255,255) : Sample(u,v);
        }
    }


	return rotate;
}

void Image::Fun()
{
	/* WORK HERE */
    int x, y;
    double r,a,rn;
    Pixel p;
    Image *fun = new Image(*this);
    for(x = 0; x < Width(); x++){
        for(y = 0; y<Height(); y++){
            // printf("%f %f\n",x/sx,y/sy);
            r = sqrt(pow((x - Width()/2),2) + pow((y - Height()/2),2));
            a = atan2( y-Height()/2,x-Width()/2);
            rn = 2*pow(r/sqrt(pow(Width(),2) + pow(Height(),2)),1.1)*r;
            double u = rn*(cos(a))+Width()/2;
            double v = rn*sin(a)+Height()/2;
            p = (u<0 || v<0 || u>Width() || v>Height()) ? Pixel(255,255,255) : Sample(u,v);
            fun->GetPixel(x,y) = p;
        }
    }
    for(x = 0; x < Width(); x++){
        for(y = 0; y<Height(); y++){
            p = fun->GetPixel(x,y);
            this->GetPixel(x,y) = p;
        }
    }

}

/**
 * Image Sample
 **/
void Image::SetSamplingMethod(int method)
{
    assert((method >= 0) && (method < IMAGE_N_SAMPLING_METHODS));
    sampling_method = method;
}


Pixel Image::Sample (double u, double v){
    /* WORK HERE */
    Pixel p;
    int x,y;
    int x0,y0,x1,y1;
    double r,g,b;
    if( sampling_method == IMAGE_SAMPLING_POINT){
        p = GetPixel((int)floor(u),(int)floor(v));
    }
    else if( sampling_method == IMAGE_SAMPLING_BILINEAR){
        x0 = floor(u);
        x1 = x0+1;
        if(x1==Width())x1-=1;
        y0 = floor(v);
        y1 = y0+1;
        if(y1==Height())y1-=1;
        // area: 0 1    (x0,y0) (x1,y0)
        //       2 3    (x0,y1) (x1,y1)
        double area0, area1,area2,area3;
        area0 = (u-x0)*(v-y0);
        area1 = (x1-u)*(v-y0);
        area2 = (u-x0)*(y1-v);
        area3 = (x1-u)*(y1-v);

        if(u!=x0 && x1!=u && y0!=v && v!=y1){
            r = GetPixel(x1,y1).r*area0+GetPixel(x0,y1).r*area1+GetPixel(x1,y0).r*area2+GetPixel(x0,y0).r*area3;
            g = GetPixel(x1,y1).g*area0+GetPixel(x0,y1).g*area1+GetPixel(x1,y0).g*area2+GetPixel(x0,y0).g*area3;
            b = GetPixel(x1,y1).b*area0+GetPixel(x0,y1).b*area1+GetPixel(x1,y0).b*area2+GetPixel(x0,y0).b*area3;
        }
        else{
            if(y0==v){
                r = GetPixel(x1,y0).r*(u-x0)+GetPixel(x0,y0).r*(x1-u);
                g = GetPixel(x1,y0).g*(u-x0)+GetPixel(x0,y0).g*(x1-u);
                b = GetPixel(x1,y0).b*(u-x0)+GetPixel(x0,y0).b*(x1-u);
            }
            else if(y1==v){
                r = GetPixel(x1,y1).r*(u-x0)+GetPixel(x0,y1).r*(x1-u);
                g = GetPixel(x1,y1).g*(u-x0)+GetPixel(x0,y1).g*(x1-u);
                b = GetPixel(x1,y1).b*(u-x0)+GetPixel(x0,y1).b*(x1-u);
            }
            else if(x0==u){
                r = GetPixel(x0,y1).r*(v-y0)+GetPixel(x0,y0).r*(y1-v);
                g = GetPixel(x0,y1).g*(v-y0)+GetPixel(x0,y0).g*(y1-v);
                b = GetPixel(x0,y1).b*(v-y0)+GetPixel(x0,y0).b*(y1-v);
            }
            else if(x1==u){
                r = GetPixel(x1,y1).r*(v-y0)+GetPixel(x1,y0).r*(y1-v);
                g = GetPixel(x1,y1).g*(v-y0)+GetPixel(x1,y0).g*(y1-v);
                b = GetPixel(x1,y1).b*(v-y0)+GetPixel(x1,y0).b*(y1-v);
            }
        }
        p.r = (int)(r+0.5);
        p.g = (int)(g+0.5);
        p.b = (int)(b+0.5);
        p.SetClamp(p.r,p.g,p.b);
        if((x0==u&&y0==v) || (x1==u&&y0==0) || (x1==u&&y1==v) || (x1==u&&y1==v))p = GetPixel((int)u,(int)v);

    }
    else if( sampling_method == IMAGE_SAMPLING_GAUSSIAN){
        double sigma = 1;
        double sum = 0.0;
        int n=4;
        double** kernel = new double*[n];
        for(x = 0;x < n; x++) kernel[x] = new double[n];
        // generate kernel
        for (x = 0; x < n; x++) {
            for (y = 0; y < n; y++) {
                kernel[x][y] = exp(-(x*x+y*y)/(2*sigma*sigma))/(2*M_PI*sigma*sigma);
                sum += kernel[x][y];
            }
        }
        // normalising
        for (x = 0; x < n; x++) {
            for (y = 0; y < n; y++) {
                kernel[x][y] /= sum;
            }
        }

        // multiply with gaussian kernel

        if(floor(u) == 0){ x0 = 0; x1 = x0 + 3;}
        else if(floor(u) == Width()-1) {x1 = floor(u); x0 = x1 - 3;}
        else if( floor(u) == Width()-2) {x0 = floor(u)-2; x1 = x0 + 3;}
        else {x0 = floor(u)-1; x1 = x0 +3;}
        // printf("test\n" );
        if(floor(v) == 0){ y0 = 0; y1 = y0 + 3;}
        else if(floor(v) == Height()-1) {y1 = floor(v); y0 = y1 - 3;}
        else if( floor(v) == Height()-2) {y0 = floor(v)-2; y1 = y0 + 3;}
        else {y0 = floor(v)-1; y1 = y0 +3;}
        // printf("test2\n" );
        // printf("%d %d %d %d\n",x0,x1,y0,y1 );
        r=0; g=0; b=0;
        for(x = x0; x <= x1; x++){
            for(y = y0; y <= y1; y++){
                r+= GetPixel(x,y).r*kernel[x-x0][y-y0];
                g+= GetPixel(x,y).g*kernel[x-x0][y-y0];
                b+= GetPixel(x,y).b*kernel[x-x0][y-y0];
            }
        }
        p.r = (int)(r+0.5);
        p.g = (int)(g+0.5);
        p.b = (int)(b+0.5);
        p.SetClamp(p.r,p.g,p.b);

    }
	return p;
}
