/* * * * * * * * * * * * * * * * * * * * * * * */
/* hw2.c * * * * * * * * * * * * * * * * * * * */
/* Created by: Jordan Bonecutter * * * * * * * */
/* 17 April 2019 * * * * * * * * * * * * * * * */
/* * * * * * * * * * * * * * * * * * * * * * * */

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <opencv2/opencv.hpp>

#include "complex.h"

#define ROWS	140
#define COLS	166
#define BSIZE	200
#define KSIZE	10
#define BRAD	2
#define SPD		(1.0)
#define SS		1
#define sqr(x)	((x)*(x))

typedef unsigned char gray_t;

int read_image(unsigned char image[ROWS][COLS], const char *ifile);
gray_t* read_ppm(const char *fname, int *rows, int *cols);
int save_ppm(const unsigned char *image, const char *fname, const int rows, const int cols);
comp *read_path(const char *fname, int *n, int *m);
void dft1(comp *im, comp *out, int c, double dir);
double mmod(double x, double y);
void subsample(comp* c, int n, int r);

int main(const int argc, const char **argv)
{
	/* Local Variables */
	int n, m, i, j;
	double t, s, r, ss;
	comp ctr, cprev;
	comp* path;
	comp* circ;
	gray_t out[ROWS][COLS];
	cv::VideoWriter writer;
	cv::Mat frame;
	comp* tmp;

	if(argc != 3)
	{
		printf("Incorrect usage: %s <input.ppm> <out.avi>\n", argv[0]);
		return 1;
	}

	/* Read in the path using the read_path function */
	path = read_path(argv[1], &n, &m);
	if(!(m%2))
	{
		m--;
	}

	/* Malloc our working arrays */
	circ = (comp*)malloc(sizeof(comp)*m/SPD + 1);
	tmp  = (comp*)malloc(sizeof(comp)*m);
	
	/* Structure for writing to a video file */
	writer = cv::VideoWriter(argv[2], cv::VideoWriter::fourcc('M', 'J', 'P', 'G'),
			60, cv::Size(2000, 2000), true);

	/* Take the forward dft of our data */
	dft1(path, tmp, m, -1.0);

	/* Take the inverse dft, but draw it out! */
	for(s = 0.; s < 2*m; s += SPD)
	{
		/* Our time variable should stay within m */
		t = fmod(s, m);

		/* Initially, set the center to the origin */
		ctr.re = ctr.im = 0.;

		/* Allocate a new frame */
		frame = NULL;
		frame = cv::Mat(2000, 2000, CV_8UC3);
		frame.setTo(cv::Scalar(0xff,0xff,0xff));

		/* Now, we need to add up the effect of each circle */
		for(j = 1; j < m-700; j++)
		{
			/* Order so that we pick low frequencies 
			 * before high frequencies */
			if(j%2 == 1) 
			{
				i = j/2 + 1;
			}
			else
			{
				i = m - j/2;
			}

			/* Draw a circle of size |F[in](u)| with its 
			 * center at the sum of the effect of all other
			 * circles */
			cv::circle(frame, cv::Point((ctr.re*9.)+1000, (ctr.im*9.)+1000), 9.*c_mag(tmp[i]),
					cv::Scalar(0xdf,0xdf,0xdf), 10);

			/* Move the center by F[in](u) */
			ctr = c_add(ctr, c_mulc(tmp[i],c_exp(mmod((2.0*M_PI*((double)i)/(double)m), M_PI)*s)));
		}

		/* Remember the circle so we can draw it out */
		circ[(int)(t/SPD)] = ctr;

		/* Draw all the circles from previous values */
		for(i = 4; i < t/SPD && i < n/SPD; i++)
		{
			cv::circle(frame, cv::Point((9.*circ[i].re)+1000, 
						(9.*circ[i].im)+1000), 8-(1./(0.25*(t/SPD-i)+0.75)), 
					cv::Scalar((int)(255./(0.25*(t/SPD-i)+0.75)),
						(int)(255./(0.25*(t/SPD-i)+0.75)),
						(int)(255./(0.25*(t/SPD-i)+0.75))),35);
		}
		if((int)(t/SPD) % 2)
			writer.write(frame);
	}

	/* Free mallocd memory */
	free(circ);
	free(tmp);
	free(path);

	return 0;
}

void subsample(comp* c, int n, int r)
{
	/* Simple nearest neighbors ss 
	 * by an integer value */
	int i;

	for(i = 0; i < n/r; i++)
	{
		c[i] = c[r*i];
	}
}

double mmod(double x, double y)
{
	/* Since frequencies greater than pi are 
	 * equivalent to that -f + 2pi we will create 
	 * a function to help us achieve this cycle:
	 *
	 * 0->pi -pi->0 */
	if((int)(x/y)%2)
	{
		return fmod(x,y) - y;
	}
	else
	{
		return fmod(x, y);
	}
	return x;
}

void dft1(comp *im, comp *out, int c, double dir)
{
	/* Simple naive implementation of dft */

	/* Local Var */
	int u, x;
	for(u = 0; u < c; u++)
	{
		/* Initialize accumlating variable to 0+j0 */
		out[u] = (comp){0,0};

		/* Implement 
		 * F[f](u) = SUM(x=0 -> c-1)[f(x)*exp(-j2pi*(u*x/M))]*/
		for(x = 0; x < c; x++)
		{
			out[u] = c_add(out[u],
					c_mulc(im[x], c_exp(dir*2.0*M_PI*((double)u*x)/(double)c)));
		}

		/* If we're taking the forward fourier, divide by c */
		if(dir < 0.)
			out[u] = c_mul(1.0/((c)), out[u]);
	}
}

gray_t* read_ppm(const char *fname, int* rows, int* cols)
{
	/* Open the file */
	FILE *fp = fopen(fname, "r");
	int ncount = 0, i, j;
	gray_t* image;

	/* Check if it exists */
	if(!fp)
	{
		return 0;
	}

	/* Get the image size and ignore magic number crap */
	/* TODO: figure out what the other stuff does */
	while(ncount < 1)
	{
		ncount += (fgetc(fp) == '\n');
	}
	fscanf(fp, "%d %d", cols, rows);
	image = (gray_t*)malloc(sizeof(gray_t)*(*rows)*(*cols));
	while(ncount < 3)
	{
		ncount += (fgetc(fp) == '\n');
	}

	/* Loop thru the image and read from the file */
	for(i = 0; i < *rows; i++)	
	{
		/* We just want a grayscale image so look 
		 * only at the green channel */
		for(j = 0; j < *cols; j++)
		{
			fgetc(fp);
			image[j+i*(*cols)] = (unsigned char)fgetc(fp);
			fgetc(fp);
		}
	}
	fclose(fp);

	return image;
}

comp* read_path(const char *fname, int *n, int *m)
{
	/* Local Var */
	gray_t *in;
	comp* ret;
	int x, y, reti, l, u, b, rows, cols;
	double dx, dy;

	/* Read the image */
	in  = read_ppm(fname, &rows, &cols);
	ret = (comp*)malloc(sizeof(comp)*rows*cols);

	/* Find a starting point */
	l = 1;
	for(x = 0; x < cols && l; x++)
	{
		for(y = 0; y <rows && l; y++)
		{
			if(in[(y)*cols+x] < 0x8f)
			{
				l = 0;	
			}
		}
	}

	/* Search for more points */
	l = -1;
	u = 1;
	for(reti = 0;; reti++)
	{
		if(in[cols*(y+u)+x] < 0x8f)
		{
			in[cols*(y+u)+x] = 0xff;

			y += u;
			ret[reti].re = x - COLS/2;
			ret[reti].im = y - ROWS/2;
		}
		else if(in[y*cols+x+l] < 0x8f)
		{
			in[y*cols+x+l] = 0xff;

			x += l;
			ret[reti].re = x - COLS/2;
			ret[reti].im = y - ROWS/2;
		}
		else if(in[cols*(y-u)+x] < 0x8f)
		{
			in[cols*(y-u)+x] = 0xff;

			y -= u;
			u *= -1;
			ret[reti].re = x - COLS/2;
			ret[reti].im = y - ROWS/2;
		}
		else if(in[y*cols+x-l] < 0x8f)
		{
			in[y*cols+x-l] = 0xff;

			x -= l;
			l *= -1;
			ret[reti].re = x - COLS/2;
			ret[reti].im = y - ROWS/2;
		}
		else if(in[cols*(y+2*u)+x] < 0x8f)
		{
			in[cols*(y+2*u)+x] = 0xff;

			y += 2*u;
			ret[reti].re = x - COLS/2;
			ret[reti].im = y - ROWS/2;
		}
		else if(in[y*cols+x+2*l] < 0x8f)
		{
			in[y*cols+x+2*l] = 0xff;

			x += 2*l;
			ret[reti].re = x - COLS/2;
			ret[reti].im = y - ROWS/2;
		}
		else if(in[cols*(y+u)+x+l] < 0x8f)
		{
			in[cols*(y+u)+x+l] = 0xff;

			x += l;
			y += u;
			ret[reti].re = x - COLS/2;
			ret[reti].im = y - ROWS/2;
		}
		else if(in[cols*(y+2*u)+x+l] < 0x8f)
		{
			in[cols*(y+2*u)+x+l] = 0xff;

			y += 2*u;
			x += l;
			ret[reti].re = x - COLS/2;
			ret[reti].im = y - ROWS/2;
		}
		else if(in[cols*(y+u)+x+2*l] < 0x8f)
		{
			in[cols*(y+u)+x+2*l] = 0xff;
			
			y += u;
			x += 2*l;
			ret[reti].re = x - COLS/2;
			ret[reti].im = y - ROWS/2;
		}
		else
		{
			b = 1;
			for(int xx = -3; xx <= 3 && b; xx++)
			{
				for(int yy = -3; yy <= 3 && b; yy++)
				{
					if(in[cols*(y+yy)+x+xx] < 0x8f)
					{
						u = yy < 0 ? -1 : 1;
						l = xx < 0 ? -1 : 1;

						y += yy;
						x += xx;
						ret[reti].re = x - COLS/2;
						ret[reti].im = y - ROWS/2;
						b = 0;
					}
				}
			}
			if(b)
			{
				break;
			}
		}
	}
	free(in);
	*n = reti;
	if(c_mag(c_add(c_mul(-1, ret[reti-1]), ret[0])) < 3.)
	{
		if(*n % 2 == 0)	
		{
			(*n)--;
			*m = *n;
		}
		return ret;
	}
	dx = (-ret[reti-1].re + ret[0].re)/120.;
	dy = (-ret[reti-1].im + ret[0].im)/120.;
	
	for(; reti < *n + 120; reti++)
	{
		ret[reti].re = ret[reti-1].re + dx;
		ret[reti].im = ret[reti-1].im + dy;
	}
	if(*n % 2 == 0)
	{	
		*m = *n + 121;
	}
	else
	{
		*m = *n + 120;
	}
	for(reti = 3; reti < *m; reti++)
	{
		ret[reti] = c_mul(0.25, c_add(ret[reti], 
					c_add(ret[reti-1], c_add(ret[reti-2], ret[reti-3]))));
	}
	return ret;
}

int save_ppm(const unsigned char *image, const char *fname, const int rows, const int cols)
{
    // Local var
    char *s, *t, c, buff[BSIZE];
    FILE *fp;
    unsigned int i, j;

    // Assert
    assert(image);
    assert(fname);

    // Insert ppm filename if it wasn't passed w/ fname
    if(!strstr(fname, ".ppm"))
    {
        s = (char *)malloc(sizeof(char) * (strlen(fname) + 5));
        memcpy(s, fname, strlen(fname));
        t = s;
        s += strlen(fname);
        *s++ = '.';
        *s++ = 'p';
        *s++ = 'p';
        *s++ = 'm';
        *s = 0;
        s = t;
    }
    else
    {
        s = (char *)fname;
    }

    // Memset buffer
    memset(buff, 0, BSIZE);

    // Open the file!
    fp = fopen(s, "w");

    // If s was alloc'd
    if(s != fname)
    {
        free(s);
    }

    // If the file couldn't be opened
    if(!fp)
    {
#ifndef NDEBUG
        printf("Invalid file name!\n");
#endif
        return -1;
    }

    // Magic string for ppm
    fputs("P6\n", fp);

    // Put image size
    sprintf(buff, "%d %d\n", cols, rows);
    fputs(buff, fp);
    memset(buff, 0, BSIZE);
    // Max pixel is 255
    fputs("255\n", fp);


    // Loop!
    for(i = 0; i < rows; i++)
    {
        for(j = 0; j < cols; j++)
        {
            c = (char)image[j + (i*cols)];
            fputc(c, fp);
            fputc(c, fp);
            fputc(c, fp);
        }
    }

    fclose(fp);

    return 1;
}

int read_image(unsigned char image[ROWS][COLS], const char *ifile)
{
    // Local var
    FILE *fp;
    int i; 

    // Assert
    assert(image);
    assert(ifile);

    // Try to open the file
    if((fp = fopen(ifile, "rb")) == NULL)
    {
        return 0;
    }

    // Read the file!
    for(i = 0; i < ROWS; i++)
    {
        if(fread(image[i], 1, COLS, fp) != COLS)
        {
            return 0;
        }
    }

    // Return success
    fclose(fp);
    return 1;
}
