#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "image.h"
#define TWOPI 6.2831853

/******************************** Resizing *****************************
  To resize we'll need some interpolation methods and a function to create
  a new image and fill it in with our interpolation methods.
************************************************************************/

float nn_interpolate(image im, float x, float y, int c)
{
    // TODO
    /***********************************************************************
      This function performs nearest-neighbor interpolation on image "im"
      given a floating column value "x", row value "y" and integer channel "c",
      and returns the interpolated value.
    ************************************************************************/
    int x_ = round(x);
    int y_ = round(y);
    if (x_ < 0) x_ = 0;
    if (x_ >= im.w) x_ = im.w - 1;
    if (y_ < 0) y_ = 0;
    if (y_ >= im.h) y_ = im.h - 1;
    return get_pixel(im, x_, y_, c);
    
}

image nn_resize(image im, int w, int h)
{
    // TODO Fill in (also fix the return line)
    /***********************************************************************
      This function uses nearest-neighbor interpolation on image "im" to a new
      image of size "w x h"
    ************************************************************************/
    float a_x = im.w * 1.0 / w;
    float a_y = im.h * 1.0 / h;
    float b_x = (im.w - w) * 1.0 / 2 / w;
    float b_y = (im.h - h) * 1.0 / 2 / h;
    image res_im = make_image(w, h, im.c);
    for (int i = 0; i < im.c; i ++) {
      for (int j = 0; j < h; j ++) {
        for (int k = 0; k < w; k ++) {
          int pos = w * h * i + w * j + k;
          float x_ = k * a_x + b_x, y_ = j * a_y + b_y;
          res_im.data[pos] = nn_interpolate(im, x_, y_, i);
        }
      }
    }
    return res_im;
}

float bilinear_interpolate(image im, float x, float y, int c)
{
    // TODO
    /***********************************************************************
      This function performs bilinear interpolation on image "im" given
      a floating column value "x", row value "y" and integer channel "c".
      It interpolates and returns the interpolated value.
    ************************************************************************/
    int low_x = floor(x);
    int low_y = floor(y);
    int up_x = ceil(x);
    int up_y = ceil(y);
    float a1 = (up_x - x) * (up_y - y);
    float a2 = (x - low_x) * (up_y - y);
    float a3 = (up_x - x) * (y - low_y);
    float a4 = (x - low_x) * (y - low_y);
    float v1 = get_pixel(im, low_x, low_y, c);
    float v2 = get_pixel(im, up_x, low_y, c);
    float v3 = get_pixel(im, low_x, up_y, c);
    float v4 = get_pixel(im, up_x, up_y, c);

    return v1 * a1 + v2 * a2 + v3 * a3 + v4 * a4;
}

image bilinear_resize(image im, int w, int h)
{
    // TODO
    /***********************************************************************
      This function uses bilinear interpolation on image "im" to a new image
      of size "w x h". Algorithm is same as nearest-neighbor interpolation.
    ************************************************************************/
    float a_x = im.w * 1.0 / w;
    float a_y = im.h * 1.0 / h;
    float b_x = (im.w - w) * 1.0 / 2 / w;
    float b_y = (im.h - h) * 1.0 / 2 / h;
    image res_im = make_image(w, h, im.c);
    for (int i = 0; i < im.c; i ++) {
      for (int j = 0; j < h; j ++) {
        for (int k = 0; k < w; k ++) {
          int pos = w * h * i + w * j + k;
          float x_ = k * a_x + b_x, y_ = j * a_y + b_y;
          res_im.data[pos] = bilinear_interpolate(im, x_, y_, i);
        }
      }
    }
    return res_im;
}


/********************** Filtering: Box filter ***************************
  We want to create a box filter. We will only use square box filters.
************************************************************************/

void l1_normalize(image im)
{
    // TODO
    /***********************************************************************
      This function divides each value in image "im" by the sum of all the
      values in the image and modifies the image in place.
    ************************************************************************/
    float s = 0.0;
    for (int i = 0; i < im.c; i ++) {
      for (int j = 0; j < im.h; j ++) {
        for (int k = 0; k < im.w; k ++) {
          s += get_pixel(im, k, j, i);
        }
      }
    }

    for (int i = 0; i < im.c; i ++) {
      for (int j = 0; j < im.h; j ++) {
        for (int k = 0; k < im.w; k ++) {
          float t = get_pixel(im, k, j, i) / s;
          set_pixel(im, k, j, i, t);
        }
      }
    }
}

image make_box_filter(int w)
{
    // TODO
    /***********************************************************************
      This function makes a square filter of size "w x w". Make an image of
      width = height = w and number of channels = 1, with all entries equal
      to 1. Then use "l1_normalize" to normalize your filter.
    ************************************************************************/
    return make_image(1,1,1);
}

image convolve_image(image im, image filter, int preserve)
{
    // TODO
    /***********************************************************************
      This function convolves the image "im" with the "filter". The value
      of preserve is 1 if the number of input image channels need to be 
      preserved. Check the detailed algorithm given in the README.  
    ************************************************************************/
    return make_image(1,1,1);
}

image make_highpass_filter()
{
    // TODO
    /***********************************************************************
      Create a 3x3 filter with highpass filter values using image.data[]
    ************************************************************************/
    return make_image(1,1,1);
}

image make_sharpen_filter()
{
    // TODO
    /***********************************************************************
      Create a 3x3 filter with sharpen filter values using image.data[]
    ************************************************************************/
    return make_image(1,1,1);
}

image make_emboss_filter()
{
    // TODO
    /***********************************************************************
      Create a 3x3 filter with emboss filter values using image.data[]
    ************************************************************************/
    return make_image(1,1,1);
}

// Question 2.3.1: Which of these filters should we use preserve when we run our convolution and which ones should we not? Why?
// Answer: TODO

// Question 2.3.2: Do we have to do any post-processing for the above filters? Which ones and why?
// Answer: TODO

image make_gaussian_filter(float sigma)
{
    // TODO
    /***********************************************************************
      sigma: a float number for the Gaussian.
      Create a Gaussian filter with the given sigma. Note that the kernel size 
      is the next highest odd integer from 6 x sigma. Return the Gaussian filter.
    ************************************************************************/
    return make_image(1,1,1);
}

image add_image(image a, image b)
{
    // TODO
    /***********************************************************************
      The input images a and image b have the same height, width, and channels.
      Sum the given two images and return the result, which should also have
      the same height, width, and channels as the inputs. Do necessary checks.
    ************************************************************************/
    return make_image(1,1,1);
}

image sub_image(image a, image b)
{
    // TODO
    /***********************************************************************
      The input image a and image b have the same height, width, and channels.
      Subtract the given two images and return the result, which should have
      the same height, width, and channels as the inputs. Do necessary checks.
    ************************************************************************/
    return make_image(1,1,1);
}

image make_gx_filter()
{
    // TODO
    /***********************************************************************
      Create a 3x3 Sobel Gx filter and return it
    ************************************************************************/
    return make_image(1,1,1);
}

image make_gy_filter()
{
    // TODO
    /***********************************************************************
      Create a 3x3 Sobel Gy filter and return it
    ************************************************************************/
    return make_image(1,1,1);
}

void feature_normalize(image im)
{
    // TODO
    /***********************************************************************
      Calculate minimum and maximum pixel values. Normalize the image by
      subtracting the minimum and dividing by the max-min difference.
    ************************************************************************/
}

image *sobel_image(image im)
{
    // TODO
    /***********************************************************************
      Apply Sobel filter to the input image "im", get the magnitude as sobelimg[0]
      and gradient as sobelimg[1], and return the result.
    ************************************************************************/
    image *sobelimg = calloc(2, sizeof(image));
    return sobelimg;
}

image colorize_sobel(image im)
{
  // TODO
  /***********************************************************************
    Create a colorized version of the edges in image "im" using the 
    algorithm described in the README.
  ************************************************************************/
  return make_image(1,1,1);
}

// EXTRA CREDIT: Median filter

/*
image apply_median_filter(image im, int kernel_size)
{
  return make_image(1,1,1);
}
*/

// SUPER EXTRA CREDIT: Bilateral filter

/*
image apply_bilateral_filter(image im, float sigma1, float sigma2)
{
  return make_image(1,1,1);
}
*/