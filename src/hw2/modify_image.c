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
    /*
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
    */
   float v1, v2, v3, v4, q1, q2, q;
    v1 = get_pixel(im, floorf(x), floorf(y), c);
    v2 = get_pixel(im, ceilf(x), floorf(y), c);
    v3 = get_pixel(im, floorf(x), ceilf(y), c);
    v4 = get_pixel(im, ceilf(x), ceilf(y), c);
    q1 = v2 * (x - floorf(x)) + v1 * (ceilf(x) - x);
    q1 = floorf(x) == ceilf(x) ? v1 : q1;
    q2 = v4 * (x - floorf(x)) + v3 * (ceilf(x) - x);
    q2 = floorf(x) == ceilf(x) ? v3 : q2;
    q = q2 * (y - floorf(y)) + q1 * (ceilf(y) - y);
    q = floorf(y) == ceilf(y) ? q1 : q;
    return q;
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

    if (s < 1e-6) return;

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
    image im = make_image(w,w,1);
    for (int i = 0; i < im.w; i ++) {
      for (int j = 0; j < im.h; j ++) {
        set_pixel(im, i, j, 0, 1.0);
      }
    }
    l1_normalize(im);
    return im;
}

image convolve_image(image im, image filter, int preserve)
{
    // TODO
    /***********************************************************************
      This function convolves the image "im" with the "filter". The value
      of preserve is 1 if the number of input image channels need to be 
      preserved. Check the detailed algorithm given in the README.  
    ************************************************************************/
    assert(im.c == filter.c || filter.c == 1);
    image res = make_image(im.w, im.h, preserve ? im.c : 1);
    int x_center = filter.w / 2;
    int y_center = filter.h / 2;
    for (int i = 0; i < im.c; i ++) {
      for (int j = 0; j < im.h; j ++) {
        for (int k = 0; k < im.w; k ++) {
          float q = 0.0;
          for (int fy = 0; fy < filter.h; fy ++) {
            for (int fx = 0; fx < filter.w; fx ++) {
              int fc = im.c == filter.c ? i : 0;
              float ival = get_pixel(im, k-x_center+fx, j-y_center+fy, i);
              float fval = get_pixel(filter, fx, fy, fc);
              q += ival * fval;
            }
          }
          if (preserve) set_pixel(res, k, j, i, q);
          else {
            float original = !i ? 0.0 : get_pixel(res, k, j, 0);
            set_pixel(res, k, j, 0, original+q);
          }
        }
      }
    }
    return res;
}

image make_highpass_filter()
{
    // TODO
    /***********************************************************************
      Create a 3x3 filter with highpass filter values using image.data[]
    ************************************************************************/
    image filter = make_image(3,3,1);
    float a[] = {0, -1, 0, -1, 4, -1, 0, -1, 0};
    int cnt = 0;
    for (int i = 0; i < filter.w; i ++) {
      for (int j = 0; j < filter.h; j ++) {
        set_pixel(filter, i, j, 0, a[cnt++]);
      }
    }
    return filter;
}

image make_sharpen_filter()
{
    // TODO
    /***********************************************************************
      Create a 3x3 filter with sharpen filter values using image.data[]
    ************************************************************************/
    image filter = make_image(3,3,1);
    float a[] = {0, -1, 0, -1, 5, -1, 0, -1, 0};
    int cnt = 0;
    for (int i = 0; i < filter.w; i ++) {
      for (int j = 0; j < filter.h; j ++) {
        set_pixel(filter, i, j, 0, a[cnt++]);
      }
    }
    return filter;
}

image make_emboss_filter()
{
    // TODO
    /***********************************************************************
      Create a 3x3 filter with emboss filter values using image.data[]
    ************************************************************************/
    image filter = make_image(3,3,1);
    float a[] = {-2, -1, 0, -1, 1, 1, 0, 1, 2};
    int cnt = 0;
    for (int i = 0; i < filter.w; i ++) {
      for (int j = 0; j < filter.h; j ++) {
        set_pixel(filter, i, j, 0, a[cnt++]);
      }
    }
    return filter;
}

// Question 2.3.1: Which of these filters should we use preserve when we run our convolution and which ones should we not? Why?
// Answer: TODO
/*
We should use preserve for the sharpen filter and the emboss filter. This is because the sharpen filter sharpens edges and lines
in an image, but we still want to see colors. The emboss filter is used to emphasize the difference of intensity neibouring pixels.
Therefore, 3 channels are needed for these two filters.
We do not use preserve for the highpass filter. This is because the highpass filter is used to detect and find edges in an image.
So, colors are not required. Also, as we notes that the for neibouring pixels with the same colors, they would result zero(black)
through the highpass filter, so one channel is enough.
*/

// Question 2.3.2: Do we have to do any post-processing for the above filters? Which ones and why?
// Answer: TODO
/*
The sharpen filter and the emboss filter require post-processing. This is because after we process covolution with these two filters,
some pixels might be overflowed in value, which results in imcompatible colors with neibouring pixels.
*/

image make_gaussian_filter(float sigma)
{
    // TODO
    /***********************************************************************
      sigma: a float number for the Gaussian.
      Create a Gaussian filter with the given sigma. Note that the kernel size 
      is the next highest odd integer from 6 x sigma. Return the Gaussian filter.
    ************************************************************************/
    int w = 6 * sigma;
    w = w % 2 ? w : w + 1;
    image filter = make_image(w, w, 1);
    for (int i = 0; i < w; i ++) {
      for (int j = 0; j < w; j ++) {
        float x = i - w/2;
        float y = j - w/2;
        float a = 1.0 / (TWOPI * sigma * sigma);
        float b = exp(-(x*x+y*y)/(2*sigma*sigma));
        set_pixel(filter, i, j, 0, a*b);
      }
    }
    l1_normalize(filter);
    return filter;
}

image add_image(image a, image b)
{
    // TODO
    /***********************************************************************
      The input images a and image b have the same height, width, and channels.
      Sum the given two images and return the result, which should also have
      the same height, width, and channels as the inputs. Do necessary checks.
    ************************************************************************/
    image im = make_image(a.w, a.h, a.c);
    for (int i = 0; i < a.c; i ++) {
      for (int j = 0; j < a.h; j ++) {
        for (int k = 0; k < a.w; k ++) {
          float a_val = get_pixel(a, k, j, i);
          float b_val = get_pixel(b, k, j, i);
          set_pixel(im, k, j, i, a_val + b_val);
        }
      }
    }
    return im;
}

image sub_image(image a, image b)
{
    // TODO
    /***********************************************************************
      The input image a and image b have the same height, width, and channels.
      Subtract the given two images and return the result, which should have
      the same height, width, and channels as the inputs. Do necessary checks.
    ************************************************************************/
    image im = make_image(a.w, a.h, a.c);
    for (int i = 0; i < a.c; i ++) {
      for (int j = 0; j < a.h; j ++) {
        for (int k = 0; k < a.w; k ++) {
          float a_val = get_pixel(a, k, j, i);
          float b_val = get_pixel(b, k, j, i);
          set_pixel(im, k, j, i, a_val - b_val);
        }
      }
    }
    return im;
}

image make_gx_filter()
{
    // TODO
    /***********************************************************************
      Create a 3x3 Sobel Gx filter and return it
    ************************************************************************/
    image filter = make_image(3,3,1);
    float a[] = {-1, 0, 1, -2, 0, 2, -1, 0, 1};
    int cnt = 0;
    for (int i = 0; i < filter.h; i ++) {
      for (int j = 0; j < filter.w; j ++) {
        set_pixel(filter, j, i, 0, a[cnt++]);
      }
    }
    return filter;
}

image make_gy_filter()
{
    // TODO
    /***********************************************************************
      Create a 3x3 Sobel Gy filter and return it
    ************************************************************************/
    image filter = make_image(3,3,1);
    float a[] = {-1, -2, -1, 0, 0, 0, 1, 2, 1};
    int cnt = 0;
    for (int i = 0; i < filter.h; i ++) {
      for (int j = 0; j < filter.w; j ++) {
        set_pixel(filter, j, i, 0, a[cnt++]);
      }
    }
    return filter;
}

void feature_normalize(image im)
{
    // TODO
    /***********************************************************************
      Calculate minimum and maximum pixel values. Normalize the image by
      subtracting the minimum and dividing by the max-min difference.
    ************************************************************************/
    float min_val = get_pixel(im, 0, 0, 0);
    float max_val = get_pixel(im, 0, 0, 0);
    for (int i = 0; i < im.c; i ++) {
      for (int j = 0; j < im.h; j ++) {
        for (int k = 0; k < im.w; k ++) {
          float t = get_pixel(im, k, j, i);
          if (min_val > t) min_val = t;
          if (max_val < t) max_val = t;
        }
      }
    }

    for (int i = 0; i < im.c; i ++) {
      for (int j = 0; j < im.h; j ++) {
        for (int k = 0; k < im.w; k ++) {
          float v = (get_pixel(im, k, j, i) - min_val) / (max_val - min_val);
          set_pixel(im, k, j, i, v);
        }
      }
    }

}

image *sobel_image(image im)
{
    // TODO
    /***********************************************************************
      Apply Sobel filter to the input image "im", get the magnitude as sobelimg[0]
      and gradient as sobelimg[1], and return the result.
    ************************************************************************/
    image *sobelimg = calloc(2, sizeof(image));
    image gx = make_gx_filter();
    image gy = make_gy_filter();
    image gx_im = convolve_image(im, gx, 0);
    image gy_im = convolve_image(im, gy, 0);

    image mag = make_image(im.w, im.h, 1);
    image dir = make_image(im.w, im.h, 1);
    for (int i = 0; i < im.h; i ++) {
      for (int j = 0; j < im.w; j ++) {
        float x = get_pixel(gx_im, j, i, 0);
        float y = get_pixel(gy_im, j, i, 0);
        set_pixel(mag, j, i, 0, sqrt(x*x + y*y));
        set_pixel(dir, j, i, 0, atan2(y, x));
      }
    }
    sobelimg[0] = mag;
    sobelimg[1] = dir;

    free_image(gx);
    free_image(gy);
    free_image(gx_im);
    free_image(gy_im);
    return sobelimg;
}

image colorize_sobel(image im)
{
  // TODO
  /***********************************************************************
    Create a colorized version of the edges in image "im" using the 
    algorithm described in the README.
  ************************************************************************/
  image* sobelimg = sobel_image(im);
  image mag = sobelimg[0];
  image dir = sobelimg[1];

  feature_normalize(mag);
  feature_normalize(dir);

  image res = make_image(im.w, im.h, 3);
  for (int i = 0; i < im.h; i ++) {
    for (int j = 0; j < im.w; j ++) {
      float m = get_pixel(mag, j, i, 0);
      float d = get_pixel(dir, j, i, 0);
      set_pixel(res, j, i, 0, d);
      set_pixel(res, j, i, 1, m);
      set_pixel(res, j, i, 2, m);
    }
  }

  hsv_to_rgb(res);
  return res;
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