#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "image.h"
#include "matrix.h"
#include <time.h>

// Frees an array of descriptors.
// descriptor *d: the array.
// int n: number of elements in array.
void free_descriptors(descriptor *d, int n)
{
    int i;
    for(i = 0; i < n; ++i){
        free(d[i].data);
    }
    free(d);
}

// Create a feature descriptor for an index in an image.
// image im: source image.
// int i: index in image for the pixel we want to describe.
// returns: descriptor for that index.
descriptor describe_index(image im, int i)
{
    int w = 5;
    descriptor d;
    d.p.x = i%im.w;
    d.p.y = i/im.w;
    d.data = calloc(w*w*im.c, sizeof(float));
    d.n = w*w*im.c;
    int c, dx, dy;
    int count = 0;
    // If you want you can experiment with other descriptors
    // This subtracts the central value from neighbors
    // to compensate some for exposure/lighting changes.
    for(c = 0; c < im.c; ++c){
        float cval = im.data[c*im.w*im.h + i];
        for(dx = -w/2; dx < (w+1)/2; ++dx){
            for(dy = -w/2; dy < (w+1)/2; ++dy){
                float val = get_pixel(im, i%im.w+dx, i/im.w+dy, c);
                d.data[count++] = cval - val;
            }
        }
    }
    return d;
}

// Marks the spot of a point in an image.
// image im: image to mark.
// ponit p: spot to mark in the image.
void mark_spot(image im, point p)
{
    int x = p.x;
    int y = p.y;
    int i;
    for(i = -9; i < 10; ++i){
        set_pixel(im, x+i, y, 0, 1);
        set_pixel(im, x, y+i, 0, 1);
        set_pixel(im, x+i, y, 1, 0);
        set_pixel(im, x, y+i, 1, 0);
        set_pixel(im, x+i, y, 2, 1);
        set_pixel(im, x, y+i, 2, 1);
    }
}

// Marks corners denoted by an array of descriptors.
// image im: image to mark.
// descriptor *d: corners in the image.
// int n: number of descriptors to mark.
void mark_corners(image im, descriptor *d, int n)
{
    int i;
    for(i = 0; i < n; ++i){
        mark_spot(im, d[i].p);
    }
}

// Creates a 1d Gaussian filter.
// float sigma: standard deviation of Gaussian.
// returns: single row image of the filter.
image make_1d_gaussian(float sigma)
{
    // TODO: make separable 1d Gaussian.
    int w = 6 * sigma;
    w = w % 2 ? w : w + 1;
    image filter = make_image(w, 1, 1);
    for (int i = 0; i < w; i ++) {
      float x = i - w/2;
      float a = 1.0 / (sqrt(TWOPI) * sigma);
      float b = exp(-(x*x)/(2*sigma*sigma));
      set_pixel(filter, i, 0, 0, a*b);
    }
    l1_normalize(filter);
    return filter;
}

// Smooths an image using separable Gaussian filter.
// image im: image to smooth.
// float sigma: std dev. for Gaussian.
// returns: smoothed image.
image smooth_image(image im, float sigma)
{
    // TODO: use two convolutions with 1d gaussian filter.
    image S = make_image(im.w, im.h, 3);
    image gau_1d = make_1d_gaussian(sigma);
    image transpose_gau = make_image(gau_1d.h, gau_1d.w, 1);
    for (int i = 0; i < gau_1d.w; i ++) {
        set_pixel(transpose_gau, 0, i, 0, get_pixel(gau_1d, i, 0, 0));
    }
    S = convolve_image(im, transpose_gau, 1);
    S = convolve_image(S, gau_1d, 1);

    free_image(gau_1d);
    free_image(transpose_gau);

    return S;
}

// Calculate the structure matrix of an image.
// image im: the input image.
// float sigma: std dev. to use for weighted sum.
// returns: structure matrix. 1st channel is Ix^2, 2nd channel is Iy^2,
//          third channel is IxIy.
image structure_matrix(image im, float sigma)
{
    image S = make_image(im.w, im.h, 3);
    // TODO: calculate structure matrix for im.
    image gx = make_gx_filter();
    image gy = make_gy_filter();
    image Ix = convolve_image(im, gx, 0);
    image Iy = convolve_image(im, gy, 0);
    for (int i = 0; i < im.h; i ++) {
      for (int j = 0; j < im.w; j ++) {
        float ix = get_pixel(Ix, j, i, 0);
        float iy = get_pixel(Iy, j, i, 0);
        set_pixel(S, j, i, 0, ix * ix);
        set_pixel(S, j, i, 1, iy * iy);
        set_pixel(S, j, i, 2, ix * iy);
      }
    }
    /* Ordinary way
    image gau = make_gaussian_filter(sigma);
    S = convolve_image(S, gau, 1);
    free_image(gau);
    */

    S = smooth_image(S, sigma);

    free_image(gx);
    free_image(gy);
    free_image(Ix);
    free_image(Iy);

    return S;
}

// Estimate the cornerness of each pixel given a structure matrix S.
// image S: structure matrix for an image.
// returns: a response map of cornerness calculations.
image cornerness_response(image S)
{
    image R = make_image(S.w, S.h, 1);
    // TODO: fill in R, "cornerness" for each pixel using the structure matrix.
    // We'll use formulation det(S) - alpha * trace(S)^2, alpha = .06.
    float alpha = .06;
    for (int i = 0; i < S.h; i ++) {
      for (int j = 0; j < S.w; j ++) {
        float ixx = get_pixel(S, j, i, 0);
        float iyy = get_pixel(S, j, i, 1);
        float ixy = get_pixel(S, j, i, 2);
        float det = ixx * iyy - ixy * ixy;
        float trace = ixx + iyy;
        set_pixel(R, j, i, 0, det - alpha * trace * trace);
      }
    }
    return R;
}

// Perform non-max supression on an image of feature responses.
// image im: 1-channel image of feature responses.
// int w: distance to look for larger responses.
// returns: image with only local-maxima responses within w pixels.
image nms_image(image im, int w)
{
    image r = copy_image(im);
    // TODO: perform NMS on the response map.
    // for every pixel in the image:
    //     for neighbors within w:
    //         if neighbor response greater than pixel response:
    //             set response to be very low (I use -999999 [why not 0??])
    for (int i = 0; i < im.h; i ++) {
      for (int j = 0; j < im.w; j ++) {
        float val = get_pixel(im, j, i, 0);
        for (int dx = -w; dx <= w; dx ++) {
          for (int dy = -w; dy <= w; dy ++) {
            float chk = get_pixel(im, j+dx, i+dy, 0);
            if (chk > get_pixel(im, j, i, 0)) set_pixel(r, j, i, 0, -1e7);
          }
        }
      }
    }
    return r;
}

// Perform harris corner detection and extract features from the corners.
// image im: input image.
// float sigma: std. dev for harris.
// float thresh: threshold for cornerness.
// int nms: distance to look for local-maxes in response map.
// int *n: pointer to number of corners detected, should fill in.
// returns: array of descriptors of the corners in the image.
descriptor *harris_corner_detector(image im, float sigma, float thresh, int nms, int *n)
{
    // Calculate structure matrix
    image S = structure_matrix(im, sigma);

    // Estimate cornerness
    image R = cornerness_response(S);

    // Run NMS on the responses
    image Rnms = nms_image(R, nms);


    //TODO: count number of responses over threshold
    int count = 0;
    for (int i = 0; i < Rnms.c; i ++) {
      for (int j = 0; j < Rnms.h; j ++) {
        for (int k = 0; k < Rnms.w; k ++) {
            if (get_pixel(Rnms, k, j, i) > thresh) count ++;
        }
      }
    }

    
    *n = count; // <- set *n equal to number of corners in image.
    descriptor *d = calloc(count, sizeof(descriptor));
    //TODO: fill in array *d with descriptors of corners, use describe_index.
    int len = Rnms.c * Rnms.h * Rnms.w;
    for (int i = 0, j = 0; i < len; i ++) {
      if (Rnms.data[i] > thresh) {
        d[j++] = describe_index(im, i);
      }
    }

    free_image(S);
    free_image(R);
    free_image(Rnms);
    return d;
}

// Find and draw corners on an image.
// image im: input image.
// float sigma: std. dev for harris.
// float thresh: threshold for cornerness.
// int nms: distance to look for local-maxes in response map.
void detect_and_draw_corners(image im, float sigma, float thresh, int nms)
{
    int n = 0;
    descriptor *d = harris_corner_detector(im, sigma, thresh, nms, &n);
    mark_corners(im, d, n);
}
