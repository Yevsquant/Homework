#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include "image.h"

float get_pixel(image im, int x, int y, int c)
{
    // TODO Fill this in
    // x: width, y: height
    if (x < 0) x = 0;
    if (x >= im.w) x = im.w - 1;
    if (y < 0) y = 0;
    if (y >= im.h) y = im.h - 1;
    return *(im.data + im.w * im.h * c + im.w * y + x);
}

void set_pixel(image im, int x, int y, int c, float v)
{
    // TODO Fill this in
    *(im.data + im.w * im.h * c + im.w * y + x) = v;
}

image copy_image(image im)
{
    image copy = make_image(im.w, im.h, im.c);
    // TODO Fill this in
    memcpy(copy.data, im.data, im.w * im.h * im.c * sizeof(float));
    return copy;
}

image rgb_to_grayscale(image im)
{
    assert(im.c == 3);
    image gray = make_image(im.w, im.h, 1);
    // TODO Fill this in
    for (int i = 0; i < im.h; i ++) {
        for (int j = 0; j < im.w; j ++) {
            float value = get_pixel(im, j, i, 0) * 0.299
                        + get_pixel(im, j, i, 1) * 0.587
                        + get_pixel(im, j, i, 2) * 0.114;
            set_pixel(gray, j, i, 0, value);
            
        }
    }
    return gray;
}

void shift_image(image im, int c, float v)
{
    // TODO Fill this in
    if (c < 0 || c >= im.c) return;
    for (int i = 0; i < im.h; i ++) {
        for (int j = 0; j < im.w; j ++) {
            float value = get_pixel(im, j, i, c) + v;
            set_pixel(im, j, i, c, value);
        }
    }
}

void clamp_image(image im)
{
    // TODO Fill this in
    for (int i = 0; i < im.c; i ++) {
        for (int j = 0; j < im.h; j ++) {
            for (int k = 0; k < im.w; k ++) {
                float val = get_pixel(im, k, j, i);
                if (val < 0) val = 0.0;
                if (val > 1) val = 1.0;
                set_pixel(im, k, j, i, val);
            }
        }
    }
}


// These might be handy
float three_way_max(float a, float b, float c)
{
    return (a > b) ? ( (a > c) ? a : c) : ( (b > c) ? b : c) ;
}

float three_way_min(float a, float b, float c)
{
    return (a < b) ? ( (a < c) ? a : c) : ( (b < c) ? b : c) ;
}

void rgb_to_hsv(image im)
{
    // TODO Fill this in
    assert(im.c == 3);
    for (int i = 0; i < im.h; i ++) {
        for (int j = 0; j < im.w; j ++) {
            float r = get_pixel(im, j, i, 0);
            float g = get_pixel(im, j, i, 1);
            float b = get_pixel(im, j, i, 2);
            float v = three_way_max(r, g, b);
            float c = v - three_way_min(r, g, b);
            float s = fabs(v) < 1e-6 ? 0.0 : c / v;
            float h = 0.0;
            if (fabs(c) >= 1e-6) {
                if (v == r) h = (g - b) / c;
                else if (v == g) h = (b - r) / c + 2;
                else h = (r - g) / c + 4;
            }
            h /= 6;
            if (h < 0) h ++;
            set_pixel(im, j, i, 0, h);
            set_pixel(im, j, i, 1, s);
            set_pixel(im, j, i, 2, v);
        }
    }
}

void hsv_to_rgb(image im)
{
    // TODO Fill this in
    assert(im.c == 3);
    void change (image im, float r, float g, float b, int i, int j);
    for (int i = 0; i < im.h; i ++) {
        for (int j = 0; j < im.w; j ++) {
            float h = get_pixel(im, j, i, 0) * 6;
            float s = get_pixel(im, j, i, 1);
            float v = get_pixel(im, j, i, 2);
            float hi = floor(h);
            float f = h - hi;
            float p = v * (1 - s);
            float q = v * (1 - f * s);
            float t = v * (1 - (1 - f) * s);
            if (fabs(hi) < 1e-6) change(im, v, t, p, i, j);
            if (fabs(hi - 1) < 1e-6) change(im, q, v, p, i, j);
            if (fabs(hi - 2) < 1e-6) change(im, p, v, t, i, j);
            if (fabs(hi - 3) < 1e-6) change(im, p, q, v, i, j);
            if (fabs(hi - 4) < 1e-6) change(im, t, p, v, i, j);
            if (fabs(hi - 5) < 1e-6) change(im, v, p, q, i, j);
        }
    }
}

// helper functions
void change (image im, float r, float g, float b, int i, int j) {
    set_pixel(im, j, i, 0, r);
    set_pixel(im, j, i, 1, g);
    set_pixel(im, j, i, 2, b);
}

void scale_image(image im, int c, float v) {
    if (c < 0 || c >= im.c) return;
    for (int i = 0; i < im.h; i ++) {
        for (int j = 0; j < im.w; j ++) {
            float value = get_pixel(im, j, i, c) * v;
            set_pixel(im, j, i, c, value);
        }
    }
}