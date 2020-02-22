#include "rasterizer.h"

using namespace std;

namespace CGL {

RasterizerImp::RasterizerImp(PixelSampleMethod psm, LevelSampleMethod lsm,
                                       size_t width, size_t height,
                                       unsigned int sample_rate) {
  this->psm = psm;
  this->lsm = lsm;
  this->width = width;
  this->height = height;
  this->sample_rate = sample_rate;

  supersample_buffer.resize(width * height * sample_rate, Color::White);
}

// Used by rasterize_point and rasterize_line
void RasterizerImp::fill_pixel(size_t x, size_t y, Color c) {
  // TODO: Task 2: You might need to this function to fix points and lines (such as the black rectangle border in test4.svg)
  // NOTE: You are not required to implement proper supersampling for points and lines
  // It is sufficient to use the same color for all supersamples of a pixel for points and lines (not triangles)

        rgb_framebuffer_target[3 * (y * width + x)    ] = (unsigned char)(c.r * 255);
        rgb_framebuffer_target[3 * (y * width + x) + 1] = (unsigned char)(c.g * 255);
        rgb_framebuffer_target[3 * (y * width + x) + 2] = (unsigned char)(c.b * 255);

}

// Optional helper function to add a sample to the supersample_buffer
void RasterizerImp::fill_supersample(size_t x, size_t y, size_t s, Color c) {
  // TODO: Task 2: You may want to implement this function. Hint: our solution uses one line
    supersample_buffer[width * y * sample_rate + x * sample_rate + s] = c;
};

// Rasterize a point: simple example to help you start familiarizing
// yourself with the starter code.
//
void RasterizerImp::rasterize_point(float x, float y, Color color) {
  // fill in the nearest pixel
  int sx = (int)floor(x);
  int sy = (int)floor(y);

  // check bounds
  if (sx < 0 || sx >= width) return;
  if (sy < 0 || sy >= height) return;


  if (sample_rate == 1) {
      fill_pixel(sx, sy, color);
  }
  else {
      for (int i = 0; i < sample_rate; i++) {
          fill_supersample(sx, sy, i, color);
      }
  }
  return;
}

// Rasterize a line.
void RasterizerImp::rasterize_line(float x0, float y0,
  float x1, float y1,
  Color color) {
  if (x0 > x1) {
    swap(x0, x1); swap(y0, y1);
  }

  float pt[] = { x0,y0 };
  float m = (y1 - y0) / (x1 - x0);
  float dpt[] = { 1,m };
  int steep = abs(m) > 1;
  if (steep) {
    dpt[0] = x1 == x0 ? 0 : 1 / abs(m);
    dpt[1] = x1 == x0 ? (y1 - y0) / abs(y1 - y0) : m / abs(m);
  }

  while (floor(pt[0]) <= floor(x1) && abs(pt[1] - y0) <= abs(y1 - y0)) {
      rasterize_point(pt[0], pt[1], color);
    pt[0] += dpt[0]; pt[1] += dpt[1];
  }
}

// Rasterize a triangle.
void RasterizerImp::rasterize_triangle(float x0, float y0,
                                       float x1, float y1,
                                       float x2, float y2,
                                       Color color) {
  // TODO: Task 1: Implement basic triangle rasterization here, no supersampling
    // changing order of points to be counterclockwise
    // https://algs4.cs.princeton.edu/91primitives/
    if (((x1 - x0) * (y2 - y0) - (x2 - x0) * (y1 - y0)) < 0) {
        swap(x0, x1);
        swap(y0, y1);
    }
    
    // declaring xtest and ytest floats
    float xtest;
    float ytest;
  
    // defining bounding box of triangle, dx and dy floats
    float xmin = min(x0, min(x1, x2));
    float xmax = max(x0, max(x1, x2));
    float ymin = min(y0, min(y1, y2));
    float ymax = max(y0, max(y1, y2));

    float dx0 = x1 - x0;
    float dx1 = x2 - x1;
    float dx2 = x0 - x2;

    float dy0 = y1 - y0;
    float dy1 = y2 - y1;
    float dy2 = y0 - y2;

    if (sample_rate == 1) {
        for (float x = floor(xmin); x <= xmax + 1; x++) {
            for (float y = floor(ymin); y <= ymax + 1; y++) {
                xtest = x + 0.5;
                ytest = y + 0.5;

                if (xtest < 0 || xtest >= width) continue;
                if (ytest < 0 || ytest >= height) continue;

                if ((-(xtest - x0) * dy0 + (ytest - y0) * dx0) >= 0) {
                    if ((-(xtest - x1) * dy1 + (ytest - y1) * dx1) >= 0) {
                        if ((-(xtest - x2) * dy2 + (ytest - y2) * dx2) >= 0) {
                            fill_pixel(xtest, ytest, color);
                        }
                    }
                }
            }
        }
    }
    // TODO: Task 2: Update to implement super-sampled rasterization
    else {
        int dim = sqrt(sample_rate);
        float d0 = (1.0 / dim) / 2; // distance from top left corner of pixel to center of top left subpixel
        float d = (1.0 / dim);      // distance between centers of subpixels
        for (float x = floor(xmin); x <= xmax + 1; x++) {
            for (float y = floor(ymin); y <= ymax + 1; y++) {
                for (int s = 0; s < sample_rate; s++) {
                    xtest = x + d0 + d*(s % dim);
                    ytest = y + d0 + d*floor(s / dim);

                    if (xtest < 0 || xtest >= width) continue;
                    if (ytest < 0 || ytest >= height) continue;

                    if ((-(xtest - x0) * dy0 + (ytest - y0) * dx0) >= 0) {
                        if ((-(xtest - x1) * dy1 + (ytest - y1) * dx1) >= 0) {
                            if ((-(xtest - x2) * dy2 + (ytest - y2) * dx2) >= 0) {
                                fill_supersample(x, y, s, color);
                            }
                        }
                    }
                }
            }
        }
    }
}


void RasterizerImp::rasterize_interpolated_color_triangle(float x0, float y0, Color c0,
                                                          float x1, float y1, Color c1,
                                                          float x2, float y2, Color c2)
{
  // TODO: Task 4: Rasterize the triangle, calculating barycentric coordinates and using them to interpolate vertex colors across the triangle
  // Hint: You can reuse code from rasterize_triangle
    Color color;
    float alpha, beta, gamma;
    
    if (((x1 - x0) * (y2 - y0) - (x2 - x0) * (y1 - y0)) < 0) {
        swap(x0, x1);
        swap(y0, y1);
    }

    // declaring xtest and ytest floats
    float xtest;
    float ytest;

    // defining bounding box of triangle, dx and dy floats
    float xmin = min(x0, min(x1, x2));
    float xmax = max(x0, max(x1, x2));
    float ymin = min(y0, min(y1, y2));
    float ymax = max(y0, max(y1, y2));

    float dx0 = x1 - x0;
    float dx1 = x2 - x1;
    float dx2 = x0 - x2;

    float dy0 = y1 - y0;
    float dy1 = y2 - y1;
    float dy2 = y0 - y2;

    if (sample_rate == 1) {
        for (float x = floor(xmin); x <= xmax + 1; x++) {
            for (float y = floor(ymin); y <= ymax + 1; y++) {
                xtest = x + 0.5;
                ytest = y + 0.5;

                if (xtest < 0 || xtest >= width) continue;
                if (ytest < 0 || ytest >= height) continue;

                if ((-(xtest - x0) * dy0 + (ytest - y0) * dx0) >= 0) {
                    if ((-(xtest - x1) * dy1 + (ytest - y1) * dx1) >= 0) {
                        if ((-(xtest - x2) * dy2 + (ytest - y2) * dx2) >= 0) {
                            alpha = (-(xtest - x1) * (y2 - y1) + (ytest - y1) * (x2 - x1)) / (-(x0 - x1) * (y2 - y1) + (y0 - y1) * (x2 - x1));
                            beta = (-(xtest - x2) * (y0 - y2) + (ytest - y2) * (x0 - x2)) / (-(x1 - x2) * (y0 - y2) + (y1 - y2) * (x0 - x2));
                            gamma = 1 - alpha - beta;
                            color = alpha * c0 + beta * c1 + gamma * c2;
                            fill_pixel(xtest, ytest, color);
                        }
                    }
                }
            }
        }
    }
    else {
        int dim = sqrt(sample_rate);
        float d0 = (1.0 / dim) / 2; // distance from top left corner of pixel to center of top left subpixel
        float d = (1.0 / dim);      // distance between centers of subpixels
        for (float x = floor(xmin); x <= xmax + 1; x++) {
            for (float y = floor(ymin); y <= ymax + 1; y++) {
                for (int s = 0; s < sample_rate; s++) {
                    xtest = x + d0 + d * (s % dim);
                    ytest = y + d0 + d * floor(s / dim);

                    if (xtest < 0 || xtest >= width) continue;
                    if (ytest < 0 || ytest >= height) continue;

                    if ((-(xtest - x0) * dy0 + (ytest - y0) * dx0) >= 0) {
                        if ((-(xtest - x1) * dy1 + (ytest - y1) * dx1) >= 0) {
                            if ((-(xtest - x2) * dy2 + (ytest - y2) * dx2) >= 0) {
                                alpha = (-(xtest - x1) * (y2 - y1) + (ytest - y1) * (x2 - x1)) / (-(x0 - x1) * (y2 - y1) + (y0 - y1) * (x2 - x1));
                                beta = (-(xtest - x2) * (y0 - y2) + (ytest - y2) * (x0 - x2)) / (-(x1 - x2) * (y0 - y2) + (y1 - y2) * (x0 - x2));
                                gamma = 1 - alpha - beta;
                                color = alpha * c0 + beta * c1 + gamma * c2;
                                fill_supersample(x, y, s, color);
                            }
                        }
                    }
                }
            }
        }
    }
}

void RasterizerImp::rasterize_textured_triangle(float x0, float y0, float u0, float v0,
    float x1, float y1, float u1, float v1,
    float x2, float y2, float u2, float v2,
    Texture& tex)
{
    // TODO: Task 5: Fill in the SampleParams struct and pass it to the tex.sample function.
    // TODO: Task 6: Set the correct barycentric differentials in the SampleParams struct.
    // Hint: You can reuse code from rasterize_triangle/rasterize_interpolated_color_triangle

    float alpha, beta, gamma, u, v;
    int level = 1;
    Vector2D uv, u1v, uv1;
    Color c;
    SampleParams sp;

    if (((x1 - x0) * (y2 - y0) - (x2 - x0) * (y1 - y0)) < 0) {
        swap(x0, x1);
        swap(y0, y1);
        swap(u0, u1);
        swap(v0, v1);
    }

    // declaring xtest and ytest floats
    float xtest;
    float ytest;

    // defining bounding box of triangle, dx and dy floats
    float xmin = min(x0, min(x1, x2));
    float xmax = max(x0, max(x1, x2));
    float ymin = min(y0, min(y1, y2));
    float ymax = max(y0, max(y1, y2));

    float dx0 = x1 - x0;
    float dx1 = x2 - x1;
    float dx2 = x0 - x2;

    float dy0 = y1 - y0;
    float dy1 = y2 - y1;
    float dy2 = y0 - y2;

    if (sample_rate == 1) {
        for (float x = floor(xmin); x <= xmax + 1; x++) {
            for (float y = floor(ymin); y <= ymax + 1; y++) {
                xtest = x + 0.5;
                ytest = y + 0.5;

                if (xtest < 0 || xtest >= width) continue;
                if (ytest < 0 || ytest >= height) continue;

                if ((-(xtest - x0) * dy0 + (ytest - y0) * dx0) >= 0) {
                    if ((-(xtest - x1) * dy1 + (ytest - y1) * dx1) >= 0) {
                        if ((-(xtest - x2) * dy2 + (ytest - y2) * dx2) >= 0) {
                            alpha = (-(xtest - x1) * (y2 - y1) + (ytest - y1) * (x2 - x1)) / (-(x0 - x1) * (y2 - y1) + (y0 - y1) * (x2 - x1));
                            beta = (-(xtest - x2) * (y0 - y2) + (ytest - y2) * (x0 - x2)) / (-(x1 - x2) * (y0 - y2) + (y1 - y2) * (x0 - x2));
                            gamma = 1 - alpha - beta;
                            u = alpha * u0 + beta * u1 + gamma * u2;
                            v = alpha * v0 + beta * v1 + gamma * v2;
                            uv = Vector2D(u, v);

                            float xp1 = min(xtest + 1, (float)width-1);
                            float yp1 = min(ytest + 1, (float)height-1);

                            alpha = (-(xp1 - x1) * (y2 - y1) + (ytest - y1) * (x2 - x1)) / (-(x0 - x1) * (y2 - y1) + (y0 - y1) * (x2 - x1));
                            beta = (-(xp1 - x2) * (y0 - y2) + (ytest - y2) * (x0 - x2)) / (-(x1 - x2) * (y0 - y2) + (y1 - y2) * (x0 - x2));
                            gamma = 1 - alpha - beta;
                            u = alpha * u0 + beta * u1 + gamma * u2;
                            v = alpha * v0 + beta * v1 + gamma * v2;
                            u1v = Vector2D(u, v);

                            alpha = (-(xtest - x1) * (y2 - y1) + (yp1 - y1) * (x2 - x1)) / (-(x0 - x1) * (y2 - y1) + (y0 - y1) * (x2 - x1));
                            beta = (-(xtest - x2) * (y0 - y2) + (yp1 - y2) * (x0 - x2)) / (-(x1 - x2) * (y0 - y2) + (y1 - y2) * (x0 - x2));
                            gamma = 1 - alpha - beta;
                            u = alpha * u0 + beta * u1 + gamma * u2;
                            v = alpha * v0 + beta * v1 + gamma * v2;
                            uv1 = Vector2D(u, v);
                            
                            sp.p_uv = uv;
                            sp.p_dx_uv = u1v;
                            sp.p_dy_uv = uv1;
                            sp.lsm = lsm;
                            sp.psm = psm;

                            c = tex.sample(sp);
                            fill_pixel(xtest, ytest, c);
                        }
                    }
                }
            }
        }
    }
    else {
        int dim = sqrt(sample_rate);
        float d0 = (1.0 / dim) / 2; // distance from top left corner of pixel to center of top left subpixel
        float d = (1.0 / dim);      // distance between centers of subpixels
        for (float x = floor(xmin); x <= xmax + 1; x++) {
            for (float y = floor(ymin); y <= ymax + 1; y++) {
                for (int s = 0; s < sample_rate; s++) {
                    xtest = x + d0 + d * (s % dim);
                    ytest = y + d0 + d * floor(s / dim);

                    if (xtest < 0 || xtest >= width) continue;
                    if (ytest < 0 || ytest >= height) continue;

                    if ((-(xtest - x0) * dy0 + (ytest - y0) * dx0) >= 0) {
                        if ((-(xtest - x1) * dy1 + (ytest - y1) * dx1) >= 0) {
                            if ((-(xtest - x2) * dy2 + (ytest - y2) * dx2) >= 0) {
                                alpha = (-(xtest - x1) * (y2 - y1) + (ytest - y1) * (x2 - x1)) / (-(x0 - x1) * (y2 - y1) + (y0 - y1) * (x2 - x1));
                                beta = (-(xtest - x2) * (y0 - y2) + (ytest - y2) * (x0 - x2)) / (-(x1 - x2) * (y0 - y2) + (y1 - y2) * (x0 - x2));
                                gamma = 1 - alpha - beta;
                                u = alpha * u0 + beta * u1 + gamma * u2;
                                v = alpha * v0 + beta * v1 + gamma * v2;
                                uv = Vector2D(u, v);
                                
                                float xp1 = min(xtest + 1, (float)width-1);
                                float yp1 = min(ytest + 1, (float)height-1);

                                alpha = (-(xp1 - x1) * (y2 - y1) + (ytest - y1) * (x2 - x1)) / (-(x0 - x1) * (y2 - y1) + (y0 - y1) * (x2 - x1));
                                beta = (-(xp1 - x2) * (y0 - y2) + (ytest - y2) * (x0 - x2)) / (-(x1 - x2) * (y0 - y2) + (y1 - y2) * (x0 - x2));
                                gamma = 1 - alpha - beta;
                                u = alpha * u0 + beta * u1 + gamma * u2;
                                v = alpha * v0 + beta * v1 + gamma * v2;
                                u1v = Vector2D(u, v);

                                alpha = (-(xtest - x1) * (y2 - y1) + (yp1 - y1) * (x2 - x1)) / (-(x0 - x1) * (y2 - y1) + (y0 - y1) * (x2 - x1));
                                beta = (-(xtest - x2) * (y0 - y2) + (yp1 - y2) * (x0 - x2)) / (-(x1 - x2) * (y0 - y2) + (y1 - y2) * (x0 - x2));
                                gamma = 1 - alpha - beta;
                                u = alpha * u0 + beta * u1 + gamma * u2;
                                v = alpha * v0 + beta * v1 + gamma * v2;
                                uv1 = Vector2D(u, v);

                                sp.p_uv = uv;
                                sp.p_dx_uv = u1v;
                                sp.p_dy_uv = uv1;
                                sp.lsm = lsm;
                                sp.psm = psm;

                                c = tex.sample(sp);
                                fill_supersample(x, y, s, c);
                            }
                        }
                    }
                }
            }
        }
    }

}

void RasterizerImp::set_sample_rate(unsigned int rate) {
  // TODO: Task 2: You may want to update this function for supersampling support

  this->sample_rate = rate;
  supersample_buffer.resize(width * height * sample_rate, Color::White);
}


void RasterizerImp::set_framebuffer_target( unsigned char* rgb_framebuffer,
                                                size_t width, size_t height )
{
  // TODO: Task 2: You may want to update this function for supersampling support

  this->width = width;
  this->height = height;
  this->rgb_framebuffer_target = rgb_framebuffer;
  
}


void RasterizerImp::clear_buffers() {
  // TODO: Task 2: You may want to update this function for supersampling support
  // Hint: With supersampling, you have an additional buffer to take care of

  std::fill(rgb_framebuffer_target, rgb_framebuffer_target + 3 * width * height, 255);
  supersample_buffer.resize(0);
  supersample_buffer.resize(width * height * sample_rate, Color::White);
}


// This function is called at the end of rasterizing all elements of the
// SVG file.  If you use a supersample buffer to rasterize SVG elements
// for antialising, you could use this call to fill the target framebuffer
// pixels from the supersample buffer data.
//
void RasterizerImp::resolve_to_framebuffer() {
  // TODO: Task 2: You will likely want to update this function for supersampling support
    if (sample_rate > 1) {
        float sr;
        float sg;
        float sb;
        Color curr;
        Color meanc;
        for (int x = 0; x < width; x++) {
            for (int y = 0; y < height; y++) {
                sr = 0;
                sg = 0;
                sb = 0;
                for (int s = 0; s < sample_rate; s++) {
                    curr = supersample_buffer[width * y *sample_rate + x*sample_rate + s];
                    sr += curr.r * curr.r;
                    sg += curr.g * curr.g;
                    sb += curr.b * curr.b;
                }
                sr = sqrt(sr / sample_rate);
                sg = sqrt(sg / sample_rate);
                sb = sqrt(sb / sample_rate);
                meanc = Color(sr, sg, sb);
                fill_pixel(x, y, meanc);
            }
        }
    }
}

Rasterizer::~Rasterizer() { }

Vector2D baryouv(float xtest, float ytest, float x0, float y0, float u0, float v0,
    float x1, float y1, float u1, float v1,
    float x2, float y2, float u2, float v2) {
    
    float alpha, beta, gamma, u, v;
    Vector2D uv;

    alpha = (-(xtest - x1) * (y2 - y1) + (ytest - y1) * (x2 - x1)) / (-(x0 - x1) * (y2 - y1) + (y0 - y1) * (x2 - x1));
    beta = (-(xtest - x2) * (y0 - y2) + (ytest - y2) * (x0 - x2)) / (-(x1 - x2) * (y0 - y2) + (y1 - y2) * (x0 - x2));
    gamma = 1 - alpha - beta;
    u = alpha * u0 + beta * u1 + gamma * u2;
    v = alpha * v0 + beta * v1 + gamma * v2;
    uv = Vector2D(u, v);
    return uv;
}

}// CGL
