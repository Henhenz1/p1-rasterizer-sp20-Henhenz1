#include "texture.h"
#include "CGL/color.h"

#include <algorithm>
#include <cmath>

namespace CGL {

Color Texture::sample(const SampleParams &sp) {
  // TODO: Task 6: Fill this in.
  // return magenta for invalid level
    float level, prop;
    Color c0, c1;
    try {
        switch (sp.lsm) {
        case 0:
            switch (sp.psm) {
            case 0:
                return sample_nearest(sp.p_uv, 0);
                break;
            case 1:
                return sample_bilinear(sp.p_uv, 0);
                break;
            }
            break;
        case 1:
            level = get_level(sp);
            switch (sp.psm) {
            case 0:
                return sample_nearest(sp.p_uv, (int)level);
                break;
            case 1:
                return sample_bilinear(sp.p_uv, (int)level);
                break;
            }
            break;
        case 2:
            level = get_level(sp);
            switch (sp.psm) {
            case 0:
                c0 = sample_nearest(sp.p_uv, floor(level));
                c1 = sample_nearest(sp.p_uv, ceil(level));
                prop = level - floor(level);
                return prop * c1 + (1 - prop) * c0;
                break;
            case 1:
                c0 = sample_bilinear(sp.p_uv, floor(level));
                c1 = sample_bilinear(sp.p_uv, ceil(level));
                prop = level - floor(level);
                return prop * c1 + (1 - prop) * c0;
                break;
            }
            break;
        }
    }
    // return magenta for invalid level
    catch (...) {
        return Color(1, 0, 1);
    }
}

float Texture::get_level(const SampleParams &sp) {
  // TODO: Task 6: Fill this in.
    Vector2D del_dxuv, del_dyuv;
    float L, D;
    switch (sp.lsm) {
    case 0:
        return 0;
        break;
    case 1:
        del_dxuv = sp.p_dx_uv - sp.p_uv;
        del_dxuv.x *= width;
        del_dxuv.y *= height;
        
        del_dyuv = sp.p_dy_uv - sp.p_uv;
        del_dyuv.x *= width;
        del_dyuv.y *= height;
        
        L = max(sqrt(pow(del_dxuv.x, 2) + pow(del_dxuv.y, 2)), sqrt(pow(del_dyuv.x, 2) + pow(del_dyuv.y, 2)));
        D = round(log2(L));
        return min((float)(mipmap.size() - 1), max((float)0, D));

        break;
    case 2:
        del_dxuv = sp.p_dx_uv - sp.p_uv;
        del_dxuv.x *= width;
        del_dxuv.y *= height;
        
        del_dyuv = sp.p_dy_uv - sp.p_uv;
        del_dyuv.x *= width;
        del_dyuv.y *= height;
        
        L = max(sqrt(pow(del_dxuv.x, 2) + pow(del_dxuv.y, 2)), sqrt(pow(del_dyuv.x, 2) + pow(del_dyuv.y, 2)));        
        D = log2(L);
        
        return min((float)(mipmap.size()-1), max((float)0, D));

        break;
    }
  return 0;
}

Color MipLevel::get_texel(int tx, int ty) {
  return Color(&texels[tx * 3 + ty * width * 3]);
}

Color Texture::sample_nearest(Vector2D uv, int level) {
  // TODO: Task 5: Fill this in.
    try {
        int tx = (int) round(mipmap[level].width * uv.x);
        int ty = (int) round(mipmap[level].height * uv.y);
        if (tx >= mipmap[level].width) {
            tx -= 1;
        }
        if (ty >= mipmap[level].height) {
            ty -= 1;
        }
        return mipmap[level].get_texel(tx, ty);
    } 
  // return magenta for invalid level
    catch (...) {
        return Color(1, 0, 1);
    }
}

Color Texture::sample_bilinear(Vector2D uv, int level) {
  // TODO: Task 5: Fill this in.
    try {
        Vector2D xy = Vector2D(uv.x * mipmap[level].width, uv.y*mipmap[level].height);
        float u0 = floor(xy.x);
        float v0 = floor(xy.y);
        float u1 = ceil(xy.x);
        float v1 = ceil(xy.y);

        if (u0 < 0) {
            u0 += 1;
        }
        else if (u1 >= mipmap[level].width) {
            u1 -= 1;
        }

        if (v0 < 0) {
            v0 += 1;
        }
        else if (v1 >= mipmap[level].height) {
            v1 -= 1;
        }
        
        float s = xy.x - u0;
        float t = xy.y - v0;

        Color u00 = mipmap[level].get_texel(u0, v0);
        Color u10 = mipmap[level].get_texel(u1, v0);
        Color u01 = mipmap[level].get_texel(u0, v1);
        Color u11 = mipmap[level].get_texel(u1, v1);

        // calculating individual RGB values for first lerp
        float r0 = u00.r + (u10.r - u00.r) * s;
        float g0 = u00.g + (u10.g - u00.g) * s;
        float b0 = u00.b + (u10.b - u00.b) * s;
        Color c0 = Color(r0, g0, b0);
        // second lerp
        float r1 = u01.r + (u11.r - u01.r) * s;
        float g1 = u01.g + (u11.g - u01.g) * s;
        float b1 = u01.b + (u11.b - u01.b) * s;
        Color c1 = Color(r1, g1, b1);
        // last lerp
        float r2 = c0.r + (c1.r - c0.r) * t;
        float g2 = c0.g + (c1.g - c0.g) * t;
        float b2 = c0.b + (c1.b - c0.b) * t;
        Color c2 = Color(r2, g2, b2);
        return c2;
    }
    // return magenta for invalid level
    catch (...) {
        return Color(1, 0, 1);
    }
}

Vector2D uvtoxy(float u, float v, float mmw, float mmh) {
    Vector2D* xy = new Vector2D(mmw * u, mmh * v);
    return *xy;
}

/****************************************************************************/

// Helpers

inline void uint8_to_float(float dst[3], unsigned char *src) {
  uint8_t *src_uint8 = (uint8_t *)src;
  dst[0] = src_uint8[0] / 255.f;
  dst[1] = src_uint8[1] / 255.f;
  dst[2] = src_uint8[2] / 255.f;
}

inline void float_to_uint8(unsigned char *dst, float src[3]) {
  uint8_t *dst_uint8 = (uint8_t *)dst;
  dst_uint8[0] = (uint8_t)(255.f * max(0.0f, min(1.0f, src[0])));
  dst_uint8[1] = (uint8_t)(255.f * max(0.0f, min(1.0f, src[1])));
  dst_uint8[2] = (uint8_t)(255.f * max(0.0f, min(1.0f, src[2])));
}

void Texture::generate_mips(int startLevel) {

  // make sure there's a valid texture
  if (startLevel >= mipmap.size()) {
    std::cerr << "Invalid start level";
  }

  // allocate sublevels
  int baseWidth = mipmap[startLevel].width;
  int baseHeight = mipmap[startLevel].height;
  int numSubLevels = (int)(log2f((float)max(baseWidth, baseHeight)));

  numSubLevels = min(numSubLevels, kMaxMipLevels - startLevel - 1);
  mipmap.resize(startLevel + numSubLevels + 1);

  int width = baseWidth;
  int height = baseHeight;
  for (int i = 1; i <= numSubLevels; i++) {

    MipLevel &level = mipmap[startLevel + i];

    // handle odd size texture by rounding down
    width = max(1, width / 2);
    // assert (width > 0);
    height = max(1, height / 2);
    // assert (height > 0);

    level.width = width;
    level.height = height;
    level.texels = vector<unsigned char>(3 * width * height);
  }

  // create mips
  int subLevels = numSubLevels - (startLevel + 1);
  for (int mipLevel = startLevel + 1; mipLevel < startLevel + subLevels + 1;
       mipLevel++) {

    MipLevel &prevLevel = mipmap[mipLevel - 1];
    MipLevel &currLevel = mipmap[mipLevel];

    int prevLevelPitch = prevLevel.width * 3; // 32 bit RGB
    int currLevelPitch = currLevel.width * 3; // 32 bit RGB

    unsigned char *prevLevelMem;
    unsigned char *currLevelMem;

    currLevelMem = (unsigned char *)&currLevel.texels[0];
    prevLevelMem = (unsigned char *)&prevLevel.texels[0];

    float wDecimal, wNorm, wWeight[3];
    int wSupport;
    float hDecimal, hNorm, hWeight[3];
    int hSupport;

    float result[3];
    float input[3];

    // conditional differentiates no rounding case from round down case
    if (prevLevel.width & 1) {
      wSupport = 3;
      wDecimal = 1.0f / (float)currLevel.width;
    } else {
      wSupport = 2;
      wDecimal = 0.0f;
    }

    // conditional differentiates no rounding case from round down case
    if (prevLevel.height & 1) {
      hSupport = 3;
      hDecimal = 1.0f / (float)currLevel.height;
    } else {
      hSupport = 2;
      hDecimal = 0.0f;
    }

    wNorm = 1.0f / (2.0f + wDecimal);
    hNorm = 1.0f / (2.0f + hDecimal);

    // case 1: reduction only in horizontal size (vertical size is 1)
    if (currLevel.height == prevLevel.height) {
      // assert (currLevel.height == 1);

      for (int i = 0; i < currLevel.width; i++) {
        wWeight[0] = wNorm * (1.0f - wDecimal * i);
        wWeight[1] = wNorm * 1.0f;
        wWeight[2] = wNorm * wDecimal * (i + 1);

        result[0] = result[1] = result[2] = 0.0f;

        for (int ii = 0; ii < wSupport; ii++) {
          uint8_to_float(input, prevLevelMem + 3 * (2 * i + ii));
          result[0] += wWeight[ii] * input[0];
          result[1] += wWeight[ii] * input[1];
          result[2] += wWeight[ii] * input[2];
        }

        // convert back to format of the texture
        float_to_uint8(currLevelMem + (3 * i), result);
      }

      // case 2: reduction only in vertical size (horizontal size is 1)
    } else if (currLevel.width == prevLevel.width) {
      // assert (currLevel.width == 1);

      for (int j = 0; j < currLevel.height; j++) {
        hWeight[0] = hNorm * (1.0f - hDecimal * j);
        hWeight[1] = hNorm;
        hWeight[2] = hNorm * hDecimal * (j + 1);

        result[0] = result[1] = result[2] = 0.0f;
        for (int jj = 0; jj < hSupport; jj++) {
          uint8_to_float(input, prevLevelMem + prevLevelPitch * (2 * j + jj));
          result[0] += hWeight[jj] * input[0];
          result[1] += hWeight[jj] * input[1];
          result[2] += hWeight[jj] * input[2];
        }

        // convert back to format of the texture
        float_to_uint8(currLevelMem + (currLevelPitch * j), result);
      }

      // case 3: reduction in both horizontal and vertical size
    } else {

      for (int j = 0; j < currLevel.height; j++) {
        hWeight[0] = hNorm * (1.0f - hDecimal * j);
        hWeight[1] = hNorm;
        hWeight[2] = hNorm * hDecimal * (j + 1);

        for (int i = 0; i < currLevel.width; i++) {
          wWeight[0] = wNorm * (1.0f - wDecimal * i);
          wWeight[1] = wNorm * 1.0f;
          wWeight[2] = wNorm * wDecimal * (i + 1);

          result[0] = result[1] = result[2] = 0.0f;

          // convolve source image with a trapezoidal filter.
          // in the case of no rounding this is just a box filter of width 2.
          // in the general case, the support region is 3x3.
          for (int jj = 0; jj < hSupport; jj++)
            for (int ii = 0; ii < wSupport; ii++) {
              float weight = hWeight[jj] * wWeight[ii];
              uint8_to_float(input, prevLevelMem +
                                        prevLevelPitch * (2 * j + jj) +
                                        3 * (2 * i + ii));
              result[0] += weight * input[0];
              result[1] += weight * input[1];
              result[2] += weight * input[2];
            }

          // convert back to format of the texture
          float_to_uint8(currLevelMem + currLevelPitch * j + 3 * i, result);
        }
      }
    }
  }
}

} // namespace CGL
