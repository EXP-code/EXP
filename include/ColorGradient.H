#ifndef _ColorGradient_H
#define _ColorGradient_H

#include <png++/png.hpp>

//! Internal class used to store colors at different points in the gradient.
class ColorGradient
{
private:
  struct ColorPoint
  {
    //! Red, green and blue values of our color.
    float r, g, b;
    //! Position of our color along the gradient (between 0 and 1).
    float val;
    ColorPoint(float red, float green, float blue, float value)
      : r(red), g(green), b(blue), val(value) {}
  };

  //! An array of color points in ascending value.
  std::vector<ColorPoint> color;
 
public:
  //! Default constructor:
  ColorGradient()  {  createDefaultGradient();  }
 
  //! Inserts a new color point into its correct position:
  void addColorPoint(float red, float green, float blue, float value)
  {
    for(int i=0; i<color.size(); i++)  {
      if(value < color[i].val) {
        color.insert(color.begin()+i, ColorPoint(red,green,blue, value));
        return;  }}
    color.push_back(ColorPoint(red,green,blue, value));
  }
 
  //! Inserts a new color point into its correct position:
  void clearGradient() { color.clear(); }
 
  //! Default gradient
  void createDefaultGradient()
  {
    createGrayGradient();
  }
  
  //! 7 color heapmap gradient into the "color" vector:
  void createGrayGradient()
  {
    color.clear();
    color.push_back(ColorPoint(0, 0, 0,   0.0000f)); // Black.
    color.push_back(ColorPoint(1, 1, 1,   1.0000f)); // White.
  }

  //! 7 color heapmap gradient into the "color" vector:
  void createSevenColorHeatMapGradient()
  {
    color.clear();
    color.push_back(ColorPoint(0, 0, 0,   0.0000f)); // Black.
    color.push_back(ColorPoint(0, 0, 1,   0.1667f)); // Blue.
    color.push_back(ColorPoint(0, 1, 1,   0.3333f)); // Cyan.
    color.push_back(ColorPoint(0, 1, 0,   0.5000f)); // Green.
    color.push_back(ColorPoint(1, 1, 0,   0.6667f)); // Yellow.
    color.push_back(ColorPoint(1, 0, 0,   0.8333f)); // Red.
    color.push_back(ColorPoint(1, 1, 1,   1.0000f)); // White.
  }

  //! 5 color heapmap gradient into the "color" vector:
  void createFiveColorHeatMapGradient()
  {
    color.clear();
    color.push_back(ColorPoint(0, 0, 1,   0.0f));  // Blue.
    color.push_back(ColorPoint(0, 1, 1,   0.25f)); // Cyan.
    color.push_back(ColorPoint(0, 1, 0,   0.5f));  // Green.
    color.push_back(ColorPoint(1, 1, 0,   0.75f)); // Yellow.
    color.push_back(ColorPoint(1, 0, 0,   1.0f));  // Red.
  }
 
  //! Inputs a (value) between 0 and 1 and outputs the (red), (green) and (blue)
  //! values representing that position in the gradient.
  png::rgb_pixel operator()(const float value)
  {
    if (color.size()==0)
      return png::rgb_pixel(0, 0, 0);
    
    float red, green, blue;

    for (int i=0; i<color.size(); i++) {
      
      ColorPoint &currC = color[i];
      if (value < currC.val) {
        ColorPoint &prevC  = color[ std::max<int>(0,i-1) ];
        float valueDiff    = (prevC.val - currC.val);
        float fractBetween = (valueDiff==0) ? 0 : (value - currC.val) / valueDiff;
        red   = (prevC.r - currC.r)*fractBetween + currC.r;
        green = (prevC.g - currC.g)*fractBetween + currC.g;
        blue  = (prevC.b - currC.b)*fractBetween + currC.b;

	return png::rgb_pixel(floor(red  *255.99),
			      floor(green*255.99),
			      floor(blue *255.99));
      }
    }
    red   = color.back().r;
    green = color.back().g;
    blue  = color.back().b;

    return png::rgb_pixel(floor(red  *255.99),
			  floor(green*255.99),
			  floor(blue *255.99));
  }
};

#endif
