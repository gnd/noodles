// current_texture_image - the texture that has the video image
//   initialized by the CPU program to contain the number of the 
//   texture unit that contains the video texture.
//   the data indicated by this is assumed to be in the form of
//   a luminance/alpha texture (two components per pixel). in the
//   shader we extract those components from each pixel and turn
//   them into RGB to draw the picture.

// texture_width - the width of the base texture for the window
//   (that is, not the size of the window or the size of the image
//   being textured onto it; but that power-of-two size texture width.)

// texel_width - the width of a texel in texture coordinates. 
//   the texture coordinates in the shader go from 0 to 1.0.
//   so each texel is (1.0 / texture_width) wide.

// texture_height - the height of the base texture for the window
//   (that is, not the size of the window or the size of the image
//   being textured onto it; but that power-of-two size texture height.)

// image_height - the height of the image being put onto the texture.
//   for instance texture_height might be 512 (a power of two) and 
//   image height might be 480.

// shader_on - a 1/0 flag for whether or not to use the YUV->RGB
//   translation code. if this is zero, fragment gets the
//   regular opengl color. 

// color_output - a 1/0 flag; if non-zero generate output in color.
// if zero, generate greyscale

//uniform float texture_height;
//uniform float image_height;
uniform sampler2D video;
float texel_width;
float texture_width;

void main()
{
  float red, green, blue;
  vec4 luma_chroma;
  float luma, chroma_u,  chroma_v;
  float pixelx, pixely;
  float xcoord, ycoord;
  vec3 yuv;
  
  // note: pixelx, pixely are 0.0 to 1.0 so "next pixel horizontally"
  //  is not just pixelx + 1; rather pixelx + texel_width.

  pixelx = gl_TexCoord[0].x;
  pixely = gl_TexCoord[0].y;

  // if pixelx is even, then that pixel contains [Y U] and the 
  //    next one contains [Y V] -- and we want the V part.
  // if  pixelx is odd then that pixel contains [Y V] and the 
  //     previous one contains  [Y U] -- and we want the U part.
  // note: only valid for images whose width is an even number of
  // pixels
	texture_width = 320.0;
	texel_width = 1.;
  xcoord = floor (pixelx * texture_width);
  luma_chroma = texture2D(video, p);
      
  // just look up the brightness
  luma = (luma_chroma.r - 0.0625) * 1.1643;
	  
  chroma_u = luma_chroma.g;
  chroma_v = texture2D(video, vec2(pixelx + texel_width, pixely)).a; 
  chroma_u = chroma_u - 0.5;
  chroma_v = chroma_v - 0.5;
  
  red = luma + 1.5958 * chroma_v;
  green = luma - 0.39173 * chroma_u - 0.81290 * chroma_v;
  blue = luma + 2.017 * chroma_u;

  // set the color based on the texture color
  gl_FragColor = vec4(red, green, blue, 1.0);
}
