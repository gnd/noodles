#ifdef GL_ES
precision mediump float;
#endif

uniform float time;
uniform vec2 resolution;
vec2 position;
vec3 color;

float counter(in int m)
{
	float divider = 1000.0 / float(m)*10.0;
	return mod(time*divider, 10.0) / 10.0;
}

vec3 linex(in float xpos, in vec3 col)
{
	return int(mod(position.x*100.0,100.0)) == int(xpos*100.0) ? col : vec3(0.0);
}

vec3 liney(in float ypos, in vec3 col)
{
        return int(mod(position.y*100.0,100.0)) == int(ypos*100.0) ? col : vec3(0.0);
}


void main( void ) {

position = ( gl_FragCoord.xy / resolution.xy );
vec3 col = vec3(1.0,1.0,1.0);
gl_FragColor = vec4(color, 1.0 );

}

