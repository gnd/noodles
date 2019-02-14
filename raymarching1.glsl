#ifdef GL_ES
precision lowp float;
#endif

uniform float time;
uniform vec2 resolution;

// counter 0 to 1 
float cnt(in int m)
{
	float divider = 1000.0 / float(m)*10.0;
	return mod(time*divider, 10.0) / 10.0;
}


void main(void) {
vec3 c;
 vec4 p_ray = vec4(gl_FragCoord.xy,0.,1.)/resolution.xyxy-0.4;
 vec4 d=p_ray;
 vec4 t;
    p_ray.z += time;
    for(float i=1.8; i>0.; i-=.1)
    {
        t = abs(mod(p_ray, 0.5)-.25);
        float x = max(t.y*2., length(t.xz)*2.);
        c = vec3(i);
        if(x<.2) break;
        p_ray -= d*x;
     }
	


c+=vec3(cnt(100),0.,0.);

gl_FragColor = vec4(c, 1.0);
}
