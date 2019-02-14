#ifdef GL_ES
precision lowp float;
#endif

uniform float time;
uniform vec2 resolution;

void main(void) {
 vec4 p_ray = vec4(gl_FragCoord.xy,0.,1.)/resolution.xyxy-0.4;
 vec4 d=p_ray;
 vec4 t;
    p_ray.z += time;
    for(float i=1.8; i>0.; i-=.1)
    {
        t = abs(mod(p_ray, 0.5)-.25);
        float x = max(t.y*2., length(t.xz)*2.);
        gl_FragColor = vec4(i);
        if(x<.2) break;
        p_ray -= d*x;
     }
	
}
