#version 330
out vec4 PixelColor;
uniform vec2 resolution;
uniform float time;
#define eps 0.0001

#define ITER 18
#define EPS 0.2
#define NEAR 1.
#define FAR 50.

float map(vec3 p){
    p=floor(p*4.);
    return fract(sin(p.x*1.3)*13.+cos(p.y*1.7)*17.+tan(p.z*19.9)*19.)*11.;
}

float trace(vec3 ro,vec3 rd) {
    float t=NEAR,d;
    for (int i=0;i<ITER;i++) {
        d=map(ro+rd*t);
        if(abs(d)<EPS||t>FAR) break;
        t+=step(d,1.)*d*0.2+d*.5;
    }
    return min(t,FAR);
}

void main(void){
    vec2 uv = gl_FragCoord.xy * 2.0 / resolution - 1.0;
    float v=1.-trace(vec3(sin(time*.1),1.+time*20.,time*2.),vec3(uv,0.02))/FAR;
    PixelColor = vec4(vec3(2.,.2,20.1)*v,1);
}