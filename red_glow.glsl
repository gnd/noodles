uniform vec2 resolution;
uniform float time;
#define MAXSTEPS 32
#define MAXDIST 4.0

float sphere( vec3 p, float r ) {
    return length(p) - r;
}

float maxcomp(in vec3 p ) {
    return max(p.x,max(p.y,p.z));
}

float box(vec3 p, vec3 b, float r) {
    vec3 d = abs(p) - b;
    return min(maxcomp(d),0.0) - r + length(max(d,0.0));
}

void main() {
    // setup some globals
    vec3 eye = vec3(cos(time), sin(time), .1);
    vec3 up = vec3(0, 1, 0);
    vec3 right = vec3(1, 0, 0);
    float u = gl_FragCoord.x * 2.0 / resolution.x - 1.0;
    float v = gl_FragCoord.y * 2.0 / resolution.y - 1.0;
    float f = length(eye);
    //vec3 rd = normalize(cross(right, up));
    vec3 rd = vec3(u,v,.0) - eye;
    vec3 ro = eye + f * rd + right * u + up * v/1.77;
    vec4 color = vec4(0.0);

    // raymarch a sphere
    float t = .0;
    for (int i = 0; i < MAXSTEPS; ++i) {
        vec3 p = ro + rd * t;
        float d = min(box(p-vec3(.0,-.5,.0), vec3(2.,.01,5.), .0) ,min(box(p, vec3(.1,.1,3.9), .0), min(box(p, vec3(0.7,0.2,0.1),0.),sphere(p,.333))));
        if (d > MAXDIST) {
            color = vec4(float(i)*(1./32.), 0.,0.,0.);
            break;
        }
        if (d < 0.01) {
            color = vec4(float(i)*(1./32.), 0.,0., 0.);
            break;
        }
        t += d;
        color = vec4(float(i)*(1./32.), 0.,0., 0.);
    }

    gl_FragColor = color;
}
