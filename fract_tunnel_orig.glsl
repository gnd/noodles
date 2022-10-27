#version 330
out vec4 PixelColor;
uniform vec2 resolution;
uniform float time;
#define eps 0.0001
#define MAXSTEPS 128
#define MAXDIST 12

// http://glslsandbox.com/e#39898.3

struct light { vec3 position; vec3 color; };
struct material { vec3 color; float reflection_ratio; float shininess; };
struct shading { float diffuse; float specular; float shadow; float aoc; float amb; };

light l1;
material m1,m2,m3;
shading s1,s2,s3;



void main()
{
    // classic scene instantiate
    vec3 ro = vec3(3., 2.5, 3.);
    // ro = vec3(cos(time/4.)*3., 2.5, sin(time/4.)*3.);
    vec3 lookat = vec3(0.,0.,0.);
    vec3 fwd = normalize(lookat-ro);
    vec3 right = normalize(vec3(fwd.z, 0., -fwd.x));
    vec3 up = normalize(cross(fwd, right));
    vec2 uv = gl_FragCoord.xy * 2.0 / resolution - 1.0;
    float aspect = resolution.x/resolution.y;
    vec3 rd = normalize(1.4*fwd + uv.x*right*aspect + uv.y*up);



    vec4 p = vec4(gl_FragCoord.xy,0.,1.) / resolution.x-.8;
    //vec4 p = vec4(uv, 0.,1.);
    vec4 r = vec4(.4);
    vec4 q=r+.99;
    p.y+=.51;
    q.zw-=time*.1;
    
    for (int i = 0; i < MAXSTEPS; ++i) {
        float d = 0;
        float s = .3;

        for (int j = 0; j < 7; j++)
            r = max(r = abs(mod(q*s+1.3,1.3)-.9), r.yzxw),
            d = max(d,(0.22-length(r*0.51)*.35)/s),
            s *= 2.9;

        q += p*d;
        
        PixelColor = vec4(0.01 * (MAXSTEPS-i));

        if (d < eps) break;
    }
    PixelColor.a = 1.;
}
