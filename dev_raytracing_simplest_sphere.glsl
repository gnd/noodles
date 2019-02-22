// do you own lightning
// write down the lightning principle (fake phong)

uniform float time;
uniform vec2 resolution;
uniform vec4 mouse;
uniform vec2 rand;

const float g_rmEpsilon = 0.005;
vec3 p,ro;

float sph(in vec3 p) {
	vec3 cen= vec3(0.,0.,-0.1);
	return length(p-cen) - 0.6;
}

float repsph(in vec3 p) {
        //p.x = mod(p.x*3., 1.) -0.3;
        //p.z = mod(p.z*3., 1.) -0.3;
        return sph(p);
}

vec3 normal(in vec3 p) {
	float x = sph(vec3(p.x+g_rmEpsilon,p.y,p.z)) - sph(vec3(p.x-g_rmEpsilon,p.y,p.z));
	float y = sph(vec3(p.x,p.y+g_rmEpsilon,p.z)) - sph(vec3(p.x,p.y-g_rmEpsilon,p.z));
	float z = sph(vec3(p.x,p.y,p.z+g_rmEpsilon)) - sph(vec3(p.x,p.y,p.z-g_rmEpsilon));
	return normalize(vec3(x,y,z));
}


void main()
{
    vec3 eye = vec3(0, 0, -1);
    vec3 up = vec3(0, 1, 0);
    vec3 right = vec3(1, 0, 0)*(resolution.x/resolution.y);

    float u = gl_FragCoord.x * 2.0 / resolution.x - 1.0;
    float v = gl_FragCoord.y * 2.0 / resolution.y - 1.0;
    ro = right * u + up * v;
    vec3 rd = normalize(cross(right, up));

    vec4 color = vec4(0.0); // Sky color

    float t = 0.0;
    const int maxSteps = 32;
    for(int i = 0; i < maxSteps; ++i)
    {
        p = ro + rd * t;
        float d = repsph(p); // Distance to sphere of radius 0.5
        if(d < g_rmEpsilon)
        {
            color = vec4(1.0); // Sphere color
            break;
        }

        t += d;
    }

    // lightning
    vec3 lightpos = vec3(sin(time/3.)*1.,cos(time/3.)*1.,1.);
    vec3 surfnorm = normal(p);

    float light = max(dot(surfnorm, normalize(lightpos-p)),0.);

	gl_FragColor = color*light;
}
