/*
 * "Frozen Lake" by Alexander Alekseev aka TDM - 2019
 * License Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.
 * Contact: tdmaav@gmail.com
 */

const float THRESHOLD 	= 0.001;
const float EPSILON 	= 5e-3;

const float HEIGHT_POWER = 5.0;
const float CRACKS_SCALE = 0.6;
const float CRACKS_THICKNESS = 0.9;
const float CRACKS_ALPHA = 0.8;
const float REFRACTION = 0.8;
const float BUBBLES_BRIGHTNESS = 0.4;
const vec3 SNOW_COLOR = vec3(0.85,0.98,1.0);
const vec3 FOG_COLOR = vec3(0.0,0.04,0.05);
const vec3 DEEP_COLOR = vec3(0.0,0.12,0.2);
const vec3 CRACKS_COLOR = vec3(0.3,0.95,1.0) * 1.2;
const vec3 CRACKS_COLOR_TOP = vec3(1.6);
const vec3 MOUNTAINS_COLOR = vec3(0.04,0.02,0.0);

// tracing

float mapCracks1(vec3 p) {
    const float SCALE = 0.1;
    p.x += sin(p.z*0.2) * 2.0;
    p.x += triangle(p.z * 0.053) * 2.0;
    p.z += triangle(p.x * 0.103) * 2.0;
    return voronoi(p.xz*SCALE).x / SCALE * 0.9;
}

float mapCracks2(vec3 p) {
    const float SCALE = 0.25;
    p.x += triangle(p.z * 0.153) * 1.5;
    p.z += triangle(p.x * 0.203) * 1.5;
    return voronoi(p.xz*SCALE).x / SCALE * 0.9;
}

vec2 traceCracks1(vec3 ori, vec3 dir, out vec3 p) {
    float t = 0.0;
    float d = 0.0;
    for(int i = 0; i < 10; i++) {
        p = ori + dir * t;
        d = mapCracks1(p);
        if(d < THRESHOLD) break;
        t += d * 0.9;
    } 
    return vec2(d,t);
}
vec2 traceCracks2(vec3 ori, vec3 dir, float s, out vec3 p) {
    float t = 0.0;
    float d = 0.0;
    for(int i = 0; i < 8; i++) {
        p = ori + dir * t;
        d = mapCracks2(p*s);
        if(d < THRESHOLD) break;
        t += d * 0.9;
    } 
    return vec2(d,t);
}
vec2 traceCracks3(vec3 ori, vec3 dir, out vec3 p) {
    float t = 0.0;
    float d = 0.0;
    for(int i = 0; i < 3; i++) {
        p = ori + dir * t;
        d = mapCracks1(p*0.7);
        if(d < THRESHOLD) break;
        t += d;
    } 
    return vec2(d,t);
}
vec2 getNormalCracks1(vec3 p) {
    float t = mapCracks1(p);
    vec2 n;
    n.x = mapCracks1(vec3(p.x+EPSILON,p.y,p.z)) - t;
    n.y = mapCracks1(vec3(p.x,p.y,p.z+EPSILON)) - t;
    return normalize(n);
}
vec2 getNormalCracks2(vec3 p) {
    float t = mapCracks2(p);
    vec2 n;
    n.x = mapCracks2(vec3(p.x+EPSILON,p.y,p.z)) - t;
    n.y = mapCracks2(vec3(p.x,p.y,p.z+EPSILON)) - t;
    return normalize(n);
}


/*
 * color
 */

// sky
vec3 getSkyColor(vec3 e, bool isReflection) {
    e.y = max(e.y,0.0);
    float yy = pow(e.y, 0.9);
    vec3 ret;
    ret.x = pow(1.0-yy-0.05,8.0) * 0.75;
    ret.y = pow(1.0-yy, 4.0) * 0.75;
    ret.z = pow(1.0-yy,2.0);
       
    
    float phi = atan(e.z,e.x) / PI;
    float h = (fbm1(phi*10.0)*0.5+0.5)*0.14-0.03;
    float mountains = isReflection ? 
        smoothstep(h+0.05,h-0.01,e.y) :
    	smoothstep(h+0.002,h,e.y);
    ret = mix(ret,MOUNTAINS_COLOR,
              mountains*(pow(e.y,0.3) * 0.15 + 0.85));
    
    h = (fbm1(phi*14.0)*0.5+0.5)*0.1-0.01;
    float mf = isReflection ? 
        smoothstep(h+0.05,h-0.01,e.y) :
    	smoothstep(h+0.002,h,e.y);
    ret = mix(ret,MOUNTAINS_COLOR,
              mf*(pow(e.y,0.5) * 0.5 + 0.5)*0.8*(1.0-mountains));
    
    
    // clouds
    vec3 p;
    intersectionPlane(vec3(0.0,300.0,0.0),e,p);
    ret = mix(ret,vec3(1.0), fbmClouds(p.xz)*(1.0-mountains)*(1.0-mf) * 0.7);
    
    //return vec3(clamp(phi,0.0,1.0));
    return ret;
}

// snow

float getSnowWindMask(in vec2 p, float t) {
    float amp = 0.5;
    float frq = 1.0;
    float wrt = t*2.0;
    p.x += sin(frq*p.y + wrt*0.9) * amp;
    p.y += cos(frq*p.x*1.5 + wrt*0.8) * amp;
    p.x += sin(frq*p.y*1.9 + wrt*0.7) * amp;
    p.y += cos(frq*p.x*1.7 + wrt*0.6) * amp;
    
    float wind = fbm2(p,t*8.0);
    wind = wind * 0.5 + 0.5;
    return wind * wind;
}

float getSnowMask(in vec2 p) {
    mat2 m = mat2( 1.6,  1.2, -1.2,  1.6 );
    float a = 1.0;
    float w = 1.0;
    float f = noise12( p );
    for(int i = 0; i < 6; i++) {
        p = m * p; a /= 1.5;
        f += a * (abs(noise12( p )));
        w += a;
    }
    f /= w;
    f = smoothstep(0.55,0.65,f);
    f = pow(f,0.25);
    f = f * 0.9;
    
#ifdef WIND
    p *= 0.02;
    p.x *= 0.5;
    float wind = getSnowWindMask(p,iTime * 0.25);
    wind = max(wind, getSnowWindMask(p,iTime * 0.25 + 1.5));
    wind *= wind;
    wind *= 0.5;    
    return max(f, wind);
#else
    return f;
#endif
}

vec3 getObjectColor(in vec3 p, const in vec3 cam, in vec3 e) {
    vec3 op = p;
    vec3 dir = e;
    const vec3 n = vec3(0.0,1.0,0.0);
    float depth = length(p - cam);
    float depth_f = max(depth*0.8, 1.0);
    p *= CRACKS_SCALE;
    
    
    // global thickness modulation
    float gth = 0.6 + 0.8 * smoothstep(0.2,0.8, noise13(p*0.05));    
    gth *= CRACKS_THICKNESS;
    
    // crack depth
    vec3 cp;
    vec3 norm = vec3(1.0,noise2(p.xz*3.)*0.2);
    norm.yz += noise2(p.xz*10.)*0.2;
    norm.x *= depth_f;    
    norm = normalize(norm.yxz);
    e.xz += norm.xz * REFRACTION;
        
    traceCracks1(p,e,cp);
    vec2 cr1_normal = getNormalCracks1(cp);
    float crack_depth = abs(cp.y - p.y);
    crack_depth = pow(max(1.0-crack_depth*0.2/gth, 0.0),HEIGHT_POWER) * 0.6;
    crack_depth *= 0.5 + 0.5 * noise13(cp*vec3(0.7,10.0,0.7));
    crack_depth *= abs(cr1_normal.x) * 0.6 + 0.4;
    
    traceCracks2(p,e,1.0,cp);
    vec2 cr2_normal = getNormalCracks2(cp);
    float crack_depth_2 = abs(cp.y - p.y);
    crack_depth_2 = pow(max(1.0 - crack_depth_2 * 0.4/gth, 0.0), HEIGHT_POWER) * 0.6;
    crack_depth_2 *= 0.5 + 0.5 * smoothstep(0.2,0.9, noise13(cp*vec3(12.0,1.0,12.0)));
    crack_depth_2 *= 0.5 + 0.5 * noise13(cp*vec3(1.0,20.0,1.0));
    crack_depth_2 *= abs(cr2_normal.x) * 0.6 + 0.4;
    
    traceCracks2(p,e,1.5,cp);
    float crack_depth_3 = abs(cp.y - p.y);
    crack_depth_3 = pow(max(1.0 - crack_depth_3 * 3.0/gth , 0.0), HEIGHT_POWER) * 0.3;
    crack_depth_3 *= 0.5 + 0.5 * smoothstep(0.3,0.9, noise13(cp*vec3(17.0,1.0,17.0)));
      
    vec2 c4n = noise2(p.xz*30.0) * 0.4;
    traceCracks3(p,e+c4n.xxy,cp);
    float crack_depth_4 = abs(cp.y - p.y + 2.0);
    crack_depth_4 = pow(max(1.0-crack_depth_4*0.2/gth, 0.0),3.0) * 0.15;
    crack_depth_4 *= 0.5 + 0.5 * noise13(cp*vec3(0.7,10.0,0.7));
    
    
    
    // base color
    vec3 col = toLinear(DEEP_COLOR);
    
    // bubbles    
    dir.xz += norm.xz * REFRACTION * 0.3;
    vec3 bp;
    intersectionPlane(cam+vec3(0.,0.5,0.),dir,bp);        
    col += pow(noise13(bp * 14.0),20.0) * BUBBLES_BRIGHTNESS * gth;
    intersectionPlane(cam+vec3(0.,1.,0.),dir,bp);        
    col += pow(noise13(bp * 15.0),20.0) * BUBBLES_BRIGHTNESS * gth;
    intersectionPlane(cam+vec3(0.,2.,0.),dir,bp);        
    col += pow(noise13(bp * 16.0),20.0) * BUBBLES_BRIGHTNESS * gth;
    
    // cracks color
    vec3 crc = toLinear(CRACKS_COLOR);
    vec3 crct = toLinear(CRACKS_COLOR_TOP);
    float a = 0.4 + 0.6 * smoothstep(0.2,0.8, noise13(p*0.07));
    a *= CRACKS_ALPHA;
    col = mix(col, mix(crc,crct,crack_depth_4), 
              crack_depth_4 * a);
    col = mix(col, mix(crc,crct,crack_depth_3), 
              crack_depth_3 * a);
    col = mix(col, mix(crc,crct,crack_depth_2), 
              crack_depth_2 * a);
    col = mix(col, mix(crc,crct,crack_depth), 
              crack_depth * a);
        
    // reflection
    float fresnel = pow(max(1.0 - dot(-e,n),0.0),5.0) * 0.9 + 0.1;
    vec3 rdir = reflect(e,norm);
    vec3 reflection = getSkyColor(rdir,true);
   
    col = mix(col,reflection,fresnel);
    
    // snow surface
    depth_f = max(depth*0.01, 1.0);
    float snow = getSnowMask(p.xz*0.1) / depth_f;
    col = mix(col,SNOW_COLOR,snow);
    return col;
}



// main
void mainImage( out vec4 fragColor, in vec2 fragCoord ) {
	vec2 iuv = fragCoord.xy / iResolution.xy * 2.0 - 1.0;
    vec2 uv = iuv;
    uv.x *= iResolution.x / iResolution.y;    
    vec2 mouse = iMouse.xy / iResolution.xy * 4.0 - 2.0;
        
    // ray
    float xst = iTime * 0.3;
    vec3 ori = vec3(sin(xst)*10.0,5.0,-iTime*5.0);
    
    vec3 ang = vec3(cos(xst+1.0)*0.1,PI*0.1, cos(xst)*0.4);
    if(iMouse.z > 0.0) ang = vec3(0.0,clamp(2.0-iMouse.y*0.01,-0.3,PI),iMouse.x*0.01);
	mat3 rot = fromEuler(ang);
    
    vec3 dir = normalize(vec3(uv.xy,-1.5)); dir.z += length(uv) * 0.05;   
    dir = normalize(dir * rot);
             
    // color
    vec3 p;
    vec3 color = getSkyColor(dir,false);
    if(intersectionPlane(ori,dir,p))
        color = getObjectColor(p,ori,dir);
    
    // post
    //color *= 1.3;
    color = pow(color,vec3(0.4545));
    
    // vignette
    float vgn = smoothstep(1.2,0.5,abs(iuv.y)) * smoothstep(1.2,0.5,abs(iuv.x));
    color *= vgn * 0.3 + 0.7;
               
	fragColor = vec4(color,1.0);
}



/*
 * "Frozen Lake" by Alexander Alekseev aka TDM - 2019
 * License Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.
 * Contact: tdmaav@gmail.com
 */

#define HASHSCALE3 vec3(.1031, .1030, .0973)
const float PI = 3.141592;

float saturate(float x) { return clamp(x,0.0,1.0); }
float mul(vec2 x) { return x.x*x.y; }
mat3 fromEuler(vec3 ang) {
	vec2 a1 = vec2(sin(ang.x),cos(ang.x));
    vec2 a2 = vec2(sin(ang.y),cos(ang.y));
    vec2 a3 = vec2(sin(ang.z),cos(ang.z));
    mat3 m;
    m[0] = vec3(a1.y*a3.y+a1.x*a2.x*a3.x,a1.y*a2.x*a3.x+a3.y*a1.x,-a2.y*a3.x);
	m[1] = vec3(-a2.y*a1.x,a1.y*a2.y,a2.x);
	m[2] = vec3(a3.y*a1.x*a2.x+a1.y*a3.x,a1.x*a3.x-a1.y*a3.y*a2.x,a2.y*a3.y);
	return m;
}

bool intersectionPlane(vec3 o, vec3 d, out vec3 p) {
    float t = o.y / d.y;
    p = o - d * t;
    return bool(step(t,0.0));
}

float hash11(float x) {
    return fract(sin(x) * 43758.5453);
}
float hash12( vec2 p ) {
	float h = dot(p,vec2(127.1,311.7));	
    return fract(sin(h)*43758.5453123);
}
vec2 hash22(vec2 p) {
	vec3 p3 = fract(vec3(p.xyx) * HASHSCALE3);
    p3 += dot(p3, p3.yzx+19.19);
    return fract((p3.xx+p3.yz)*p3.zy);
}
float hash13(in vec3 p) {
    p  = fract( p*0.3183099+.1 );
	p *= 17.0;
    return fract( p.x*p.y*p.z*(p.x+p.y+p.z) );
}
float noise11(in float p) {
    float i = floor( p );
    float f = fract( p );	
	float u = f*f*(3.0-2.0*f);
    return -1.0+2.0*mix(hash11(i),hash11(i+1.0),u);
}
float noise12( in vec2 p ) {
    vec2 i = floor( p );
    vec2 f = fract( p );	
	vec2 u = f*f*(3.0-2.0*f);
    return -1.0+2.0*mix( mix( hash12( i + vec2(0.0,0.0) ), 
                     hash12( i + vec2(1.0,0.0) ), u.x),
                mix( hash12( i + vec2(0.0,1.0) ), 
                     hash12( i + vec2(1.0,1.0) ), u.x), u.y);
}
vec2 noise2( in vec2 p ) {
    vec2 i = floor( p );
    vec2 f = fract( p );	
	vec2 u = f*f*(3.0-2.0*f);
    return -1.0+2.0*mix( mix( hash22( i + vec2(0.0,0.0) ), 
                     hash22( i + vec2(1.0,0.0) ), u.x),
                mix( hash22( i + vec2(0.0,1.0) ), 
                     hash22( i + vec2(1.0,1.0) ), u.x), u.y);
}

float noise13(in vec3 p) {
    vec3 i = floor( p );
    vec3 f = fract( p );	
	vec3 u = f*f*(3.0-2.0*f);
    
    float a = hash13( i + vec3(0.0,0.0,0.0) );
	float b = hash13( i + vec3(1.0,0.0,0.0) );    
    float c = hash13( i + vec3(0.0,1.0,0.0) );
	float d = hash13( i + vec3(1.0,1.0,0.0) ); 
    float v1 = mix(mix(a,b,u.x), mix(c,d,u.x), u.y);
    
    a = hash13( i + vec3(0.0,0.0,1.0) );
	b = hash13( i + vec3(1.0,0.0,1.0) );    
    c = hash13( i + vec3(0.0,1.0,1.0) );
	d = hash13( i + vec3(1.0,1.0,1.0) );
    float v2 = mix(mix(a,b,u.x), mix(c,d,u.x), u.y);
        
    return abs(mix(v1,v2,u.z));
}

float fbm1(in float p) {
    float m = 2.0;
    float a = 1.0;
    float w = 1.0;
    float f = noise11( p );
    for(int i = 0; i < 8; i++) {
        p *= m; a /= 1.8;
        f += a*noise11( p );
        w += a;
    }
    return f / w;
}

float fbm2(in vec2 p, float t) {
    float m = 2.0;
    float a = 1.0;
    float w = 1.0;
    float f = noise12( p );
    for(int i = 0; i < 8; i++) {
        p *= m; a /= 1.5;
        f += a*noise12( p+t );
        w += a;
    }
    return f / w;
}

vec2 fbm22(in vec2 p) {
    float m = 2.0;
    float a = 1.0;
    float w = 1.0;
    vec2 f = noise2( p );
    for(int i = 0; i < 8; i++) {
        p *= m; a /= 1.2;
        f += a*noise2(p);
        w += a;
    }
    return f / w;
}

float fbmClouds(in vec2 p) {
    p *= 0.001;
    float m = 2.0;
    float a = 1.0;
    float w = 1.0;
    float f = noise12( p );
    for(int i = 0; i < 4; i++) {
        p *= m; a /= 1.5;
        f += a* abs(noise12( p ));
        w += a;
    }
    f /= w;
    //f = pow(max(f,0.0001),5.0);
    f = max((f - 0.4) / (1.0 - 0.4), 1e-4);
    f = sqrt(f);
    return f;
}

// iq's voronoi
vec3 voronoi( in vec2 x ) {
    vec2 n = floor(x);
    vec2 f = fract(x);

    //----------------------------------
    // first pass: regular voronoi
    //----------------------------------
	vec2 mg, mr;

    float md = 8.0;
    for( int j=-1; j<=1; j++ )
    for( int i=-1; i<=1; i++ )
    {
        vec2 g = vec2(float(i),float(j));
		vec2 o = hash22( n + g );
        vec2 r = g + o - f;
        float d = dot(r,r);

        if( d<md )
        {
            md = d;
            mr = r;
            mg = g;
        }
    }

    //----------------------------------
    // second pass: distance to borders
    //----------------------------------
    md = 8.0;
    for( int j=-2; j<=2; j++ )
    for( int i=-2; i<=2; i++ )
    {
        vec2 g = mg + vec2(float(i),float(j));
		vec2 o = hash22( n + g );
        vec2 r = g + o - f;

        if( dot(mr-r,mr-r)>0.00001 )
        md = min( md, dot( 0.5*(mr+r), normalize(r-mr) ) );
    }

    return vec3( md, mr );
}

// http://www.iquilezles.org/www/articles/distfunctions/distfunctions.htm
float sdEllipsoid( in vec3 p, in vec3 r ) {
    float k0 = length(p/r);
    float k1 = length(p/(r*r));
    return k0*(k0-1.0)/k1;
}


float triangle(float x) {
	return abs(1.0 - mod(abs(x), 2.0)) * 2.0 - 1.0;
}

// gamma correction
const float GAMMA = 2.2;
const float iGAMMA = 1.0 / GAMMA;
float toLinear(float c) { return pow(c,GAMMA); }
vec2 toLinear(vec2 c) { return pow(c,vec2(GAMMA)); }
vec3 toLinear(vec3 c) { return pow(c,vec3(GAMMA)); }
float toSRGB(float c) { return pow(c,iGAMMA); }
vec2 toSRGB(vec2 c) { return pow(c,vec2(iGAMMA)); }
vec3 toSRGB(vec3 c) { return pow(c,vec3(iGAMMA)); }
