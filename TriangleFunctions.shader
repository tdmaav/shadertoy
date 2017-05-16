// "Triangle functions" by Alexander Alekseev aka TDM - 2016
// License Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.

const int NUM_STEPS = 32;
const int AO_SAMPLES = 4;
const vec2 AO_PARAM = vec2(1.2, 3.8);
const vec2 CORNER_PARAM = vec2(0.25, 40.0);
const float INV_AO_SAMPLES = 1.0 / float(AO_SAMPLES);
const float TRESHOLD     = 0.2;
const float EPSILON     = 1e-4;
const float LIGHT_INTENSITY = 0.3;
const vec3 WHITE     = vec3(1.2,1.07,0.98) * LIGHT_INTENSITY;

vec3 tri[3];

// math
float saturate(float x) { return clamp(x,0.,1.); }
float boolUnion(float a,float b) { return min(a,b); }    

/******************************************************************
 triangle distance
 ******************************************************************/

float triangleDistance(vec3 p,vec3 v0,vec3 v1,vec3 v2) {
    vec3 e0 = v1-v0, e1 = v2-v0, e2 = v2-v1;
    vec3 pe0 = p-v0, pe1 = p-v1;
    vec3 normal = cross(e0,e1);
                      
    // distance to plane  
    if(dot(normal,cross(pe0,e0)) < 0.0 &&
       dot(normal,cross(pe1,e2)) < 0.0 &&
       dot(normal,cross(pe0,e1)) > 0.0) {
        
        return abs(dot(p-v0,normalize(normal)));
        
    // distance to edges
    } else {
        vec3 dp0 = e0 * saturate(dot(e0,pe0) / dot(e0,e0)) - pe0;
        vec3 dp1 = e1 * saturate(dot(e1,pe0) / dot(e1,e1)) - pe0;
        vec3 dp2 = e2 * saturate(dot(e2,pe1) / dot(e2,e2)) - pe1;                                 
        return sqrt(min(min(dot(dp0,dp0),dot(dp1,dp1)),dot(dp2,dp2)));
    }
}

#if 1
/******************************************************************
 triangle intersection - based on Cramer's rule
     27 multiplications (6*2+3*4+3)
     1 division
 ******************************************************************/

float triangleIntersection(vec3 o,vec3 d, vec3 v0,vec3 v1,vec3 v2) {
    vec3 e0 = v1-v0, e1 = v2-v0, r = o - v0;
    vec3 ce0d = cross(e0,d);
    vec3 ce1r = cross(e1,r);
    float idet = 1.0 / dot(e1,ce0d);
    
    // check intersection
       float u = dot(d,ce1r) * idet;    
    if(u < 0.0 || u > 1.0) return -1.0;  
       float v = dot(r,ce0d) * idet;    
    if(v < 0.0 || (u+v) > 1.0) return -1.0;  
    
    // distance to triangle
    return dot(e0,ce1r) * idet;
}

#else

/******************************************************************
 triangle intersection - simple approach
     45 multiplications (6×4+3×6+3)
     1 division
 ******************************************************************/

float triangleIntersection(vec3 o,vec3 d, vec3 v0,vec3 v1,vec3 v2) {
    vec3 e0 = v1-v0, e1 = v2-v0;
    vec3 normal = cross(e0,e1);
    
    // distance to plane    
    float t = -(dot(normal,o) - dot(normal,v0)) / dot(normal,d);
    
    // point inside triangle
    vec3 point = o + d * t;
    if(dot(normal, cross(e0,point-v0)) < 0.0) return -1.0;
    if(dot(normal, cross(e1,v2-point)) < 0.0) return -1.0;
    if(dot(normal, cross(v2-v1,point-v1)) < 0.0) return -1.0;   
    
    return t;
}

#endif


// lighting
float diffuse(vec3 n,vec3 l,float p) { return pow(max(dot(n,l),0.0),p); }
float specular(vec3 n,vec3 l,vec3 e,float s) {    
    float nrm = (s + 8.0) / (3.1415 * 8.0);
    return pow(max(dot(reflect(e,n),l),0.0),s) * nrm;
}

float plane(vec3 gp, vec4 p) {
    return dot(p.xyz,gp+p.xyz*p.w);
}

// world
float map_tri(vec3 p) {        
    return triangleDistance(p, tri[0],tri[1],tri[2]);
}
float map(vec3 p) {
    float d = map_tri(p);
    d = boolUnion(d,plane(p,vec4(0.0,1.0,0.0,1.3)));
    return d;
}

// tracing
vec3 getNormal(vec3 p, float dens) {
    vec3 n;
    n.x = map(vec3(p.x+EPSILON,p.y,p.z));
    n.y = map(vec3(p.x,p.y+EPSILON,p.z));
    n.z = map(vec3(p.x,p.y,p.z+EPSILON));
    return normalize(n-map(p));
}
float getOcclusion(vec3 p, vec3 n) {
    float r = 0.0;
    for(int i = 0; i < AO_SAMPLES; i++) {
        float f = float(i)*INV_AO_SAMPLES;
        float h = 0.01+f*AO_PARAM.x;
        float d = map(p + n * h) - TRESHOLD;
        r += clamp(h-d,0.0,1.0) * (1.0-f);
    }    
    return pow(clamp(1.0-r*INV_AO_SAMPLES*AO_PARAM.y,0.0,1.0),0.25);
}
vec2 spheretracing(vec3 ori, vec3 dir, out vec3 p) {
    vec2 td = vec2(0.0);
    for(int i = 0; i < NUM_STEPS; i++) {
        p = ori + dir * td.x;
        td.y = map(p);
        if(td.y < TRESHOLD) break;
        td.x += td.y-TRESHOLD;
    }
    return td;
}

// object
vec3 getObjectColor(vec3 p, vec3 l, vec3 n, vec3 e) {
    vec3 color = vec3(0.0,0.5,0.5);
    color += vec3(diffuse(n,l,1.0) * WHITE);
    color += vec3(specular(n,l,e,20.0) * WHITE); 
    color *= (n.y * 0.2 + 0.8);
    return color;
}

// main
void mainImage( out vec4 fragColor, in vec2 fragCoord ) {
    vec2 iuv = fragCoord.xy / iResolution.xy * 2.0 - 1.0;
    vec2 uv = iuv;
    uv.x *= iResolution.x / iResolution.y;    
        
    // ray    
    vec3 ori = vec3(0.0,0.1,4.0);
    vec3 dir = normalize(vec3(uv.xy,-2.0));
    vec2 sc = vec2(cos(iGlobalTime),sin(iGlobalTime));
    dir.xz = vec2(dir.x * sc.x - dir.z * sc.y, dir.x * sc.y + dir.z * sc.x);
    ori.xz = vec2(ori.x * sc.x - ori.z * sc.y, ori.x * sc.y + ori.z * sc.x);
    
    // triangle    
    vec2 sc2 = vec2(sin(iGlobalTime),cos(iGlobalTime)) * 0.3;
    tri[0] = vec3(-1.0+sc2.x,    -0.4+sc2.y,    -sc2.x*2.0);
    tri[1] = vec3( 0.0+sc2.y,     1.2+sc2.x,    sc2.y);
    tri[2] = vec3( 1.0-sc2.x,    -0.4-sc2.y,    sc2.x*2.0);        
        
    // tracing
    vec3 p;
    vec2 td = spheretracing(ori,dir,p);
    vec3 n = getNormal(p,td.y);
    float occ = getOcclusion(p,n);
    vec3 light = normalize(-dir); 
         
    // color
    vec3 color = vec3(1.0);    
    if(map_tri(p) - EPSILON <= td.y) color = getObjectColor(p,light,n,dir);
    
    float t = triangleIntersection(ori, dir, tri[0],tri[1],tri[2]);
    if(t > 0.0) color *= 0.75;    
    
    color *= occ;    
    color = sqrt(color);
    
    fragColor = vec4(color,1.0);
}
