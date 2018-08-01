/*
"Raytracing primitives" by Alexander Alekseev aka TDM - 2014
License Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.
Contact: tdmaav@gmail.com
*/

//#define SURFACE_FUNCTION() intersectionRayBox(ori,dir,vec3(1.0),r0,normal)
//#define SURFACE_FUNCTION() intersectionRayBoxWithHole(ori,dir,vec3(1.0),0.8,r0,normal)
//#define SURFACE_FUNCTION() intersectionRaySphere(ori,dir,1.0,r0,normal)
//#define SURFACE_FUNCTION() intersectionRayCylinder(ori,dir,vec2(1.,0.1),r0,normal)
#define SURFACE_FUNCTION() intersectionRayCylinderWithHole(ori,dir,vec3(1.0,0.5,0.5),r0,normal)


// math
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

/****************************************
 * intersections
 ****************************************/

bool intersectionRayQuad(vec3 o, vec3 d, vec2 size, out vec3 p) {
    p = o - d * o.z / d.z;
    return bool(mul(step(abs(p.xy),size)));
}

float intersectionRayBox(vec3 o, vec3 d, vec3 ext, out vec3 r0, out vec3 nrm) {
    vec3 t0 = (-o - ext) / d; 
    vec3 t1 = (-o + ext) / d;    
    vec3 n = min(t0,t1); n.x = max(max(n.x,n.y),n.z);
    vec3 f = max(t0,t1); f.x = min(min(f.x,f.y),f.z);
    r0 = o + d * n.x;    
    nrm = abs(r0/ext);
    nrm = step(vec3(max(max(nrm.x,nrm.y),nrm.z)),nrm) * sign(r0);     
    return step(n.x,f.x);
}

float intersectionRayBoxWithHole(vec3 o, vec3 d, vec3 ext, float hole, out vec3 r0, out vec3 nrm) {
    vec3 t0 = (-o - ext) / d; 
    vec3 t1 = (-o + ext) / d;
    vec3 n = min(t0,t1); n.x = max(max(n.x,n.y),n.z);
    vec3 f = max(t0,t1); f.x = min(min(f.x,f.y),f.z);
    r0 = o + d * n.x;
    nrm = abs(r0/ext);
    nrm = step(vec3(max(max(nrm.x,nrm.y),nrm.z)),nrm) * sign(r0);
    float ret0 = step(n.x,f.x);
    
    // inner
    vec3 ext2 = vec3(hole,ext.y,hole);
    t0 = (-o - ext2) / d; 
    t1 = (-o + ext2) / d;
    n = min(t0,t1); n.x = max(max(n.x,n.y),n.z);
    f = max(t0,t1); f.x = min(min(f.x,f.y),f.z);
    float ret1 = step(n.x,f.x);
    
    // caps
    float EPS = ext2.y - 1e-4;
    vec3 r1n = o + d * n.x;
    vec3 r1f = o + d * f.x;
    float nocap0 = step(EPS,abs(r1n.y));
    float nocap1 = step(abs(r1f.y),EPS);
    float nocap2 = step(EPS,abs(r1f.y));    
    ret0 *= 1.0 - ret1 * nocap0 * nocap2;
    ret1 *= nocap0 * nocap1;
    
    if(ret1 > 0.5) {
        r0 = r1f;
        nrm = abs(r0/ext2);
        nrm = step(vec3(max(max(nrm.x,nrm.y),nrm.z)),nrm) * -sign(r0);
    }

    return ret0;
}

float intersectionRaySphere(vec3 o, vec3 d, float radius, out vec3 r0, out vec3 nrm) {
    vec3 rp = o;
    float b = dot(rp,d);
    float dist = b * b - (dot(rp,rp) - radius * radius);
    if(dist <= 0.0) return -1.0;
    float t = -b - sqrt(dist);
    r0 = o + d * t;
    nrm = normalize(r0);
    return 1.0;
}

float intersectionRayCylinder(vec3 o, vec3 d, vec2 rh, out vec3 r0, out vec3 nrm) {
    float a = dot(d.xz,d.xz);
    float b = 2.0 * dot(o.xz,d.xz);
    float c = dot(o.xz,o.xz) - rh.x * rh.x;
    float disc = b * b - 4.0 * a * c;
    if(disc < 0.) return -1.0;
    a *= 2.0;
    float t = (-b-sqrt(disc)) / a;
    r0 = o + d * t;
    
    // caps    
    if(abs(r0.y) > rh.y) {
        nrm = vec3(0.,sign(r0.y),0.);        
        float t1 = (-b + sqrt(disc)) / a;
        vec3 r1 = o + d * t1;
        if(r1.y * nrm.y > rh.y) return -1.0;        
        r0 = o + d * ((rh.y*nrm.y-o.y) / d.y); 
        
    // side
    } else {    
        nrm.xz = normalize(r0.xz);
        nrm.y = 0.0;
    }   
    return 1.0;
}

float intersectionRayCylinderWithHole(vec3 o, vec3 d, vec3 rrh, out vec3 r0, out vec3 nrm) {
    float a = dot(d.xz,d.xz);
    float b = 2.0 * dot(o.xz,d.xz);
    float c = dot(o.xz,o.xz) - rrh.x * rrh.x;
    float disc = b * b - 4.0 * a * c;
    if(disc < 0.) return -1.0;
    a *= 2.0;
    float t = (-b-sqrt(disc)) / a;
    r0 = o + d * t;
    
    // caps    
    if(abs(r0.y) > rrh.z) {
        nrm = vec3(0.,sign(r0.y),0.);        
        float t1 = (-b + sqrt(disc)) / a;
        vec3 r1 = o + d * t1;
        if(r1.y * nrm.y > rrh.z) return -1.0;        
        r0 = o + d * ((rrh.z*nrm.y-o.y) / d.y);
        
        // inner side
        if(dot(r0.xz,r0.xz) < rrh.y * rrh.y) {
            a = dot(d.xz,d.xz);
            c = dot(o.xz,o.xz) - rrh.y * rrh.y;
            disc = b * b - 4.0 * a * c;
            t = (-b+sqrt(disc)) / (a * 2.0);
            r0 = o + d * t;
            if(abs(r0.y) > rrh.z) return -1.0;
                
            nrm.xz = -normalize(r0.xz);
            nrm.y = 0.0;            
            return 1.0;
        }
        
    // side
    } else {    
        nrm.xz = normalize(r0.xz);
        nrm.y = 0.0;
    }   
    return 1.0;
}

/****************************************
 * main
 ****************************************/

void mainImage( out vec4 fragColor, in vec2 fragCoord ) {
    vec2 iuv = fragCoord.xy / iResolution.xy * 2.0 - 1.0;
    vec2 uv = iuv;
    uv.x *= iResolution.x / iResolution.y;    
    vec2 mouse = iMouse.xy / iResolution.xy * 4.0 - 2.0;
        
    // ray
    vec3 ang = vec3(0.0,sin(iTime)*0.75,cos(iTime*1.5)*0.75);
    if(iMouse.z > 0.0) ang = vec3(0.0,-mouse.y,mouse.x);
    mat3 rot = fromEuler(ang);
    
    vec3 ori = vec3(0.0,0.0,5.0);
    vec3 dir = normalize(vec3(uv.xy,-3.0));
    ori = ori * rot;
    dir = dir * rot;
             
    // color
    vec3 p, r0, normal;
    vec3 color = vec3(0.0);
    if(SURFACE_FUNCTION() > 0.) {        
        vec3 refl = reflect(-dir,normal);
        
        float factor = pow(dot(-dir,normal),2.0);
        //factor *= texture(iChannel1,r0.xy).x * 2.0;
        
        color = mix(texture(iChannel0,refl).xyz,
                    texture(iChannel0,r0).xxx,
                    min(factor,1.0));
        
        //color = normal * 0.5 + 0.5;
    } else {
        color = vec3(0.0);
    }
               
    fragColor = vec4(color,1.0);
}
