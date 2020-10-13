/*
 * "Crystal Field" by Alexander Alekseev aka TDM - 2020
 * License Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.
 * Contact: tdmaav@gmail.com
 */


mat4 fromEuler(vec3 ang) {
	vec2 a1 = vec2(sin(ang.x),cos(ang.x));
    vec2 a2 = vec2(sin(ang.y),cos(ang.y));
    vec2 a3 = vec2(sin(ang.z),cos(ang.z));
    mat4 m;
    m[0] = vec4(a1.y*a3.y+a1.x*a2.x*a3.x,a1.y*a2.x*a3.x+a3.y*a1.x,-a2.y*a3.x,0.0);
	m[1] = vec4(-a2.y*a1.x,a1.y*a2.y,a2.x,0.0);
	m[2] = vec4(a3.y*a1.x*a2.x+a1.y*a3.x,a1.x*a3.x-a1.y*a3.y*a2.x,a2.y*a3.y,0.0);
	m[3] = vec4(0.0,0.0,0.0,1.0);
	return m;
}
vec2 rotateZ(vec2 v, float a) {
    vec2 sc = vec2(sin(a),cos(a));
    vec2 ret = v;
    ret.x = v.x * sc.y - v.y * sc.x;
    ret.y = v.x * sc.x + v.y * sc.y;
    return ret;
}
vec3 rotate(vec3 v, mat4 m) {
    return vec3(dot(v,m[0].xyz),dot(v,m[1].xyz),dot(v,m[2].xyz));
}
float smix(float a, float b, float t) {
    t = clamp(t,0.0,1.0);
    return mix(a,b, t*t*(3.0-2.0*t));
}
vec4 smix(float a, float b, vec4 t) {
    t = clamp(t,0.0,1.0);
    return mix(vec4(a),vec4(b), t*t*(3.0-2.0*t));
}
float tri(float x) {
    return 1.0 - abs(fract(x) - 0.5) * 4.0;
}

float intersectionRayBox(vec3 o, vec3 d, vec3 ext, out vec3 r0, out vec3 nrm, out vec2 t) {
    vec3 t0 = (-o - ext) / d;
    vec3 t1 = (-o + ext) / d;
    vec3 n = min(t0,t1); t.x = max(max(n.x,n.y),n.z);
    vec3 f = max(t0,t1); t.y = min(min(f.x,f.y),f.z);
    r0 = o + d * t.x;
    nrm = abs(r0/ext);
    nrm = step(vec3(max(max(nrm.x,nrm.y),nrm.z)),nrm) * sign(r0);
    return step(t.x,t.y);
}

// SDF
float rbox(vec3 p,vec3 s) {
	p = abs(p)-s;
    return length(p-min(p,0.0));
}
float octahedron(vec3 p, float s) {
    p = abs(p);
    float m = p.x + p.y + p.z - s;
	return m*0.57735027;
}

// hash
float hash12(in vec2 p) {
	float h = dot(p,vec2(127.1,311.7));	
    return fract(sin(h)*43758.5453123);
}
vec3 hash32(in vec2 p) {
	vec3 p3 = fract(vec3(p.xyx) * vec3(.1031, .1030, .0973));
    p3 += dot(p3, p3.yxz+33.33);
    return fract((p3.xxy+p3.yzz)*p3.zyx);
}

// comment some of these to get better performance
#define AA
#define DISPERSION

const int NUM_STEPS = 24;
const int NUM_REFR_ITER = 30;
const float TRESHOLD 	= 1e-4;
const float EPSILON 	= 1e-4;
const float PI 			= 3.141592;
const float IOR = 0.7;

float isCrystalf(vec2 c) {
    return step(0.8,hash12(c));
}
bool isCrystal(vec2 c) {
    return hash12(c) > 0.8;
}
float cellHeight(vec2 c) {
    vec3 rnd = hash32(-c*16.7);
    float h = rnd.x * 1.8;
    h += sin((c.x+c.y)*0.2+iTime*1. + rnd.y) * 0.5;
    h += sin((c.x-c.y)*0.4+iTime*.5 + rnd.z) * 0.5;
    h += isCrystalf(c);
    h += pow(dot(c,c) * 0.003,3.0);
    return h;
}

float crystalSDF(vec3 p, in vec2 cell) {
    float h = cellHeight(cell);
    p.y -= h - 100.0;
    
    //return rbox(-p,vec2(0.4,100.0).xyx);
    
    p.xz = rotateZ(p.xz,sin(iTime)+(h-2.0));
    vec3 rp = vec3(rotateZ(p.xz,PI*0.25),p.y).xzy;
    float d = rbox(-rp,vec2(0.45,100.0).xyx);
	d = max(d, rbox(-p,vec2(0.45,100.0).xyx));
    d = max(d, octahedron(p,100.0));
    return d;
}

vec3 getCrystalNormal(vec3 p, in vec2 cell, float dens) {
    vec3 n;
    n.x = crystalSDF(vec3(p.x+EPSILON,p.y,p.z), cell);
    n.y = crystalSDF(vec3(p.x,p.y+EPSILON,p.z), cell);
    n.z = crystalSDF(vec3(p.x,p.y,p.z+EPSILON), cell);
    n = n - crystalSDF(p, cell);   
    float cr = (tri(p.y*1.4) + tri(p.y*2.7)) * 0.4;
    n.y += cr * 1e-4;
    return normalize(n);
}

float getAO(vec3 p, vec3 n, in vec2 c) { 
    vec2 cell_coord = fract(p.xz);
    
    float dhc = cellHeight(c);
    vec4 dh = vec4(cellHeight(c-vec2(1.0,0.0)),
                   cellHeight(c+vec2(1.0,0.0)),
                   cellHeight(c-vec2(0.0,1.0)),
                   cellHeight(c+vec2(0.0,1.0))) - p.y;
    vec4 side = max(-dh,0.0) / (dhc - min(dh,0.0) - p.y) * 0.75;
    dh = max(dh, 1e-5);    
    
    vec4 aot = min(smix(0.5,1.0,vec4(cell_coord,1.0-cell_coord).xzyw/dh), vec4(1.0));
    vec4 aos = min(mix(vec4(0.5),vec4(1.0), 1.0-(1.0-side)*max(n.xxzz*vec4(-1.,1.,-1.,1.),0.0)), 1.0);
    float hao = 1.0 - (dhc - p.y) * 0.4;
    //aot = 1.0 - (1.0 - aot) * n.y;        
    
    float ao = aot.x * aot.y * aot.z * aot.w * 
        	   aos.x * aos.y * aos.z * aos.w * hao;
    return ao;//clamp(ao,0.0,1.0);
}


/*
 * grid tracing
 */

vec3 marchCrystal(vec3 ori, vec3 dir, vec2 t_range, vec3 cell_origin, vec2 cell, out vec3 p) {
    float t = t_range.x;
    float d = 0.0;
    for(int i = 0; i < NUM_STEPS; i++) {
        p = ori + dir * t;
        d = crystalSDF(p - cell_origin, cell);
        if(abs(d) <= TRESHOLD || t >= t_range.y) break;
        t += d - TRESHOLD;
    } 
    return vec3(d, t, step(t,t_range.y-EPSILON));
}

float advance(in vec3 ori, in vec3 dir, float t) {
    vec2 dir2 = normalize(dir.xz);
    float cosa = dot(dir.xz,dir2);
    
    // cell id & internal offset
    vec3 p = ori + dir * t;
    vec2 offset = fract(p.xz);

    // next t 
    vec2 t0 = offset / -dir2;
    vec2 t1 = (1.0-offset) / dir2;
    t0 = max(t0,t1);
    float nt = min(t0.x,t0.y);
    return t + nt / cosa;
}

vec3 trace(vec3 ori, vec3 dir, vec2 uv, int max_iter,  out vec3 p, out vec3 n, out vec2 c) {
    vec2 dir2 = normalize(dir.xz);
    float cosa = dot(dir.xz,dir2);
    
    float t = 0.0;
    for(int i = 0; i < 60 && i < max_iter; i++) {
        
        // cell id & internal offset
        p = ori + dir * t;
    	c = p.xz;
        vec2 offset = fract(c);
        c = floor(c);
        
        // next t 
        vec2 t0 = offset / -dir2;
        vec2 t1 = (1.0-offset) / dir2;
        t0 = max(t0,t1);
        float nt = min(t0.x,t0.y);
        nt = t + nt / cosa;
        
        // march the cell
        vec3 cell_origin = vec3((c+0.5),0.0).xzy;
        if(isCrystal(c)) {            
            vec3 st_res = marchCrystal(ori,dir,vec2(t,nt),cell_origin,c,p);
            if(st_res.z > 0.5) {
                n = getCrystalNormal(p - cell_origin, c, st_res.x); 
                return st_res;
            }
        } else {
            vec2 nf;
    		float h = cellHeight(c);
            vec3 ori2 = ori - cell_origin;
    		ori2.y -= h - 100.0;
            if(intersectionRayBox(ori2,dir,vec2(0.45,100.0).xyx,p,n,nf) > 0.5) {
        		p = ori + dir * nf.x;
                return vec3(0.0,nf.x,1.0);
            }
        }
        
        t = nt + EPSILON;
    }
    return vec3(0.0);
}

/*
 * color
 */

vec3 getCellColor(vec2 cell) {
    vec3 rnd = hash32(cell);   
    rnd = vec3(step(0.9,rnd.x), step(0.7,rnd.x), step(0.8,rnd.x));
    return rnd * 0.5 + 0.5;
}

vec3 getBlockColor(vec3 p, vec3 n, vec2 cell) {
	vec3 localy_p = p - vec3(0.0,cellHeight(cell),0.0);
    vec2 tuv = mix(mix(localy_p.xz * 1.2,
                       localy_p.yz * 1.13,
                       abs(n.x)),
                   localy_p.yx * 0.92,
                   abs(n.z));

    vec3 color = texture(iChannel0,tuv).xyy * vec3(1.0,1.3,1.3);
    color = pow(color,vec3(2.2));// * getCellColor(cell);

    vec3 l0 = normalize(vec3(1.0,0.5,0.4));
    color += color * max(dot(n,l0),0.0);
    
    color *= getAO(p,n,cell);
    
    return color;
}

vec3 addFog(vec3 color, vec3 dir, vec3 traceinfo) {
    vec3 bg = vec3(1.5) * (dir.y * 0.25 + 0.75) * smoothstep(-1.0,0.0,dir.y);
    float fog = clamp(log(traceinfo.y*0.08)/log(10.0), 0.0,1.0);
    fog = 1.0 - (1.0-fog) * traceinfo.z;
    return mix(color,bg, fog);
}

vec3 getRefractedColor(vec3 ori, vec3 dir, vec3 n, vec3 p, vec2 c, float eta) {
    vec3 rdir = refract(dir,n,eta);
    ori = p + rdir * advance(p,rdir,0.0);   

    vec3 cr_p, cr_n;
    vec2 cr_cell;
    vec3 cell_origin = vec3((c+0.5),0.0).xzy;           
    vec3 st_res = marchCrystal(ori,-rdir,vec2(0.0,1.5),cell_origin,c,cr_p);
    cr_n = getCrystalNormal(cr_p - cell_origin, c, st_res.x);                               
    
    //vec3 refl_dir = reflect(rdir,-cr_n);    
    rdir = refract(rdir,-cr_n,eta);
    //float fresnel = pow(1.0 - max(dot(cr_n,rdir),0.0), 5.0);
    ori = cr_p + rdir * advance(cr_p,rdir,0.0);

    vec3 cr_ti = trace(ori,rdir,vec2(0.0),NUM_REFR_ITER, cr_p,cr_n,cr_cell);
    vec3 color = getBlockColor(cr_p, cr_n, cr_cell);
    
    // internal reflection    
    /*ori = cr_p + rdir * advance(cr_p,rdir,0.0);
    cr_ti = trace(ori,refl_dir,vec2(0.0),NUM_REFR_ITER, cr_p,cr_n,cr_cell);
    vec3 refl_color = getBlockColor(cr_p, cr_n, cr_cell);        
    color = mix(color,refl_color,fresnel);*/
    
    color = addFog(color,rdir,cr_ti);

    return color;
}

vec3 getPixel(vec2 uv) {
    // ray
    vec3 ang = vec3(0.0,0.56 - sin(iTime*0.2)*0.05,iTime * 0.1);
    mat4 rot = fromEuler(ang);
    
    vec3 ori = vec3(0.0,-4.0,15.0);
    vec3 dir = normalize(vec3(uv.xy,-1.7)); dir.z += length(uv) * 0.08;    
    
    // DOF
#ifdef AA    
    vec2 offset = (hash32(uv*825.7).xy * 2.0 - 1.0) * 0.01;
    dir.xy += offset;
    ori.xy -= offset * 5.;    
#endif
    
    dir = normalize(dir);
    ori = rotate(ori,rot);
    dir = rotate(dir,rot);
    dir.xz = rotateZ(dir.xz,sin(iTime*0.17)*0.5);
    
    // tracing
    vec3 p,n;    
    vec2 cell;
    vec3 traceinfo = trace(ori,dir,uv,100, p,n,cell);        
    vec3 color;
    
    // crystal
    if(isCrystal(cell)) {     
    	#ifdef DISPERSION
            color.z = getRefractedColor(ori,dir,n,p,cell,mix(IOR,1.0,0.75)).z;
            color.yz += getRefractedColor(ori,dir,n,p,cell,mix(IOR,1.0,0.5)).yz;
            color.xy += getRefractedColor(ori,dir,n,p,cell,mix(IOR,1.0,0.25)).xy;
            color.x += getRefractedColor(ori,dir,n,p,cell,IOR).x;
        #else
        	color = getRefractedColor(ori,dir,n,p,cell,IOR);
        #endif
        color *= getCellColor(cell);
        
        
    // regular block
    } else {
        color = getBlockColor(p, n, cell);        
    }    
    
    // reflection 
    float fresnel = pow(1.0 - max(dot(n,-dir),0.0), 4.0);
    vec3 rdir = reflect(dir,n);
    ori = p + rdir * advance(p,rdir,0.0); 
    vec3 r_traceinfo = trace(ori,rdir,vec2(0.0),NUM_REFR_ITER, p,n,cell);
    vec3 refl = getBlockColor(p, n, cell); 
    refl = addFog(refl,rdir,r_traceinfo); 
    color = mix(color,refl, fresnel);        
        
    color = addFog(color,dir,traceinfo);
    return color;
}

/*
 * main
 */

void mainImage( out vec4 fragColor, in vec2 fragCoord ) { 
#ifdef AA
    vec3 color = vec3(0.0);
    for(int i = -1; i <= 1; i++) {
        for(int j = -1; j <= 1; j++) {
            vec2 uv = (fragCoord+vec2(i,j)/3.0) / iResolution.xy;
            uv = uv * 2.0 - 1.0;
            uv.x *= iResolution.x / iResolution.y; 
    		color += getPixel(uv);
        }
    }
    color /= 9.0;
#else    
    vec2 uv = fragCoord.xy / iResolution.xy;
    uv = uv * 2.0 - 1.0;
    uv.x *= iResolution.x / iResolution.y; 
    vec3 color = getPixel(uv);
#endif
    
    color *= vec3(1.4,1.7,1.7);
    color = pow(color,vec3(1.0/2.2));
    
	fragColor = vec4(color,1.0);
}

