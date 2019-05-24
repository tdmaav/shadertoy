/*
	"Hexagon cell uv" by Alexander Alekseev aka TDM - 2018
	License Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.
*/

// ret.x  - distance to border, 
// ret.y  - distance to center
// ret.zw - cell uv
// id - cell coordinates
vec4 hex(vec2 uv, out vec2 id) {
    uv *= mat2(1.1547,0.0,-0.5773503,1.0);
    vec2 f = fract(uv);
    float triid = 1.0;
	if((f.x+f.y) > 1.0) { f = 1.0 - f; triid = -1.0; }
    
    vec2 co = step(f.yx,f) * step(1.0-f.x-f.y,max(f.x,f.y));
    id = floor(uv) + (triid < 0.0 ? 1.0 - co : co);
    co = (f - co) * triid * mat2(0.866026,0.0,0.5,1.0);    
    
    uv = abs(co);    
    return vec4(0.5-max(uv.y,abs(dot(vec2(0.866026,0.5),uv))),length(co),co);
}


// misc stuff
const float PI	 	= 3.141592;
float diffuse(vec3 n,vec3 l) {
    return max(dot(n,l), 0.0);
}
float specular(vec3 n,vec3 l,vec3 e,float s) {    
    float nrm = (s + 8.0) / (PI * 8.0);
    return pow(max(dot(reflect(e,n),l),0.0),s) * nrm;
}
float hash( float p ) {
	float h = p * 127.1;	
    return fract(sin(h)*43758.5453123);
}
float hash( vec2 p ) {
	float h = dot(p,vec2(127.1,311.7));	
    return fract(sin(h)*43758.5453123);
}
float vnoise(in float p) {
    float i = floor( p );
    float f = fract( p );	
	float u = f*f*(3.0-2.0*f);
    return mix( hash( i ), hash( i + 1.0 ), u);
}
float vnoise(in vec2 p) {
    vec2 i = floor( p );
    vec2 f = fract( p );	
	vec2 u = f*f*(3.0-2.0*f);
    return mix( mix( hash( i + vec2(0.0,0.0) ), 
                     hash( i + vec2(1.0,0.0) ), u.x),
                mix( hash( i + vec2(0.0,1.0) ), 
                     hash( i + vec2(1.0,1.0) ), u.x), u.y);
}

void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    vec2 uv = fragCoord/iResolution.xy;
    uv.x -= 0.5;
    uv.x *= iResolution.x / iResolution.y;
    uv.x += 0.5;

    // get hexagon info
    vec2 id0, id1;
    vec4 h = hex(uv*4.0, id0);
    vec4 hl = hex(uv*2.0, id1);
    
    // color
    vec3 n = vec3(h.zw * 1.8, 1.0);
    n.z = sqrt(abs(1.0 - dot(n.xy,n.xy)));
    vec3 l = normalize(vec3(sin(iTime*0.3),cos(iTime*0.3),1.0));
    
    // roughness modifier using polar coords
    float cell_rnd = hash(id0) * 4.0 - 2.0;
    float rough = smoothstep(-0.2,0.2, sin(iTime*cell_rnd+atan(h.w,h.z)*3.0+h.y*20.0));
    float gm = rough * smoothstep(-0.4,0.4, sin(iTime*cell_rnd+atan(h.w,h.z)*3.0+h.y*20.0+2.5));
    rough = 0.3 + 0.7 * rough;
    
    vec3 col = vec3(0.01,0.01,0.02) * (1.0 - rough);
    col += diffuse(n,l) * 0.04;
    col += specular(n,l,normalize(vec3(uv,-1.0)),3.0*rough) * 0.2;
    col += specular(n,l,normalize(vec3(uv,-1.0)),20.0*rough) * 0.05;    
    
    // glow
    float ga = min(pow(min(h.x*8.0,1.0), 1.0/8.0), 1.0-gm*0.1);
    ga = 1.0 - (1.0 - ga) * pow(vnoise(uv*2.0+vec2(iTime,0.0))*1.1, 4.0) *
        smoothstep(0.1,0.15,hl.x);    
    col = mix(vec3(1.4,0.0,0.3), col, ga);
    col += 5.0 * pow(1.0-ga, 2.0);
    
    // ao
    col *= 1.0 - (1.0 - pow(h.x*8.0, 1.0/3.0)) * ga;

    // Output to screen
    fragColor = vec4(pow(col,vec3(1.0/2.2)),1.0);
}
