#define PI 3.1415926535897932384626433832795

vec3 obj_pos = vec3(0.0,0.0,-10.0);
float obj_size = 5.0;

float sphere(vec3 dir, vec3 center, float radius) {
    vec3 rp = -center;
    float b = dot(rp,dir);
    float dist = b * b - (dot(rp,rp) - radius * radius);
    if(dist <= 0.0) return -1.0;
    return -b - sqrt(dist);
}

float somestep(float t) {
    return pow(t,4.0);
}

vec3 getFishEye(vec2 uv, float level) {
    float len = length(uv);
    float a = len * level;
    return vec3(uv / len * sin(a), -cos(a));
}

vec3 textureAVG(samplerCube tex, vec3 tc) {
    const float diff0 = 0.35;
    const float diff1 = 0.12;
     vec3 s0 = texture(tex,tc).xyz;
    vec3 s1 = texture(tex,tc+vec3(diff0)).xyz;
    vec3 s2 = texture(tex,tc+vec3(-diff0)).xyz;
    vec3 s3 = texture(tex,tc+vec3(-diff0,diff0,-diff0)).xyz;
    vec3 s4 = texture(tex,tc+vec3(diff0,-diff0,diff0)).xyz;
    
    vec3 s5 = texture(tex,tc+vec3(diff1)).xyz;
    vec3 s6 = texture(tex,tc+vec3(-diff1)).xyz;
    vec3 s7 = texture(tex,tc+vec3(-diff1,diff1,-diff1)).xyz;
    vec3 s8 = texture(tex,tc+vec3(diff1,-diff1,diff1)).xyz;
    
    return (s0 + s1 + s2 + s3 + s4 + s5 + s6 + s7 + s8) * 0.111111111;
}

vec3 textureBlured(samplerCube tex, vec3 tc) {
       vec3 r = textureAVG(tex,vec3(1.0,0.0,0.0));
    vec3 t = textureAVG(tex,vec3(0.0,1.0,0.0));
    vec3 f = textureAVG(tex,vec3(0.0,0.0,1.0));
    vec3 l = textureAVG(tex,vec3(-1.0,0.0,0.0));
    vec3 b = textureAVG(tex,vec3(0.0,-1.0,0.0));
    vec3 a = textureAVG(tex,vec3(0.0,0.0,-1.0));
        
    float kr = dot(tc,vec3(1.0,0.0,0.0)) * 0.5 + 0.5; 
    float kt = dot(tc,vec3(0.0,1.0,0.0)) * 0.5 + 0.5;
    float kf = dot(tc,vec3(0.0,0.0,1.0)) * 0.5 + 0.5;
    float kl = 1.0 - kr;
    float kb = 1.0 - kt;
    float ka = 1.0 - kf;
    
    kr = somestep(kr);
    kt = somestep(kt);
    kf = somestep(kf);
    kl = somestep(kl);
    kb = somestep(kb);
    ka = somestep(ka);    
    
    float d;
    vec3 ret;
    ret  = f * kf; d  = kf;
    ret += a * ka; d += ka;
    ret += l * kl; d += kl;
    ret += r * kr; d += kr;
    ret += t * kt; d += kt;
    ret += b * kb; d += kb;
    
    return ret / d;
}

float phong(vec3 l, vec3 e, vec3 n, float power) {
    float nrm = (power + 8.0) / (PI * 8.0);
    return pow(max(dot(l,reflect(e,n)),0.0), power) * nrm;
}

// GGX code from https://www.shadertoy.com/view/MlB3DV
float G1V ( float dotNV, float k ) {
    return 1.0 / (dotNV*(1.0 - k) + k);
}
float GGX(vec3 N, vec3 V, vec3 L, float roughness, float F0) {
        float alpha = roughness*roughness;
    vec3 H = normalize (V + L);

    float dotNL = clamp (dot (N, L), 0.0, 1.0);
    float dotNV = clamp (dot (N, V), 0.0, 1.0);
    float dotNH = clamp (dot (N, H), 0.0, 1.0);
    float dotLH = clamp (dot (L, H), 0.0, 1.0);

    float D, vis;
    float F;

    // NDF : GGX
    float alphaSqr = alpha*alpha;
    float pi = 3.1415926535;
    float denom = dotNH * dotNH *(alphaSqr - 1.0) + 1.0;
    D = alphaSqr / (pi * denom * denom);

    // Fresnel (Schlick)
    float dotLH5 = pow (1.0 - dotLH, 5.0);
    F = F0 + (1.0 - F0)*(dotLH5);

    // Visibility term (G) : Smith with Schlick's approximation
    float k = alpha / 2.0;
    vis = G1V (dotNL, k) * G1V (dotNV, k);

    return /*dotNL */ D * F * vis;
}

vec3 getColor(vec3 ray) {
    float dist = sphere(ray,obj_pos,obj_size);    
    if(dist > 0.0) {
        
        vec3 point = ray * dist;
        vec3 normal = point - obj_pos;
        normal = normalize(normal);
        
        // material
        float metallic = 0.04;
        float roughness = step(fract(normal.x * 2.02), 0.5) + 0.1;
        float fresnel_pow = mix(5.0, 3.5,metallic);
        //const vec3 color_mod = vec3(1.000, 0.766, 0.336);
        vec3 color_mod = vec3(1.0);
        vec3 light_color = pow(texture(iChannel0,vec3(1.0,0.0,0.0)).xyz * 1.2, vec3(2.2));
                
                
        // IBL
        vec3 ibl_diffuse = pow(textureBlured(iChannel0,normal), vec3(2.2));
        vec3 ibl_reflection = pow(textureBlured(iChannel0,reflect(ray,normal)), vec3(2.2));
        
        // fresnel
        float fresnel = max(1.0 - dot(normal,-ray), 0.0);
        fresnel = pow(fresnel,fresnel_pow);    
        
        // reflection        
        vec3 refl = pow(texture(iChannel0,reflect(ray,normal)).xyz, vec3(2.2));
        refl = mix(refl,ibl_reflection,(1.0-fresnel)*roughness);
        refl = mix(refl,ibl_reflection,roughness);
        
        // specular
        vec3 light = normalize(vec3(-0.5,1.0,0.0));
        float power = 1.0 / max(roughness * 0.4,0.01);
        //vec3 spec = light_color * phong(light,ray,normal,power);
        vec3 spec = light_color * GGX(normal,-ray,light,roughness*0.7, 0.2);
        refl -= spec;
        
        // diffuse
        vec3 diff = ibl_diffuse * vec3(1.0,0.2,0.0);
        diff = mix(diff * color_mod,refl,fresnel);        

        vec3 color = mix(diff,refl * color_mod,metallic) + spec;
        return pow(color, vec3(1.0/2.2));
        
    } else {      
        
        return texture(iChannel0,ray).xyz;
    }
}

void mainImage( out vec4 fragColor, in vec2 fragCoord ) {   
    vec2 uv = fragCoord.xy / iResolution.xy;
    uv = uv * 2.0 - 1.0;
    uv.x *= iResolution.x / iResolution.y;
    vec3 dir = getFishEye(uv,1.2);
    
    // rotation
    float c = cos(iTime);
    float s = sin(iTime);
    dir.xz = vec2(dir.x * c - dir.z * s, dir.x * s + dir.z * c);
    obj_pos.xz = vec2(obj_pos.x * c - obj_pos.z * s, obj_pos.x * s + obj_pos.z * c);
    
    // color
    float fish_eye = smoothstep(2.0,1.6,length(uv)) * 0.25 + 0.75;
    fragColor = vec4(getColor(dir) * fish_eye, 1.0);
}

