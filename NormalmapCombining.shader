/*
 * |===============================================|
 * |      ORIGINAL      |   DERIVATIVES ADDITION   |
 * |===============================================|
 * |  NORMALS ADDITION  |   DERIVATIVES BLENDING   |
 * |===============================================|
 */

// normalmap texture
vec3 textureNormal(vec2 uv) {
    uv = fract(uv) * 3.0 - 1.5;    
        
    vec3 ret;
    ret.xy = sqrt(uv * uv) * sign(uv);
    ret.z = sqrt(abs(1.0 - dot(ret.xy,ret.xy)));
    ret = ret * 0.5 + 0.5;    
    return mix(vec3(0.5,0.5,1.0), ret, smoothstep(1.0,0.98,dot(uv,uv)));
}

// normals combine: normals addition
vec3 combineNormals0(vec3 n0, vec3 n1) {
    n0 = n0 * 2.0 - 1.0;
    n1 = n1 * 2.0 - 1.0;
    return normalize(n0 + n1) * 0.5 + 0.5;
}

// normals combine: derivatives addition
vec3 combineNormals1(vec3 n0, vec3 n1) {
    n0 = n0 * 2.0 - 1.0;
    n1 = n1 * 2.0 - 1.0;
    n0 = vec3(n0.xy + n1.xy, n0.z * n1.z);
    return normalize(n0) * 0.5 + 0.5;
}

// normals combine: derivatives blending
vec3 combineNormals2(vec3 n0, vec3 n1) {
    n0 = n0 * 2.0 - 1.0;
    n1 = n1 * 2.0 - 1.0;    
	n0  = vec3(n0.xy * n1.z + n1.xy * n0.z, n0.z * n1.z);    
    return normalize(n0) * 0.5 + 0.5;
}

//
void mainImage( out vec4 fragColor, in vec2 fragCoord ) {
	vec2 uv = fragCoord.xy / iResolution.xy;
    uv = (uv * 2.0 - 1.0);
    uv.x *= iResolution.x / iResolution.y;
    
    vec3 color;
    vec2 time = vec2(sin(iGlobalTime * 0.1), cos(iGlobalTime * 0.1));
    vec3 n0 = textureNormal(uv + time);
    vec3 n1 = textureNormal(uv - time + vec2(0.25));
    
    // combine normals
    if(uv.x < 0.0) {
        if(uv.y > 0.0) {
    		color = n0;
        } else {
            color = combineNormals0(n0,n1);
        }
    } else {
        if(uv.y > 0.0) {
        	color = combineNormals1(n0,n1);
        } else {            
        	color = combineNormals2(n0,n1);
        }
    }
    
    // borders
    color += max(smoothstep(0.01,0.005,abs(uv.x)),0.0);
    color += max(smoothstep(0.01,0.005,abs(uv.y)),0.0);
    
	fragColor = vec4(color,1.0);
}
