float2x2 rot(float a)
{
    float s = sin(a);
    float c = cos(a);
    return float2x2(c, -s, s, c);
}

float3x3 rotate_x(float a) { float sa = sin(a); float ca = cos(a); return float3x3(float3(1., .0, .0), float3(.0, ca, sa), float3(.0, -sa, ca)); }
float3x3 rotate_y(float a) { float sa = sin(a); float ca = cos(a); return float3x3(float3(ca, .0, sa), float3(.0, 1., .0), float3(-sa, .0, ca)); }
float3x3 rotate_z(float a) { float sa = sin(a); float ca = cos(a); return float3x3(float3(ca, sa, .0), float3(-sa, ca, .0), float3(.0, .0, 1.)); }