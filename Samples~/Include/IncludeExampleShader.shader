Shader "StijnVanCoillie/RayMarching/Helpers/IncludeExampleShader"
{
    Properties
    {
    }
    SubShader
    {
        Tags { "RenderType"="Opaque" }
        LOD 100

        Pass
        {
            CGPROGRAM
            #pragma vertex vert
            #pragma fragment frag

            #include "UnityCG.cginc"

            // How to include the helpers.
            // TODO: Automate this process...
            #include "Packages/com.stijnvancoillie.raymarching.helpers/Runtime/SDF_2D.cginc"

            struct appdata
            {
                float4 vertex : POSITION;
                float2 uv : TEXCOORD0;
            };

            struct v2f
            {
                float2 uv : TEXCOORD0;
                float4 vertex : SV_POSITION;
            };

            v2f vert (appdata v)
            {
                v2f o;
                o.vertex = UnityObjectToClipPos(v.vertex);
                o.uv = v.uv;
                return o;
            }

            fixed4 frag (v2f i) : SV_Target
            {
                return step( sdCircle(i.uv - 0.5, 0.5), 0);
            }
            ENDCG
        }
    }
}
