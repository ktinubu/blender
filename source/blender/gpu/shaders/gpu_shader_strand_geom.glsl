#define MAX_CURVE_VERTS 64

layout(points) in;
layout(line_strip, max_vertices = MAX_CURVE_VERTS) out;

in mat3 vRotation[];
in uvec4 v_control_index[];
in vec4 v_control_weight[];
in vec2 v_root_distance[];
in vec3 vColor[];

out vec3 fPosition;
out vec3 fTangent;
out vec3 fColor;

struct Samplers {
	usamplerBuffer control_curves;
	samplerBuffer control_points;
	samplerBuffer control_normals;
	samplerBuffer control_tangents;
};

uniform Samplers samplers;

bool is_valid_index(uint index)
{
	return index < uint(0xFFFFFFFF);
}

void emit_vertex(in vec3 location, in vec3 tangent)
{
	gl_Position = gl_ProjectionMatrix * gl_ModelViewMatrix * vec4(location, 1.0);
	fPosition = (gl_ModelViewMatrix * vec4(location, 1.0)).xyz;
	fTangent = gl_NormalMatrix * tangent;
	fColor = vColor[0];
	EmitVertex();
}

void main()
{
	vec3 root = gl_in[0].gl_Position.xyz;
	
	int index[4] = int[4](int(v_control_index[0].x),
	                      int(v_control_index[0].y),
	                      int(v_control_index[0].z),
	                      int(v_control_index[0].w));
	bool valid[4] = bool[4](is_valid_index(v_control_index[0].x),
		                    is_valid_index(v_control_index[0].y),
		                    is_valid_index(v_control_index[0].z),
		                    is_valid_index(v_control_index[0].w));
	float weight[4] = float[4](v_control_weight[0].x,
	                           v_control_weight[0].y,
	                           v_control_weight[0].z,
	                           v_control_weight[0].w);

	int cvert_begin[4];
	int num_cverts[4];
	float fnum_verts = 0.0;
	vec3 croot[4];
	vec3 offset[4];
	for (int k = 0; k < 4; ++k) {
		if (!valid[k])
			continue;

		uvec4 curve = texelFetch(samplers.control_curves, index[k]);
		cvert_begin[k] = int(curve.x);
		num_cverts[k] = int(curve.y);
		
		fnum_verts += weight[k] * float(num_cverts[k]);
		
		croot[k] = texelFetch(samplers.control_points, cvert_begin[k]).xyz;
		offset[k] = root - croot[k];
	}
	int num_verts = max(int(ceil(fnum_verts)), 2);

	float dt[4];
	for (int k = 0; k < 4; ++k) {
		dt[k] = float(num_cverts[k] - 1) / float(num_verts - 1);
	}

	float t[4] = float[4](0.0, 0.0, 0.0, 0.0);
	for (int i = 0; i < num_verts; ++i) {
		vec3 loc = vec3(0.0, 0.0, 0.0);
		vec3 nor = vec3(0.0, 0.0, 0.0);

		vec3 cloc[4], cnor[4], ctang[4];
		for (int k = 0; k < 4; ++k) {
			if (!valid[k])
				continue;

			interpolate_control_curve(samplers.control_points, samplers.control_normals, samplers.control_tangents,
			                          t[0], cvert_begin[0], num_cverts[0],
			                          cloc[k], cnor[k], ctang[k]);

			loc += weight[k] * (cloc[k] + offset[k]);
			nor += weight[k] * cnor[k];
			t[k] += dt[k];
		}

		mat3 cframe0 = mat3_from_vectors(cnor[0], ctang[0]);
		displace_vertex(loc, nor, float(i) / float(num_verts - 1), 1.0, cloc[0], cframe0);

		emit_vertex(loc, nor);
	}
	
	EndPrimitive();
}