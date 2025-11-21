#pragma once

#include <gemmi/mmread.hpp>
#include <gemmi/structure.hpp>

#include <godot_cpp/classes/mesh_instance3d.hpp>
#include <godot_cpp/classes/multi_mesh.hpp>
#include <godot_cpp/classes/multi_mesh_instance3d.hpp>
#include <godot_cpp/classes/node3d.hpp>
#include <godot_cpp/classes/surface_tool.hpp>
#include <godot_cpp/variant/color.hpp>
#include <godot_cpp/variant/string.hpp>
#include <godot_cpp/variant/vector3.hpp>

#include <string>
#include <vector>

namespace godot {

class ProteinStructure : public Node3D {
	GDCLASS(ProteinStructure, Node3D)

public:
	enum RepresentationMode {
		REPRESENTATION_ATOMS = 0,
		REPRESENTATION_STICKS = 1,
		REPRESENTATION_CARTOON = 2,
	};

	ProteinStructure() = default;
	~ProteinStructure() override = default;

	void load_file(const String &path);
	void set_representation(int mode);
	void apply_rotamer(int chain_idx, int residue_idx, float chi_angle);

protected:
	static void _bind_methods();

private:
	static constexpr float SCALE_FACTOR = 0.1f;
	static constexpr float DEFAULT_STICK_RADIUS = 0.07f;
	static constexpr int CARTOON_SUBDIVISIONS = 5;

	gemmi::Structure model_data;
	MultiMeshInstance3D *atom_mesh_instance = nullptr;
	MultiMeshInstance3D *bond_mesh_instance = nullptr;
	MeshInstance3D *ribbon_mesh_instance = nullptr;
	RepresentationMode current_mode = REPRESENTATION_ATOMS;

	void clear_visuals();
	bool ensure_structure_loaded() const;
	Vector3 to_vec3(const gemmi::Position &pos) const;
	Vector3 compute_center_of_mass() const;
	float get_element_mass(const std::string &symbol) const;
	float get_vdw_radius(const std::string &symbol) const;
	Color get_cpk_color(const std::string &symbol) const;

	void make_atom_multimesh();
	void make_bond_multimesh();
	void make_cartoon_mesh();

	void show_current_mode();
};

} // namespace godot
