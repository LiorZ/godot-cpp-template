#pragma once

#include <gemmi/elem.hpp>
#include <gemmi/metadata.hpp>
#include <gemmi/mmread.hpp>
#include <gemmi/model.hpp>

#include <godot_cpp/classes/mesh_instance3d.hpp>
#include <godot_cpp/classes/multi_mesh.hpp>
#include <godot_cpp/classes/multi_mesh_instance3d.hpp>
#include <godot_cpp/classes/node3d.hpp>
#include <godot_cpp/classes/surface_tool.hpp>
#include <godot_cpp/classes/standard_material3d.hpp>
#include <godot_cpp/core/type_info.hpp>
#include <godot_cpp/variant/color.hpp>
#include <godot_cpp/variant/string.hpp>
#include <godot_cpp/variant/vector3.hpp>

#include <string>
#include <unordered_map>
#include <vector>

namespace godot {

class ProteinStructure : public Node3D {
	GDCLASS(ProteinStructure, Node3D)

public:
	enum RepresentationMode {
		REPRESENTATION_ATOMS = 0,
		REPRESENTATION_STICKS = 1,
		REPRESENTATION_CARTOON = 2,
		REPRESENTATION_BALL_AND_STICK = 3,
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
	static constexpr float DEFAULT_STICK_RADIUS = 0.035f;
	static constexpr float BALL_STICK_RADIUS_SCALE = 0.35f;
	static constexpr float BOND_TOLERANCE = 0.45f;
	static constexpr int CARTOON_SUBDIVISIONS = 200;

	enum class SecondaryType {
		COIL = 0,
		HELIX = 1,
		SHEET = 2,
	};

	std::unordered_map<std::string, SecondaryType> secondary_lookup;
	gemmi::Structure model_data;
	MultiMeshInstance3D *atom_mesh_instance = nullptr;
	MultiMeshInstance3D *bond_mesh_instance = nullptr;
	MultiMeshInstance3D *ball_mesh_instance = nullptr;
	MeshInstance3D *ribbon_mesh_instance = nullptr;
	std::vector<std::pair<Vector3, Vector3>> bond_segments;
	RepresentationMode current_mode = REPRESENTATION_ATOMS;

	void clear_visuals();
	bool ensure_structure_loaded() const;
	Vector3 to_vec3(const gemmi::Position &pos) const;
	Vector3 compute_center_of_mass() const;
	float get_element_mass(const std::string &symbol) const;
	float get_vdw_radius(const std::string &symbol) const;
	float get_covalent_radius(const gemmi::Element &el) const;
	Color get_cpk_color(const std::string &symbol) const;
	Color get_secondary_color(SecondaryType type) const;
	float get_secondary_width(SecondaryType type) const;
	float get_secondary_thickness(SecondaryType type) const;
	std::string residue_key(const std::string &chain_name, const gemmi::SeqId &seqid) const;
	void build_secondary_lookup();
	SecondaryType get_secondary_for_residue(const std::string &chain_name, const gemmi::SeqId &seqid) const;
	Vector3 compute_residue_normal(const gemmi::Residue &residue, const Vector3 &tangent_hint) const;
	void compute_bonds();
	Ref<StandardMaterial3D> make_standard_material(bool double_sided, float roughness = 0.35f, float metallic = 0.05f) const;

	void make_atom_multimesh();
	void make_ballstick_multimesh();
	void make_bond_multimesh();
	void make_cartoon_mesh();

	void show_current_mode();
};

} // namespace godot

VARIANT_ENUM_CAST(godot::ProteinStructure::RepresentationMode);
