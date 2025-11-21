#include "protein_structure.h"

#include <algorithm>
#include <cmath>

#include <godot_cpp/core/math.hpp>
#include <godot_cpp/classes/array_mesh.hpp>
#include <godot_cpp/classes/cylinder_mesh.hpp>
#include <godot_cpp/classes/file_access.hpp>
#include <godot_cpp/classes/sphere_mesh.hpp>
#include <godot_cpp/variant/utility_functions.hpp>

using namespace godot;

namespace {

Vector3 catmull_rom(const Vector3 &p0, const Vector3 &p1, const Vector3 &p2, const Vector3 &p3, float t) {
	const float t2 = t * t;
	const float t3 = t2 * t;
	return 0.5f * ((2.0f * p1) + (-p0 + p2) * t + (2.0f * p0 - 5.0f * p1 + 4.0f * p2 - p3) * t2 +
				   (-p0 + 3.0f * p1 - 3.0f * p2 + p3) * t3);
}

} // namespace

void ProteinStructure::_bind_methods() {
	ClassDB::bind_method(D_METHOD("load_file", "path"), &ProteinStructure::load_file);
	ClassDB::bind_method(D_METHOD("set_representation", "mode"), &ProteinStructure::set_representation);
	ClassDB::bind_method(D_METHOD("apply_rotamer", "chain_idx", "residue_idx", "chi_angle"), &ProteinStructure::apply_rotamer);

	BIND_ENUM_CONSTANT(REPRESENTATION_ATOMS);
	BIND_ENUM_CONSTANT(REPRESENTATION_STICKS);
	BIND_ENUM_CONSTANT(REPRESENTATION_CARTOON);
}

Vector3 ProteinStructure::to_vec3(const gemmi::Position &pos) const {
	return Vector3(static_cast<float>(pos.x), static_cast<float>(pos.y), static_cast<float>(pos.z));
}

float ProteinStructure::get_element_mass(const std::string &symbol) const {
	if (symbol == "H") {
		return 1.008f;
	}
	if (symbol == "C") {
		return 12.011f;
	}
	if (symbol == "N") {
		return 14.007f;
	}
	if (symbol == "O") {
		return 15.999f;
	}
	if (symbol == "S") {
		return 32.06f;
	}
	if (symbol == "P") {
		return 30.974f;
	}
	return 1.0f;
}

float ProteinStructure::get_vdw_radius(const std::string &symbol) const {
	if (symbol == "H") {
		return 1.2f;
	}
	if (symbol == "C") {
		return 1.7f;
	}
	if (symbol == "N") {
		return 1.55f;
	}
	if (symbol == "O") {
		return 1.52f;
	}
	if (symbol == "S") {
		return 1.8f;
	}
	if (symbol == "P") {
		return 1.8f;
	}
	return 1.5f;
}

Color ProteinStructure::get_cpk_color(const std::string &symbol) const {
	if (symbol == "H") {
		return Color(1.0f, 1.0f, 1.0f);
	}
	if (symbol == "C") {
		return Color(0.2f, 0.2f, 0.2f);
	}
	if (symbol == "N") {
		return Color(0.1f, 0.1f, 0.8f);
	}
	if (symbol == "O") {
		return Color(0.8f, 0.1f, 0.1f);
	}
	if (symbol == "S") {
		return Color(0.95f, 0.85f, 0.1f);
	}
	if (symbol == "P") {
		return Color(1.0f, 0.5f, 0.0f);
	}
	return Color(0.6f, 0.6f, 0.6f);
}

Vector3 ProteinStructure::compute_center_of_mass() const {
	Vector3 sum(0, 0, 0);
	float total_mass = 0.0f;
	for (const gemmi::Model &model : model_data.models) {
		for (const gemmi::Chain &chain : model.chains) {
			for (const gemmi::Residue &residue : chain.residues) {
				for (const gemmi::Atom &atom : residue.atoms) {
					const std::string symbol(atom.element.name());
					const float mass = get_element_mass(symbol);
					sum += to_vec3(atom.pos) * mass;
					total_mass += mass;
				}
			}
		}
	}

	if (total_mass <= 0.0f) {
		return Vector3(0, 0, 0);
	}
	return sum / total_mass;
}

void ProteinStructure::clear_visuals() {
	if (atom_mesh_instance != nullptr) {
		atom_mesh_instance->queue_free();
		atom_mesh_instance = nullptr;
	}
	if (bond_mesh_instance != nullptr) {
		bond_mesh_instance->queue_free();
		bond_mesh_instance = nullptr;
	}
	if (ribbon_mesh_instance != nullptr) {
		ribbon_mesh_instance->queue_free();
		ribbon_mesh_instance = nullptr;
	}
}

bool ProteinStructure::ensure_structure_loaded() const {
	if (!model_data.models.empty()) {
		return true;
	}
	UtilityFunctions::push_warning("ProteinStructure: No structure loaded. Call load_file() first.");
	return false;
}

void ProteinStructure::load_file(const String &path) {
	if (path.is_empty()) {
		UtilityFunctions::printerr("ProteinStructure: Provided path is empty.");
		return;
	}

	if (!FileAccess::file_exists(path)) {
		UtilityFunctions::printerr(vformat("ProteinStructure: File does not exist: %s", path));
		return;
	}

	const std::string native_path = path.utf8().get_data();
	try {
		model_data = gemmi::read_structure(native_path);
	} catch (const std::exception &e) {
		UtilityFunctions::printerr(vformat("ProteinStructure: Failed to parse structure: %s", String(e.what())));
		return;
	}

	if (model_data.models.empty()) {
		UtilityFunctions::printerr("ProteinStructure: Parsed structure is empty.");
		return;
	}

	// Center structure to avoid floating point jitter in Godot.
	const Vector3 center = compute_center_of_mass();
	for (gemmi::Model &model : model_data.models) {
		for (gemmi::Chain &chain : model.chains) {
			for (gemmi::Residue &residue : chain.residues) {
				for (gemmi::Atom &atom : residue.atoms) {
					atom.pos.x -= center.x;
					atom.pos.y -= center.y;
					atom.pos.z -= center.z;
				}
			}
		}
	}

	clear_visuals();
	make_atom_multimesh();
	make_bond_multimesh();
	make_cartoon_mesh();
	show_current_mode();
}

void ProteinStructure::show_current_mode() {
	if (atom_mesh_instance != nullptr) {
		atom_mesh_instance->set_visible(current_mode == REPRESENTATION_ATOMS);
	}
	if (bond_mesh_instance != nullptr) {
		bond_mesh_instance->set_visible(current_mode == REPRESENTATION_STICKS);
	}
	if (ribbon_mesh_instance != nullptr) {
		ribbon_mesh_instance->set_visible(current_mode == REPRESENTATION_CARTOON);
	}
}

void ProteinStructure::set_representation(int mode) {
	RepresentationMode requested = static_cast<RepresentationMode>(mode);
	switch (requested) {
		case REPRESENTATION_ATOMS:
		case REPRESENTATION_STICKS:
		case REPRESENTATION_CARTOON:
			current_mode = requested;
			break;
		default:
			UtilityFunctions::push_warning("ProteinStructure: Invalid representation mode, keeping previous value.");
			break;
	}

	show_current_mode();
}

void ProteinStructure::apply_rotamer(int chain_idx, int residue_idx, float chi_angle) {
	// Placeholder for rotamer edits. Updating MultiMeshes incrementally requires gemmi topology data.
	UtilityFunctions::push_warning(vformat(
			"ProteinStructure: apply_rotamer not implemented yet (chain %d residue %d chi %f).", chain_idx, residue_idx, chi_angle));
}

void ProteinStructure::make_atom_multimesh() {
	if (!ensure_structure_loaded()) {
		return;
	}

	int atom_count = 0;
	for (const gemmi::Model &model : model_data.models) {
		for (const gemmi::Chain &chain : model.chains) {
			for (const gemmi::Residue &residue : chain.residues) {
				atom_count += static_cast<int>(residue.atoms.size());
			}
		}
	}

	if (atom_count == 0) {
		UtilityFunctions::push_warning("ProteinStructure: No atoms found for MultiMesh generation.");
		return;
	}

	Ref<MultiMesh> multimesh;
	multimesh.instantiate();
	multimesh->set_transform_format(MultiMesh::TRANSFORM_3D);
	multimesh->set_use_colors(true);
	multimesh->set_instance_count(atom_count);

	Ref<SphereMesh> sphere;
	sphere.instantiate();
	sphere->set_radius(0.5f);
	multimesh->set_mesh(sphere);

	int idx = 0;
	for (const gemmi::Model &model : model_data.models) {
		for (const gemmi::Chain &chain : model.chains) {
			for (const gemmi::Residue &residue : chain.residues) {
				for (const gemmi::Atom &atom : residue.atoms) {
					const std::string symbol(atom.element.name());
					const float radius = get_vdw_radius(symbol) * SCALE_FACTOR;
					const Color color = get_cpk_color(symbol);
					Transform3D transform;
					transform.origin = to_vec3(atom.pos) * SCALE_FACTOR;
					transform.basis.scale(Vector3(radius, radius, radius));
					multimesh->set_instance_transform(idx, transform);
					multimesh->set_instance_color(idx, color);
					++idx;
				}
			}
		}
	}

	atom_mesh_instance = memnew(MultiMeshInstance3D);
	atom_mesh_instance->set_multimesh(multimesh);
	add_child(atom_mesh_instance);
}

void ProteinStructure::make_bond_multimesh() {
	if (!ensure_structure_loaded()) {
		return;
	}

	struct AtomInfo {
		Vector3 pos;
		std::string name;
	};

	std::vector<std::pair<Vector3, Vector3>> bonds;
	for (const gemmi::Model &model : model_data.models) {
		for (const gemmi::Chain &chain : model.chains) {
			std::vector<AtomInfo> residue_atoms;
			std::vector<Vector3> backbone_atoms;
			backbone_atoms.reserve(chain.residues.size());

			for (const gemmi::Residue &residue : chain.residues) {
				residue_atoms.clear();
				Vector3 ca_pos;
				bool has_ca = false;

				for (const gemmi::Atom &atom : residue.atoms) {
					const Vector3 pos = to_vec3(atom.pos) * SCALE_FACTOR;
					residue_atoms.push_back({pos, atom.name});
					if (atom.name == "CA") {
						ca_pos = pos;
						has_ca = true;
					}
				}

				// Bonds within residue.
				for (size_t i = 0; i < residue_atoms.size(); ++i) {
					for (size_t j = i + 1; j < residue_atoms.size(); ++j) {
						const float distance = residue_atoms[i].pos.distance_to(residue_atoms[j].pos);
						if (distance <= 1.8f * SCALE_FACTOR) {
							bonds.emplace_back(residue_atoms[i].pos, residue_atoms[j].pos);
						}
					}
				}

				if (has_ca) {
					backbone_atoms.push_back(ca_pos);
				} else {
					backbone_atoms.push_back(Vector3());
				}
			}

			// Connect CA atoms between neighboring residues to highlight backbone continuity.
			for (size_t i = 0; i + 1 < backbone_atoms.size(); ++i) {
				if (backbone_atoms[i] != Vector3() && backbone_atoms[i + 1] != Vector3()) {
					bonds.emplace_back(backbone_atoms[i], backbone_atoms[i + 1]);
				}
			}
		}
	}

	if (bonds.empty()) {
		return;
	}

	Ref<MultiMesh> multimesh;
	multimesh.instantiate();
	multimesh->set_transform_format(MultiMesh::TRANSFORM_3D);
	multimesh->set_use_colors(true);
	multimesh->set_instance_count(static_cast<int>(bonds.size()));

	Ref<CylinderMesh> cylinder;
	cylinder.instantiate();
	cylinder->set_height(1.0f);
	cylinder->set_top_radius(0.5f);
	cylinder->set_bottom_radius(0.5f);
	multimesh->set_mesh(cylinder);

	for (int i = 0; i < static_cast<int>(bonds.size()); ++i) {
		const Vector3 a = bonds[i].first;
		const Vector3 b = bonds[i].second;
		const Vector3 direction = b - a;
		const float length = direction.length();
		if (length <= 0.0f) {
			continue;
		}

		const Vector3 center = a + direction * 0.5f;
		Vector3 up = direction.normalized();
		Vector3 basis_x = up.cross(Vector3(0, 1, 0));
		if (basis_x.length_squared() < CMP_EPSILON) {
			basis_x = up.cross(Vector3(1, 0, 0));
		}
		basis_x.normalize();
		const Vector3 basis_z = basis_x.cross(up).normalized();

		Basis basis(basis_x, up, basis_z);
		basis.scale(Vector3(DEFAULT_STICK_RADIUS, length * 0.5f, DEFAULT_STICK_RADIUS));
		const Transform3D transform(basis, center);

		multimesh->set_instance_transform(i, transform);
		multimesh->set_instance_color(i, Color(0.7f, 0.7f, 0.7f));
	}

	bond_mesh_instance = memnew(MultiMeshInstance3D);
	bond_mesh_instance->set_multimesh(multimesh);
	add_child(bond_mesh_instance);
	bond_mesh_instance->set_visible(false);
}

void ProteinStructure::make_cartoon_mesh() {
	if (!ensure_structure_loaded()) {
		return;
	}

	Ref<SurfaceTool> st;
	st.instantiate();
	st->begin(Mesh::PRIMITIVE_TRIANGLES);

	const float ribbon_half_width = 0.4f * SCALE_FACTOR;
	const float ribbon_thickness = 0.05f * SCALE_FACTOR;

	for (const gemmi::Model &model : model_data.models) {
		for (const gemmi::Chain &chain : model.chains) {
			std::vector<Vector3> ca_positions;
			std::vector<Vector3> up_vectors;
			for (const gemmi::Residue &residue : chain.residues) {
				Vector3 ca_pos;
				Vector3 oxygen_pos;
				bool has_ca = false;
				bool has_oxygen = false;
				for (const gemmi::Atom &atom : residue.atoms) {
					if (atom.name == "CA") {
						ca_pos = to_vec3(atom.pos) * SCALE_FACTOR;
						has_ca = true;
					} else if (atom.name == "O") {
						oxygen_pos = to_vec3(atom.pos) * SCALE_FACTOR;
						has_oxygen = true;
					}
				}

				if (!has_ca) {
					continue;
				}

				ca_positions.push_back(ca_pos);
				if (has_oxygen) {
					Vector3 up = (oxygen_pos - ca_pos);
					if (up.length_squared() < CMP_EPSILON) {
						up = Vector3(0, 1, 0);
					}
					up_vectors.push_back(up.normalized());
				} else {
					up_vectors.push_back(Vector3(0, 1, 0));
				}
			}

			if (ca_positions.size() < 2) {
				continue;
			}

			// Pad endpoints for Catmull-Rom.
			std::vector<Vector3> padded = ca_positions;
			padded.insert(padded.begin(), ca_positions.front());
			padded.push_back(ca_positions.back());

			for (int i = 0; i < static_cast<int>(ca_positions.size()) - 1; ++i) {
				const Vector3 p0 = padded[std::max(0, i)];
				const Vector3 p1 = padded[i + 1];
				const Vector3 p2 = padded[i + 2];
				const Vector3 p3 = padded[std::min(static_cast<int>(padded.size()) - 1, i + 3)];

				const Vector3 up0 = up_vectors[i];
				const Vector3 up1 = up_vectors[i + 1];

				for (int step = 0; step < CARTOON_SUBDIVISIONS; ++step) {
					const float t0 = static_cast<float>(step) / CARTOON_SUBDIVISIONS;
					const float t1 = static_cast<float>(step + 1) / CARTOON_SUBDIVISIONS;

					const Vector3 a = catmull_rom(p0, p1, p2, p3, t0);
					const Vector3 b = catmull_rom(p0, p1, p2, p3, t1);

					Vector3 up = (up0.lerp(up1, t0)).normalized();
					if (up.length_squared() < CMP_EPSILON) {
						up = Vector3(0, 1, 0);
					}

					Vector3 tangent = (b - a).normalized();
					if (tangent.length_squared() < CMP_EPSILON) {
						continue;
					}

					Vector3 side = tangent.cross(up).normalized() * ribbon_half_width;
					up = tangent.cross(side).normalized() * ribbon_thickness;

					const Vector3 v0 = a + side + up;
					const Vector3 v1 = a - side + up;
					const Vector3 v2 = b + side + up;
					const Vector3 v3 = b - side + up;

					const Vector3 normal = (side.cross(tangent)).normalized();

					st->set_normal(normal);
					st->add_vertex(v0);
					st->set_normal(normal);
					st->add_vertex(v1);
					st->set_normal(normal);
					st->add_vertex(v2);

					st->set_normal(normal);
					st->add_vertex(v2);
					st->set_normal(normal);
					st->add_vertex(v1);
					st->set_normal(normal);
					st->add_vertex(v3);
				}
			}
		}
	}

	Ref<ArrayMesh> ribbon_mesh = st->commit();
	if (ribbon_mesh.is_null()) {
		return;
	}

	ribbon_mesh_instance = memnew(MeshInstance3D);
	ribbon_mesh_instance->set_mesh(ribbon_mesh);
	add_child(ribbon_mesh_instance);
	ribbon_mesh_instance->set_visible(false);
}
