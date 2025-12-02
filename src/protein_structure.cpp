#include "protein_structure.h"

#include <algorithm>
#include <cmath>
#include <unordered_map>

#include <gemmi/elem.hpp>
#include <gemmi/metadata.hpp>
#include <gemmi/seqid.hpp>
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
	BIND_ENUM_CONSTANT(REPRESENTATION_BALL_AND_STICK);
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

float ProteinStructure::get_covalent_radius(const gemmi::Element &el) const {
	const float cov_r = el.covalent_r();
	if (cov_r > 0.0f) {
		return cov_r;
	}
	// Fallback to a light default when data is missing.
	return 0.7f;
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

Color ProteinStructure::get_secondary_color(SecondaryType type) const {
	switch (type) {
		case SecondaryType::HELIX:
			return Color(0.8f, 0.1f, 0.2f);
		case SecondaryType::SHEET:
			return Color(0.95f, 0.8f, 0.2f);
		case SecondaryType::COIL:
		default:
			return Color(0.6f, 0.6f, 0.7f);
	}
}

float ProteinStructure::get_secondary_width(SecondaryType type) const {
	switch (type) {
		case SecondaryType::HELIX:
			return 0.8f * SCALE_FACTOR;
		case SecondaryType::SHEET:
			return 1.2f * SCALE_FACTOR;
		case SecondaryType::COIL:
		default:
			return 0.6f * SCALE_FACTOR;
	}
}

float ProteinStructure::get_secondary_thickness(SecondaryType type) const {
	switch (type) {
		case SecondaryType::HELIX:
			return 0.22f * SCALE_FACTOR;
		case SecondaryType::SHEET:
			return 0.18f * SCALE_FACTOR;
		case SecondaryType::COIL:
		default:
			return 0.15f * SCALE_FACTOR;
	}
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

std::string ProteinStructure::residue_key(const std::string &chain_name, const gemmi::SeqId &seqid) const {
	return chain_name + ":" + seqid.str();
}

void ProteinStructure::build_secondary_lookup() {
	secondary_lookup.clear();

	auto seq_in_range = [](const gemmi::SeqId &value, const gemmi::SeqId &start, const gemmi::SeqId &end) {
		return !(value < start) && !(end < value);
	};

	for (const gemmi::Helix &helix : model_data.helices) {
		for (const gemmi::Model &model : model_data.models) {
			for (const gemmi::Chain &chain : model.chains) {
				if (chain.name != helix.start.chain_name) {
					continue;
				}
				for (const gemmi::Residue &residue : chain.residues) {
					if (seq_in_range(residue.seqid, helix.start.res_id.seqid, helix.end.res_id.seqid)) {
						secondary_lookup[residue_key(chain.name, residue.seqid)] = SecondaryType::HELIX;
					}
				}
			}
		}
	}

	for (const gemmi::Sheet &sheet : model_data.sheets) {
		for (const gemmi::Sheet::Strand &strand : sheet.strands) {
			for (const gemmi::Model &model : model_data.models) {
				for (const gemmi::Chain &chain : model.chains) {
					if (chain.name != strand.start.chain_name) {
						continue;
					}
					for (const gemmi::Residue &residue : chain.residues) {
						if (!seq_in_range(residue.seqid, strand.start.res_id.seqid, strand.end.res_id.seqid)) {
							continue;
						}

						const std::string key = residue_key(chain.name, residue.seqid);
						if (secondary_lookup.find(key) == secondary_lookup.end()) {
							secondary_lookup[key] = SecondaryType::SHEET;
						}
					}
				}
			}
		}
	}
}

ProteinStructure::SecondaryType ProteinStructure::get_secondary_for_residue(const std::string &chain_name, const gemmi::SeqId &seqid) const {
	const std::string key = residue_key(chain_name, seqid);
	const auto found = secondary_lookup.find(key);
	if (found != secondary_lookup.end()) {
		return found->second;
	}
	return SecondaryType::COIL;
}

Vector3 ProteinStructure::compute_residue_normal(const gemmi::Residue &residue, const Vector3 &tangent_hint) const {
	Vector3 c_pos;
	Vector3 o_pos;
	Vector3 ca_pos;
	bool has_c = false;
	bool has_o = false;
	bool has_ca = false;

	for (const gemmi::Atom &atom : residue.atoms) {
		if (atom.name == "C") {
			c_pos = to_vec3(atom.pos);
			has_c = true;
		} else if (atom.name == "O") {
			o_pos = to_vec3(atom.pos);
			has_o = true;
		} else if (atom.name == "CA") {
			ca_pos = to_vec3(atom.pos);
			has_ca = true;
		}
	}

	if (has_c && has_o) {
		Vector3 dir = (o_pos - c_pos);
		if (dir.length_squared() > CMP_EPSILON) {
			return dir.normalized();
		}
	}

	if (has_ca && has_c) {
		Vector3 dir = (c_pos - ca_pos);
		if (dir.length_squared() > CMP_EPSILON) {
			return dir.normalized();
		}
	}

	if (tangent_hint.length_squared() > CMP_EPSILON) {
		Vector3 dir = tangent_hint.cross(Vector3(0, 1, 0));
		if (dir.length_squared() < CMP_EPSILON) {
			dir = tangent_hint.cross(Vector3(1, 0, 0));
		}
		if (dir.length_squared() > CMP_EPSILON) {
			return dir.normalized();
		}
	}

	return Vector3(0, 1, 0);
}

Ref<StandardMaterial3D> ProteinStructure::make_standard_material(bool double_sided, float roughness, float metallic) const {
	Ref<StandardMaterial3D> mat;
	mat.instantiate();
	mat->set_shading_mode(StandardMaterial3D::SHADING_MODE_PER_PIXEL);
	mat->set_flag(BaseMaterial3D::FLAG_ALBEDO_FROM_VERTEX_COLOR, true);
	mat->set_roughness(roughness);
	mat->set_metallic(metallic);
	mat->set_transparency(StandardMaterial3D::TRANSPARENCY_DISABLED);
	mat->set_cull_mode(double_sided ? BaseMaterial3D::CULL_DISABLED : BaseMaterial3D::CULL_BACK);
	return mat;
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
	if (ball_mesh_instance != nullptr) {
		ball_mesh_instance->queue_free();
		ball_mesh_instance = nullptr;
	}
	if (ribbon_mesh_instance != nullptr) {
		ribbon_mesh_instance->queue_free();
		ribbon_mesh_instance = nullptr;
	}
	bond_segments.clear();
	secondary_lookup.clear();
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
		model_data = gemmi::read_structure_file(native_path);
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
	build_secondary_lookup();
	compute_bonds();
	make_atom_multimesh();
	make_ballstick_multimesh();
	make_bond_multimesh();
	make_cartoon_mesh();
	show_current_mode();
}

void ProteinStructure::show_current_mode() {
	if (atom_mesh_instance != nullptr) {
		atom_mesh_instance->set_visible(current_mode == REPRESENTATION_ATOMS);
	}
	if (ball_mesh_instance != nullptr) {
		ball_mesh_instance->set_visible(current_mode == REPRESENTATION_BALL_AND_STICK);
	}
	if (bond_mesh_instance != nullptr) {
		bond_mesh_instance->set_visible(current_mode == REPRESENTATION_STICKS || current_mode == REPRESENTATION_BALL_AND_STICK);
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
		case REPRESENTATION_BALL_AND_STICK:
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

void ProteinStructure::compute_bonds() {
	bond_segments.clear();

	struct RenderAtom {
		Vector3 pos;
		gemmi::Element element;

		RenderAtom(const Vector3 &p, const gemmi::Element &el) : pos(p), element(el) {}
	};

	std::vector<RenderAtom> atoms;
	atoms.reserve(1024);

	for (const gemmi::Model &model : model_data.models) {
		for (const gemmi::Chain &chain : model.chains) {
			for (const gemmi::Residue &residue : chain.residues) {
				for (const gemmi::Atom &atom : residue.atoms) {
					atoms.emplace_back(to_vec3(atom.pos), atom.element);
				}
			}
		}
	}

	const int atom_count = static_cast<int>(atoms.size());
	for (int i = 0; i < atom_count; ++i) {
		const float cov_r1 = get_covalent_radius(atoms[i].element);
		for (int j = i + 1; j < atom_count; ++j) {
			const float cov_r2 = get_covalent_radius(atoms[j].element);
			const float max_dist = cov_r1 + cov_r2 + BOND_TOLERANCE;

			const Vector3 delta = atoms[i].pos - atoms[j].pos;
			const float dist_sq = delta.length_squared();
			if (dist_sq <= CMP_EPSILON) {
				continue;
			}

			if (dist_sq <= max_dist * max_dist) {
				bond_segments.emplace_back(atoms[i].pos * SCALE_FACTOR, atoms[j].pos * SCALE_FACTOR);
			}
		}
	}
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
	sphere->set_radial_segments(24);
	sphere->set_rings(16);
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
	atom_mesh_instance->set_material_override(make_standard_material(false));
	add_child(atom_mesh_instance);
}

void ProteinStructure::make_ballstick_multimesh() {
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
	sphere->set_radial_segments(24);
	sphere->set_rings(16);
	multimesh->set_mesh(sphere);

	int idx = 0;
	for (const gemmi::Model &model : model_data.models) {
		for (const gemmi::Chain &chain : model.chains) {
			for (const gemmi::Residue &residue : chain.residues) {
				for (const gemmi::Atom &atom : residue.atoms) {
					const std::string symbol(atom.element.name());
					const float radius = get_vdw_radius(symbol) * SCALE_FACTOR * BALL_STICK_RADIUS_SCALE;
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

	ball_mesh_instance = memnew(MultiMeshInstance3D);
	ball_mesh_instance->set_multimesh(multimesh);
	ball_mesh_instance->set_material_override(make_standard_material(false, 0.25f, 0.05f));
	add_child(ball_mesh_instance);
	ball_mesh_instance->set_visible(false);
}

void ProteinStructure::make_bond_multimesh() {
	if (!ensure_structure_loaded()) {
		return;
	}

	if (bond_segments.empty()) {
		compute_bonds();
	}

	if (bond_segments.empty()) {
		return;
	}

	Ref<MultiMesh> multimesh;
	multimesh.instantiate();
	multimesh->set_transform_format(MultiMesh::TRANSFORM_3D);
	multimesh->set_use_colors(true);
	multimesh->set_instance_count(static_cast<int>(bond_segments.size()));

	Ref<CylinderMesh> cylinder;
	cylinder.instantiate();
	cylinder->set_height(1.0f);
	cylinder->set_top_radius(0.5f);
	cylinder->set_bottom_radius(0.5f);
	multimesh->set_mesh(cylinder);

	for (int i = 0; i < static_cast<int>(bond_segments.size()); ++i) {
		const Vector3 a = bond_segments[i].first;
		const Vector3 b = bond_segments[i].second;
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
	bond_mesh_instance->set_material_override(make_standard_material(false, 0.2f, 0.08f));
	add_child(bond_mesh_instance);
	bond_mesh_instance->set_visible(false);
}

void ProteinStructure::make_cartoon_mesh() {
	if (!ensure_structure_loaded()) {
		return;
	}

	struct Slice {
		Vector3 center;
		Vector3 side;
		Vector3 normal;
		Color color;
		float width = 0.0f;
		float thickness = 0.0f;
	};

	int vertex_cursor = 0;

	Ref<SurfaceTool> st;
	st.instantiate();
	st->begin(Mesh::PRIMITIVE_TRIANGLES);
	st->set_smooth_group(-1);
	st->set_uv(Vector2(0, 0));

	auto catmull_rom_tangent = [](const Vector3 &p0, const Vector3 &p1, const Vector3 &p2, const Vector3 &p3, float t) {
		const float t2 = t * t;
		return 0.5f * ((p2 - p0) + (2.0f * (2.0f * p0 - 5.0f * p1 + 4.0f * p2 - p3)) * t +
					   (3.0f * (-p0 + 3.0f * p1 - 3.0f * p2 + p3)) * t2);
	};

	// helper to store vertex indices for strip stitching
	struct SliceIndices {
		int v0; // +side +up
		int v1; // -side +up
		int v2; // +side -up
		int v3; // -side -up
	};

	for (const gemmi::Model &model : model_data.models) {
		for (const gemmi::Chain &chain : model.chains) {
			struct CartoonResidue {
				Vector3 pos;
				Vector3 normal;
				SecondaryType sec = SecondaryType::COIL;
				const gemmi::Residue *res_ptr = nullptr;
			};

			std::vector<CartoonResidue> residues;
			residues.reserve(chain.residues.size());

			for (const gemmi::Residue &residue : chain.residues) {
				Vector3 ca_pos;
				bool has_ca = false;
				for (const gemmi::Atom &atom : residue.atoms) {
					if (atom.name == "CA") {
						ca_pos = to_vec3(atom.pos) * SCALE_FACTOR;
						has_ca = true;
						break;
					}
				}

				if (!has_ca) {
					continue;
				}

				CartoonResidue entry;
				entry.pos = ca_pos;
				entry.sec = get_secondary_for_residue(chain.name, residue.seqid);
				entry.res_ptr = &residue;
				residues.push_back(entry);
			}

			if (residues.size() < 2) {
				continue;
			}

			for (size_t i = 0; i < residues.size(); ++i) {
				const Vector3 prev = i > 0 ? residues[i - 1].pos : residues[i].pos;
				const Vector3 next = i + 1 < residues.size() ? residues[i + 1].pos : residues[i].pos;
				const Vector3 tangent_hint = next - prev;

				Vector3 normal = compute_residue_normal(*residues[i].res_ptr, tangent_hint);
				if (i > 0 && residues[i - 1].normal.length_squared() > CMP_EPSILON && normal.length_squared() > CMP_EPSILON) {
					if (normal.dot(residues[i - 1].normal) < 0.0f) {
						normal = -normal;
					}
				}
				if (normal.length_squared() < CMP_EPSILON) {
					normal = Vector3(0, 1, 0);
				}
				residues[i].normal = normal.normalized();
			}

			// Pad endpoints for Catmull-Rom.
			std::vector<Vector3> padded;
			padded.reserve(residues.size() + 2);
			padded.push_back(residues.front().pos);
			for (const CartoonResidue &res : residues) {
				padded.push_back(res.pos);
			}
			padded.push_back(residues.back().pos);

			std::vector<Slice> samples;

			for (int i = 0; i < static_cast<int>(residues.size()) - 1; ++i) {
				const Vector3 p0 = padded[std::max(0, i)];
				const Vector3 p1 = padded[i + 1];
				const Vector3 p2 = padded[i + 2];
				const Vector3 p3 = padded[std::min(static_cast<int>(padded.size()) - 1, i + 3)];

				const Vector3 normal0 = residues[i].normal;
				const Vector3 normal1 = residues[i + 1].normal;
				const float width0 = get_secondary_width(residues[i].sec);
				const float width1 = get_secondary_width(residues[i + 1].sec);
				const float thickness0 = get_secondary_thickness(residues[i].sec);
				const float thickness1 = get_secondary_thickness(residues[i + 1].sec);
				const Color color0 = get_secondary_color(residues[i].sec);
				const Color color1 = get_secondary_color(residues[i + 1].sec);

				for (int step = 0; step <= CARTOON_SUBDIVISIONS; ++step) {
					const float t = static_cast<float>(step) / CARTOON_SUBDIVISIONS;

					const Vector3 center = catmull_rom(p0, p1, p2, p3, t);
					Vector3 normal = normal0.lerp(normal1, t).normalized();
					if (normal.length_squared() < CMP_EPSILON) {
						normal = Vector3(0, 1, 0);
					}

					Vector3 tangent = catmull_rom_tangent(p0, p1, p2, p3, t).normalized();
					if (tangent.length_squared() < CMP_EPSILON) {
						continue;
					}

					Vector3 side = tangent.cross(normal);
					if (side.length_squared() < CMP_EPSILON) {
						side = tangent.cross(Vector3(0, 1, 0));
					}
					if (side.length_squared() < CMP_EPSILON) {
						side = tangent.cross(Vector3(1, 0, 0));
					}
					side.normalize();
					normal = side.cross(tangent).normalized();

					Slice slice;
					slice.center = center;
					slice.side = side;
					slice.normal = normal;
					slice.width = Math::lerp(width0, width1, t);
					slice.thickness = Math::lerp(thickness0, thickness1, t);
					slice.color = color0.lerp(color1, t);
					samples.push_back(slice);
				}
			}

			if (samples.size() < 2) {
				continue;
			}

			std::vector<SliceIndices> indices;
			indices.reserve(samples.size());

			for (const Slice &slice : samples) {
				const Vector3 side_offset = slice.side * slice.width;
				const Vector3 up_offset = slice.normal * (slice.thickness * 0.5f);

				const Vector3 v0 = slice.center + side_offset + up_offset;
				const Vector3 v1 = slice.center - side_offset + up_offset;
				const Vector3 v2 = slice.center + side_offset - up_offset;
				const Vector3 v3 = slice.center - side_offset - up_offset;

				SliceIndices idxs;
				st->set_color(slice.color);
				idxs.v0 = vertex_cursor++;
				st->add_vertex(v0);
				st->set_color(slice.color);
				idxs.v1 = vertex_cursor++;
				st->add_vertex(v1);
				st->set_color(slice.color);
				idxs.v2 = vertex_cursor++;
				st->add_vertex(v2);
				st->set_color(slice.color);
				idxs.v3 = vertex_cursor++;
				st->add_vertex(v3);
				indices.push_back(idxs);
			}

			for (size_t i = 0; i + 1 < indices.size(); ++i) {
				const SliceIndices &a = indices[i];
				const SliceIndices &b = indices[i + 1];

				// Top
				st->add_index(a.v0);
				st->add_index(a.v1);
				st->add_index(b.v0);
				st->add_index(b.v0);
				st->add_index(a.v1);
				st->add_index(b.v1);

				// Bottom
				st->add_index(a.v2);
				st->add_index(a.v3);
				st->add_index(b.v2);
				st->add_index(b.v2);
				st->add_index(a.v3);
				st->add_index(b.v3);

				// Positive side
				st->add_index(a.v0);
				st->add_index(a.v2);
				st->add_index(b.v0);
				st->add_index(b.v0);
				st->add_index(a.v2);
				st->add_index(b.v2);

				// Negative side
				st->add_index(a.v1);
				st->add_index(b.v1);
				st->add_index(a.v3);
				st->add_index(a.v3);
				st->add_index(b.v1);
				st->add_index(b.v3);
			}
		}
	}

	st->index();
	st->generate_normals();
	st->generate_tangents();
	Ref<ArrayMesh> ribbon_mesh = st->commit();
	if (ribbon_mesh.is_null()) {
		return;
	}

	ribbon_mesh_instance = memnew(MeshInstance3D);
	ribbon_mesh_instance->set_mesh(ribbon_mesh);
	ribbon_mesh_instance->set_material_override(make_standard_material(true, 0.28f, 0.03f));
	add_child(ribbon_mesh_instance);
	ribbon_mesh_instance->set_visible(false);
}
