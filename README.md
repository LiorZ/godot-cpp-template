This is a comprehensive project specification document designed to be handed to a specialized coding agent (or an experienced C++ developer). It outlines the architecture, technical constraints, and implementation steps for creating a high-performance Protein Visualization GDExtension for Godot 4.

-----

# Project Specification: Godot Bio-Molecular Plugin (GDExtension)

## 1\. Executive Summary

**Goal:** specific Develop a high-performance Godot 4 GDExtension using C++ that acts as a bridge between raw biological structural data (PDB/mmCIF files) and real-time game engine visualization.
**Inspiration:** "Molecular Nodes" (Blender), but optimized for real-time runtime interaction (games/simulations) rather than offline rendering.
**Core Philosophy:** Parse once using industry-standard bioinformatics libraries (Gemmi), then render efficiently using GPU instancing (MultiMesh) and procedural geometry (SurfaceTool).

-----

## 2\. Technical Stack & Constraints

  * **Engine:** Godot 4.2+ (Stable).
  * **Language:** C++17 or C++20 (Required for GDExtension).
  * **Dependencies:**
      * **[godot-cpp](https://github.com/godotengine/godot-cpp):** Standard bindings for GDExtension.
      * **[gemmi](https://github.com/project-gemmi/gemmi):** A header-only (or static) C++ library for parsing PDB/CIF, symmetry, and structural biology math.
  * **Build System:** SCons.
  * **Target Platforms:** Windows, Linux, macOS (Architecture must support cross-compilation).

-----

## 3\. System Architecture

The system revolves around a single "manager" node that holds the data and manages child nodes for rendering.

### 3.1 Class Structure

**`ProteinStructure` (Inherits `godot::Node3D`)**
This is the main API class exposed to the Godot Editor.

  * **Private Members:**
      * `gemmi::Structure model_data`: The loaded biological data.
      * `MultiMeshInstance3D* atom_mesh_instance`: Handles Space-filling/Ball representations.
      * `MultiMeshInstance3D* bond_mesh_instance`: Handles Stick representations.
      * `MeshInstance3D* ribbon_mesh_instance`: Handles Ribbon/Cartoon generation.
      * `RepresentationMode current_mode`: Enum (Atoms, Sticks, Cartoon).
  * **Exposed Methods:**
      * `load_file(String path)`: Parses PDB/CIF.
      * `set_representation(int mode)`: Toggles visibility and generation logic.
      * `apply_rotamer(int chain_idx, int residue_idx, float chi_angle)`: Mutates data and triggers a visual rebuild.

-----

## 4\. Implementation Steps (End-to-End)

### Phase 1: Build System & Dependency Injection

**Objective:** Setup a compiling GDExtension environment that links Gemmi.

1.  **Folder Structure:**
    ```text
    /project_root
      /godot_project (The actual game project)
      /src (C++ source)
      /godot-cpp (Submodule)
      /gemmi (Submodule/Include)
      SConstruct (Build script)
    ```
2.  **SCons Configuration:**
      * Modify `SConstruct` to include the `gemmi/include` path.
      * Ensure standard C++ exceptions are enabled (Gemmi uses them).
      * *Note for Agent:* Gemmi is header-only by default, but building a static lib for it reduces compile times for the plugin.

### Phase 2: The Data Layer (Parsing)

**Objective:** Read a file and store it in memory using Gemmi.

1.  **Implement `load_file(String path)`:**
      * Convert Godot `String` to `std::string`.
      * Call `gemmi::read_structure(path)`.
      * **Crucial Step - Centering:** Calculate the Center of Mass of the protein and subtract this vector from all atom coordinates. PDB coordinates can be far from (0,0,0), causing floating-point jitter in Godot.
      * **Coordinate Scaling:** Map 1 Angstrom (PDB unit) to 0.1 or 1.0 Godot Units (meters). Define a `SCALE_FACTOR` constant.

### Phase 3: Visualization A - Atoms (Space Filling)

**Objective:** Render thousands of atoms with zero CPU overhead per frame.

1.  **Technique:** `MultiMeshInstance3D`.
2.  **Mesh:** Use a low-poly `SphereMesh` resource.
3.  **Algorithm:**
      * Create a `MultiMesh` with `transform_format = TRANSFORM_3D` and `use_colors = true`.
      * Set `instance_count` to `model_data.atoms.size()`.
      * Iterate through Gemmi atoms:
          * **Position:** `atom.pos` \* `SCALE_FACTOR`.
          * **Scale:** Determine radius based on Element (C=1.7A, H=1.2A, etc.).
          * **Color:** Use a standard CPK mapping (Carbon=Black/Grey, Oxygen=Red, Nitrogen=Blue, Sulfur=Yellow).
      * Populate the MultiMesh buffer.

### Phase 4: Visualization B - Sticks (Bonds)

**Objective:** Visualize connectivity.

1.  **Technique:** `MultiMeshInstance3D` with a `CylinderMesh`.
2.  **Bond Calculation:**
      * Use `gemmi::calculate_bonds` or Gemmi's neighbor search to find atoms within bonding distance (\~1.6 Angstroms).
      * *Constraint:* Only bond if atoms are in the same residue or adjacent residues (Peptide bond).
3.  **Transform Logic:**
      * For every bond (Atom A -\> Atom B):
          * **Position:** Midpoint between A and B.
          * **Orientation:** Use `Transform3D.looking_at()` so the cylinder Z-axis aligns with the vector (B - A).
          * **Height:** Distance(A, B).

### Phase 5: Visualization C - Cartoons (Ribbons)

**Objective:** Procedural mesh generation for the protein backbone. This is the most mathematically complex task.

1.  **Technique:** `SurfaceTool` (Godot API).
2.  **Data Extraction:**
      * Iterate the Gemmi model to find all Alpha-Carbon (CA) atoms in a continuous chain.
3.  **Spline Generation (Catmull-Rom):**
      * Do not connect CA atoms with straight lines. Generate interpolated points (e.g., 5-10 segments per residue).
      * Use `gemmi::Spline` or implement a standard cubic Hermite spline helper.
4.  **Geometry Extrusion:**
      * You cannot just extrude a tube; protein ribbons are flat.
      * **Orientation Vector:** You must calculate the "Peptide Plane" normal. A rough approximation for the "Up" vector at any CA point is the vector pointing toward the Oxygen atom of the backbone carbonyl group.
      * **Mesh Construction:** At each step along the spline, generate a "profile" (a flat rectangle or ellipse) oriented by the tangent (direction of chain) and the normal (Up vector). Bridge these profiles using `SurfaceTool.add_vertex`.

### Phase 6: Interaction (The "Game" Part)

**Objective:** Real-time modification.

1.  **Rotamer Switching:**
      * Implement `set_rotamer(residue_id, chi_1_angle)`.
      * Use Gemmi's topology functions to rotate the specific atoms of the sidechain relative to the backbone.
      * **Optimization:** Do not rebuild the entire MultiMesh. Update only the instances (transforms) corresponding to the modified atoms in the MultiMesh buffer.

-----

## 5\. Coding Guidelines for the Agent

1.  **Memory Safety:** Use smart pointers where possible, but remember Godot objects (`Node3D*`, `MultiMesh*`) are managed by Godot's internal reference counting or scene tree deletion. Do not `delete` a Node manually if it is in the tree; use `queue_free()`.
2.  **Error Handling:** PDB parsing fails often. Wrap Gemmi calls in `try/catch` blocks and use `godot::UtilityFunctions::printerr()` to report errors to the Godot console.
3.  **Threading:** Heavy parsing/meshing should eventually be moved to a `godot::WorkerThreadPool` task to avoid freezing the game during loading, but keep it on the main thread for the MVP (Minimum Viable Product).

## 6\. Deliverables

1.  `SConstruct` file configured for godot-cpp + gemmi.
2.  `register_types.cpp` (Entry point for GDExtension).
3.  `protein_structure.h` and `protein_structure.cpp` implementing the specification above.
4.  A compiled `.gdextension` file and dynamic library (`.dll`/`.so`) ready to drop into a Godot project.

## 7\. Reference Resources

  * **Gemmi Doc:** [https://gemmi.readthedocs.io/en/latest/](https://gemmi.readthedocs.io/en/latest/)
  * **Godot GDExtension C++:** [https://docs.godotengine.org/en/stable/tutorials/scripting/gdextension/index.html](https://docs.godotengine.org/en/stable/tutorials/scripting/gdextension/index.html)
  * **PDB Format:** [https://www.wwpdb.org/documentation/file-format-content/format33/v3.3.html](https://www.wwpdb.org/documentation/file-format-content/format33/v3.3.html)
