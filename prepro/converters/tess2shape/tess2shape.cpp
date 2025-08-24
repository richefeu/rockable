#include "tess2shape.hpp"

/**
 * @brief Convert a polyhTool::poly into a Shape.
 *
 * This function takes a polyhTool::poly and fills a Shape with the same
 * vertices, edges and faces. The Shape is cleared before filling.
 *
 * @param [in] block : the polyhTool::poly to convert.
 * @param [in/out] shape : the Shape to fill.
 *
 * @return 1 if everything is okay, 0 otherwise.
 */
int Tess2Shape::PolyToShape(const polyhTool::poly& block, Shape& shape) {
  shape.vertex.clear();
  shape.edge.clear();
  shape.face.clear();

  size_t nv = block.nodes.size();
  size_t ne = block.edges.size();
  size_t nf = block.faces.size();

  // vetexes
  for (size_t iv = 0; iv < nv; iv++) {
    vec3r v = vec3r();
    v = block.nodes[iv].pos;
    shape.vertex.push_back(v);
  }

  // edges
  for (size_t ie = 0; ie < ne; ie++) {
    std::pair<size_t, size_t> e;
    e.first = block.edges[ie].node0;
    e.second = block.edges[ie].node1;
    shape.edge.push_back(e);
  }

  // faces
  for (size_t ifc = 0; ifc < nf; ifc++) {
    polyhTool::pface f = block.faces[ifc];
    std::vector<size_t> ids;
    for (auto it = f.nodes.begin(); it != f.nodes.end(); it++) {
      size_t node_face = *it;
      ids.push_back(node_face);
    }
    shape.face.push_back(ids);
  }

  return 1;
}

/**
 * @brief Build a cube as a polyhTool::poly.
 *
 * This function creates a polyhTool::poly from a cube with the given
 * min and max vertices.
 *
 * @param [in] min_v : the minimum vertex of the cube (bottom left front).
 * @param [in] max_v : the maximum vertex of the cube (top right back).
 *
 * @return a polyhTool::poly representing the cube.
 */
polyhTool::poly Tess2Shape::cube(vec3r min_v, vec3r max_v) {
  polyhTool::poly p;

  p.nodes.resize(8);
  p.edges.resize(12);
  p.faces.resize(6);

  p.nodes[0].pos = vec3r(min_v.x, min_v.y, min_v.z);
  p.nodes[1].pos = vec3r(max_v.x, min_v.y, min_v.z);
  p.nodes[2].pos = vec3r(max_v.x, min_v.y, max_v.z);
  p.nodes[3].pos = vec3r(min_v.x, min_v.y, max_v.z);
  p.nodes[4].pos = vec3r(min_v.x, max_v.y, min_v.z);
  p.nodes[5].pos = vec3r(max_v.x, max_v.y, min_v.z);
  p.nodes[6].pos = vec3r(max_v.x, max_v.y, max_v.z);
  p.nodes[7].pos = vec3r(min_v.x, max_v.y, max_v.z);

  p.edges[0].node0 = 0;
  p.edges[0].node1 = 1;
  p.edges[1].node0 = 1;
  p.edges[1].node1 = 2;
  p.edges[2].node0 = 2;
  p.edges[2].node1 = 3;
  p.edges[3].node0 = 3;
  p.edges[3].node1 = 0;
  p.edges[4].node0 = 4;
  p.edges[4].node1 = 5;
  p.edges[5].node0 = 5;
  p.edges[5].node1 = 6;
  p.edges[6].node0 = 6;
  p.edges[6].node1 = 7;
  p.edges[7].node0 = 7;
  p.edges[7].node1 = 4;
  p.edges[8].node0 = 0;
  p.edges[8].node1 = 4;
  p.edges[9].node0 = 5;
  p.edges[9].node1 = 1;
  p.edges[10].node0 = 2;
  p.edges[10].node1 = 6;
  p.edges[11].node0 = 3;
  p.edges[11].node1 = 7;

  p.faces[0].nodes = {0, 1, 2, 3};
  p.faces[1].nodes = {4, 5, 6, 7};
  p.faces[2].nodes = {0, 1, 5, 4};
  p.faces[3].nodes = {1, 5, 6, 2};
  p.faces[4].nodes = {2, 6, 7, 3};
  p.faces[5].nodes = {0, 4, 7, 3};

  return p;
}

/**
 * @brief Read a tessellation file and fill the corresponding members of the
 *        class.
 *
 * This function reads a file in the tessellation format and fills the
 * following members of the class:
 *   - number_of_cells
 *   - seed_positions
 *   - vertex
 *   - edge
 *   - face
 *   - polyhedron
 *
 * The file is assumed to be in the tessellation format, i.e. it starts with
 * a line containing the string "**cell", followed by a line containing the
 * number of cells, then a line containing the string "**seed", followed by
 * the positions of the seeds, and so on.
 *
 * The function reports to the standard output the number of cells, vertices,
 * edges and faces read from the file.
 *
 * @return nothing
 */
void Tess2Shape::readTess() {
  std::cout << "> Read tesselation\n";
  std::ifstream file(tessFileName.c_str(), std::ios::in);

  if (!file.is_open()) {
    std::cerr << "@Tess2Shape::readTess, cannot open the file" << tessFileName << std::endl;
  }

  std::string line;
  int tmp_int;
  double tmp_double;

  // Cells
  do {
    std::getline(file, line);
  } while (line != std::string(" **cell"));

  file >> number_of_cells;
  std::cout << "  number_of_cells " << number_of_cells << std::endl;

  // Positions
  do {
    std::getline(file, line);
  } while (line != std::string("  *seed"));

  seed_positions.resize(number_of_cells);
  for (size_t i = 0; i < seed_positions.size(); i++) {
    file >> tmp_int >> seed_positions[i].x >> seed_positions[i].y >> seed_positions[i].z >> tmp_double;
  }

  // Vertices
  do {
    std::getline(file, line);
  } while (line != std::string(" **vertex"));

  file >> tmp_int;
  std::cout << "  number of vertices " << tmp_int << std::endl;
  vertex.resize(tmp_int);
  for (size_t i = 0; i < vertex.size(); i++) {
    file >> vertex[i].id >> vertex[i].v.x >> vertex[i].v.y >> vertex[i].v.z >> tmp_int;
  }

  // Edges
  do {
    std::getline(file, line);
  } while (line != std::string(" **edge"));

  file >> tmp_int;
  std::cout << "  number of edges " << tmp_int << std::endl;
  edge.resize(tmp_int);
  for (size_t i = 0; i < edge.size(); i++) {
    file >> edge[i].id >> edge[i].ver_1 >> edge[i].ver_2 >> tmp_int;
  }

  // Face
  do {
    std::getline(file, line);
  } while (line != std::string(" **face"));

  file >> tmp_int;
  std::cout << "  number of faces " << tmp_int << std::endl;
  face.resize(tmp_int);
  for (size_t i = 0; i < face.size(); i++) {
    file >> face[i].id;
    file >> tmp_int;               // number of vertices
    face[i].ver_.resize(tmp_int);  // cout << face[i].ver_.size() << std::endl;
    for (size_t j = 0; j < face[i].ver_.size(); j++) {
      file >> face[i].ver_[j];
    }

    file >> tmp_int;
    face[i].edge_.resize(tmp_int);
    for (size_t j = 0; j < face[i].edge_.size(); j++) {
      file >> face[i].edge_[j];
    }

    file >> face[i].face_eq_d >> face[i].face_eq_a >> face[i].face_eq_b >> face[i].face_eq_c;
    file >> tmp_int >> tmp_int;
    file >> tmp_double >> tmp_double >> tmp_double;
  }

  // Polyhedron
  do {
    std::getline(file, line);
  } while (line != std::string(" **polyhedron"));

  file >> tmp_int;
  std::cout << "  number of polyhedra " << tmp_int << std::endl;
  polyhedron.resize(tmp_int);
  for (size_t i = 0; i < polyhedron.size(); i++) {
    file >> tmp_int;
    file >> tmp_int;  // number_of_faces
    // std::cout << "     polyhedron " << i << " has " << tmp_int << " faces" << std::endl;
    polyhedron[i].face_.resize(tmp_int);
    for (size_t j = 0; j < polyhedron[i].face_.size(); j++) {
      file >> polyhedron[i].face_[j];
    }
  }
}

/**
 * @brief Create the Rockable-shapes from the polyhedra.
 *
 * This function iterates over the Shapes that have been previously extracted
 * and it removes a slice of thickness 'Minskowski radius'.
 *
 *
 * @return nothing
 */
void Tess2Shape::cutShapes() {
  std::cout << "> Cut polyhedra\n";
  for (size_t is = 0; is < Shapes.size(); is++) {
    size_t nsub = 0;
    Shape& shp = Shapes[is];
    AABB aabb(shp.vertex);
    polyhTool::poly block = cube(aabb.min, aabb.max);

    polyhTool::poly sub_block;

    for (size_t ip = 0; ip < shape_plans[is].size(); ip++) {
      polyhTool::plan p = shape_plans[is][ip];
      p.normal.normalize();
      p.normal *= -1.0;
      p.pos += radius * p.normal;
      if (polyhTool::cut_poly(block, p, sub_block) == 1) {
        block = sub_block;
        nsub++;
      }
    }

    if (nsub != 0) {
      PolyToShape(block, shp);
    }
  }
}

/**
 * @brief Create a list of Rockable-shapes from the polyhedron list.
 *
 * Each Shape is associated with a polyhedron and contains the vertices, edges
 * and faces of the polyhedron. The Shape also contains the shape plans, i.e.
 * the list of planes that define the polyhedron.
 *
 * @return nothing
 */
void Tess2Shape::extractShapes() {
  std::cout << "> Transfert polyhedron cells to Rockable-shapes\n";

  size_t np = polyhedron.size();

  Shapes.resize(np);
  face_plans.resize(np);
  shape_plans.resize(np);

  for (size_t ip = 0; ip < polyhedron.size(); ip++) {
    Shape& shp_id = Shapes[ip];
    //vec3r& position = seed_positions[ip];
    std::vector<plan_>& face_plan = face_plans[ip];
    std::vector<polyhTool::plan>& shape_plan = shape_plans[ip];

    shp_id.name = "shape_" + std::to_string(ip);
    shp_id.radius = radius;
    shp_id.volume = 1.0;
    shp_id.inertia_mass = vec3r(1., 1., 1.);
    shp_id.obb.e[0] = vec3r::unit_x();
    shp_id.obb.e[1] = vec3r::unit_y();
    shp_id.obb.e[2] = vec3r::unit_z();
    shp_id.orientation = quat(0., 0., 0., 1.0);
    shp_id.position.reset();
    shp_id.MCnstep = MCnstep;
    shp_id.isSurface = false;

    std::vector<vertex_> V;

    // Cherche, classe et supprime les doublons dans la liste des vertex
    for (size_t if1 = 0; if1 < polyhedron[ip].face_.size(); if1++) {  // loop over faces
      int face_id = abs(polyhedron[ip].face_[if1]);
      face_ ff;
      ff.id = face_id;
      int face_index = getIndex(face, ff);
      for (size_t iv1 = 0; iv1 < face[face_index].ver_.size(); iv1++) {  // loop over vertex
        vertex_ vv;
        vv.id = face[face_index].ver_[iv1];
        int vertex_index = getIndex(vertex, vv);
        V.push_back(vertex[vertex_index]);
      }
    }

    sort(V.begin(), V.end());
    V.erase(unique(V.begin(), V.end()), V.end());  // remove pairs of identical data

    // Cherche la liste des edges
    std::vector<edge_> E;                                             // Edge id
    for (size_t if1 = 0; if1 < polyhedron[ip].face_.size(); if1++) {  // loop over faces
      int face_id = abs(polyhedron[ip].face_[if1]);
      face_ ff;
      ff.id = face_id;
      int face_index = getIndex(face, ff);
      for (size_t ie1 = 0; ie1 < face[face_index].edge_.size(); ie1++) {  // loop over vertex
        edge_ ee;
        ee.id = face[face_index].edge_[ie1];
        int edge_index = getIndex(edge, ee);
        E.push_back(edge[edge_index]);
      }
    }

    sort(E.begin(), E.end());
    E.erase(unique(E.begin(), E.end()), E.end());  // supprime les doublons

    shp_id.vertex.resize(V.size());

    for (size_t iv1 = 0; iv1 < V.size(); iv1++) {
      vec3r v = vec3r();
      v = V[iv1].v;
      shp_id.vertex[iv1] = v;
    }

     AABB aabb(shp_id.vertex);
     vec3r extent = 0.5*(aabb.max - aabb.min);
     vec3r center = 0.5*(aabb.max + aabb.min);
     shp_id.obb.extent = extent;
     shp_id.obb.center = center;

    shp_id.edge.resize(E.size());

    for (size_t ie1 = 0; ie1 < E.size(); ie1++) {
      std::pair<size_t, size_t> e;

      vertex_ vv1;
      vv1.id = E[ie1].ver_1;
      vertex_ vv2;
      vv2.id = E[ie1].ver_2;
      int vertex1_index = getIndex(V, vv1);
      int vertex2_index = getIndex(V, vv2);

      e.first = vertex1_index;
      e.second = vertex2_index;
      shp_id.edge[ie1] = e;
    }

    shp_id.face.resize(polyhedron[ip].face_.size());
    face_plan.resize(shp_id.face.size());
    shape_plan.resize(shp_id.face.size());

    for (size_t if1 = 0; if1 < polyhedron[ip].face_.size(); if1++) {
      // loop over faces
      int face_id = polyhedron[ip].face_[if1];

      face_ ff;
      ff.id = abs(face_id);
      int face_index = getIndex(face, ff);

      shp_id.face[if1].resize(face[face_index].ver_.size());

      face_plan[if1].d = face[face_index].face_eq_d;
      face_plan[if1].normal.x = face[face_index].face_eq_a;
      face_plan[if1].normal.y = face[face_index].face_eq_b;
      face_plan[if1].normal.z = face[face_index].face_eq_c;

      shape_plan[if1].normal = face_plan[if1].normal;
      if (face_id < 0) {
        shape_plan[if1].normal *= -1.0;
      }

      vec3r pos = vec3r();
      for (size_t iv1 = 0; iv1 < face[face_index].ver_.size(); iv1++) {  // loop over vertex
        vertex_ vv;
        vv.id = face[face_index].ver_[iv1];
        int vertex_index = getIndex(V, vv);
        shp_id.face[if1][iv1] = vertex_index;
        pos += V[vertex_index].v;
      }

      int nv = (int)face[face_index].ver_.size();
      if (nv != 0) {
        pos /= (double)nv;
      }

      shape_plan[if1].pos = pos;
    }
  }

}

/**
 * @brief Read a Command
 *
 *
 * @param [in] name : the name of the command file to read.
 *
 * @return nothing.
 */
void Tess2Shape::readCommands(const char* name) {
  std::ifstream Com(name);
  if (!Com.is_open()) {
    std::cerr << msg::warn() << "@Tess2Shape::readCommands, Cannot read " << name << msg::normal() << std::endl;
    return;
  }

  kwParser parser;

  parser.kwMap["tessFileName"] = __GET__(Com, tessFileName);
  parser.kwMap["inputFileName"] = __GET__(Com, inputFileName);
  parser.kwMap["shapesFileName"] = __GET__(Com, shapesFileName);

  parser.kwMap["MinskowskiRadius"] = __GET__(Com, radius);

  parser.kwMap["ParticlesGroup"] = __GET__(Com, ParticlesGroup);
  parser.kwMap["ParticlesCluster"] = __GET__(Com, ParticlesCluster);
  parser.kwMap["ParticlesHomothety"] = __GET__(Com, ParticlesHomothety);

  parser.kwMap["MCnstep"] = __GET__(Com, MCnstep);

  // This single line actually parses the file
  parser.parse(Com);
}

/**
 * @brief Write a Rockable shape file.
 *
 * This function writes a Rockable shape file
 * containing the shapes defined in the Shapes vector.
 *
 * @return nothing
 */
void Tess2Shape::writeRockableShapeFile() {
  std::cout << "> Writting '" << shapesFileName << "'\n";

  std::ofstream shape_file(shapesFileName.c_str());
  for (size_t i = 0; i < Shapes.size(); i++) {
    Shapes[i].write(shape_file);
  }
}

/**
 * @brief Write the Rockable particles data to an output stream.
 *
 * This function iterates over all the Shapes and writes their properties to
 * the provided output stream in a specific format. Each line in the output
 * represents a Shape's particle data including its name, particle group,
 * cluster, homothety, position, orientation, and other fixed values.
 *
 * @param os The output stream to which the particle data is written.
 */

void Tess2Shape::writeRockableParticles(std::ostream& os) {
  for (size_t i = 0; i < Shapes.size(); i++) {
    os << Shapes[i].name << ' ' << ParticlesGroup << ' ' << ParticlesCluster << ' ' << ParticlesHomothety << ' '
       << Shapes[i].position * ParticlesHomothety << ' ' << "0 0 0" << ' ' << "0 0 0" << ' ' << Shapes[i].orientation << ' ' << "0 0 0"
       << ' ' << "0 0 0" << '\n';
  }
}

/**
 * @brief Write a Rockable configuration file.
 *
 * This function writes a Rockable configuration file (with extension .txt)
 * containing the particles defined in the Shapes vector.
 *
 * The function only writes the following lines to the file:
 *   - 'Rockable 29-11-2018'
 *   - 'shapeFile <shapesFileName>'
 *   - 'Particles <number of particles>'
 *   - <particles definitions>
 *
 * The particles definitions are written by calling writeRockableParticles.
 *
 * @return nothing
 */
void Tess2Shape::writeRockableInputFile() {
  std::cout << "> Writting '" << inputFileName << "'\n";

  std::ofstream conf(inputFileName.c_str());
  conf << "Rockable 29-11-2018\n";
  conf << "shapeFile " << shapesFileName << '\n';
  conf << "Particles " << Shapes.size() << '\n';
  writeRockableParticles(conf);
}

/**
 * @brief Main function of the Tess2Shape class.
 *
 * This function is the main entry point of the Tess2Shape class. It reads a
 * tessellation file, extracts the Shapes from it, cuts the Shapes, precomputes
 * everythings and writes out the Rockable configuration file and the shape
 * file.
 *
 * @return nothing
 */
void Tess2Shape::exec() {
  readTess();
  extractShapes();
  cutShapes();

  // Precompute everythings
  for (size_t i = 0; i < Shapes.size(); i++) {
    Shapes[i].massProperties();
    Shapes[i].preCompDone = 'y';
  }

  writeRockableShapeFile();
  writeRockableInputFile();
}
