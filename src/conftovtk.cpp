//  Copyright or © or Copr. Rockable
//
//  vincent.richefeu@3sr-grenoble.fr
//
//  This software is a computer program whose purpose is
//    (i)  to hold sphero-polyhedral shapes,
//    (ii) to manage breakable interfaces.
//  It is developed for an ACADEMIC USAGE
//
//  This software is governed by the CeCILL-B license under French law and
//  abiding by the rules of distribution of free software.  You can  use,
//  modify and/ or redistribute the software under the terms of the CeCILL-B
//  license as circulated by CEA, CNRS and INRIA at the following URL
//  "http://www.cecill.info".
//
//  As a counterpart to the access to the source code and  rights to copy,
//  modify and redistribute granted by the license, users are provided only
//  with a limited warranty  and the software's author,  the holder of the
//  economic rights,  and the successive licensors  have only  limited
//  liability.
//
//  In this respect, the user's attention is drawn to the risks associated
//  with loading,  using,  modifying and/or developing or reproducing the
//  software by the user in light of its specific status of free software,
//  that may mean  that it is complicated to manipulate,  and  that  also
//  therefore means  that it is reserved for developers  and  experienced
//  professionals having in-depth computer knowledge. Users are therefore
//  encouraged to load and test the software's suitability as regards their
//  requirements in conditions enabling the security of their systems and/or
//  data to be ensured and,  more generally, to use and operate it in the
//  same conditions as regards security.
//
//  The fact that you are presently reading this means that you have had
//  knowledge of the CeCILL-B license and that you accept its terms.

#include "conftovtk.hpp"

bool tryToReadConf(int num) {
  char file_name[256];
  sprintf(file_name, "conf%d", num);
  if (fileTool::fileExists(file_name)) {
    std::cout << "Read " << file_name << std::endl;
    box.clearMemory();
    box.loadConf(file_name);
    complexityNumber = 0;
    for (size_t i = 0; i < box.Particles.size(); ++i) {
      complexityNumber += box.Particles[i].shape->vertex.size();
    }
    confNum = box.iconf;
    box.computeAABB();
  } else {
    std::cout << file_name << " does not exist" << std::endl;
    return false;
  }
  return true;
}

std::string TruncateLongString(std::string const str, unsigned int maxLength) {
  std::string str_return = str;

  if (str.length() > maxLength) return str_return.substr(0, maxLength);
  return str_return;
}

void SeparateParticlesByType() {
  Sphers.clear();
  Joncs.clear();
  Polyrs.clear();

  SphersId.clear();
  JoncsId.clear();
  PolyrsId.clear();

  for (size_t id = 0; id < box.Particles.size(); id++) {
    Particle* P = &box.Particles[id];

    if (P->shape->vertex.size() == 1) {
      Sphers.push_back(P);
      SphersId.push_back(id);
    }
    if (P->shape->vertex.size() == 2) {
      Joncs.push_back(P);
      JoncsId.push_back(id);
    }
    if (P->shape->vertex.size() > 2) {
      Polyrs.push_back(P);
      PolyrsId.push_back(id);
    }
  }
}

void writeVTKSPHER(int num) {
  unsigned int i;
  char vtk_file[200];
  FILE* sortie_vtk;

  unsigned int nbgrains = (unsigned int)Sphers.size();

  sprintf(vtk_file, "Shepres%.4i.vtk", num);
  sortie_vtk = fopen(vtk_file, "w");

  fprintf(sortie_vtk, "# vtk DataFile Version 3.0\n");
  fprintf(sortie_vtk, "Sortie Particles\n");
  fprintf(sortie_vtk, "ASCII\n");
  fprintf(sortie_vtk, "DATASET UNSTRUCTURED_GRID\n");
  fprintf(sortie_vtk, "POINTS %i float\n", nbgrains);

  for (i = 0; i < nbgrains; i++) {
    vec3r& X = Sphers[i]->pos;
    fprintf(sortie_vtk, "%lf %lf %lf\n", X.x, X.y, X.z);
  }

  fprintf(sortie_vtk, "POINT_DATA %i\n", nbgrains);

  fprintf(sortie_vtk, "VECTORS Radius float\n");

  for (i = 0; i < nbgrains; i++) {
    fprintf(sortie_vtk, "0. 0. %lf\n", Sphers[i]->MinskowskiRadius());
  }

  fprintf(sortie_vtk, "VECTORS Vel float\n");

  for (i = 0; i < nbgrains; i++) {
    vec3r& V = Sphers[i]->vel;
    fprintf(sortie_vtk, "%lf %lf %lf\n", V.x, V.y, V.z);
  }

  fprintf(sortie_vtk, "SCALARS Material float\n");
  fprintf(sortie_vtk, "LOOKUP_TABLE default\n");

  for (i = 0; i < nbgrains; i++) {
    fprintf(sortie_vtk, "%i\n", Sphers[i]->group);
  }

  fprintf(sortie_vtk, "SCALARS Id int\n");
  fprintf(sortie_vtk, "LOOKUP_TABLE default\n");

  for (i = 0; i < nbgrains; i++) {
    fprintf(sortie_vtk, "%li\n", SphersId[i]);
  }

  fclose(sortie_vtk);
}

void writeVTKPOLYR(int num) {
  char vtk_file[200];
  FILE* sortie_vtk;
  int vertex_compt = 0;
  int nb_points = 0;
  int nb_faces = 0;
  int nb_face_points = 0;
  unsigned int nbgrains = (unsigned int)Polyrs.size();
  for (size_t i = 0; i < nbgrains; i++) {
    nb_points += Polyrs[i]->shape->vertex.size();
    nb_faces += Polyrs[i]->shape->face.size();

    for (size_t j = 0; j < Polyrs[i]->shape->face.size(); j++) {
      nb_face_points += Polyrs[i]->shape->face[j].size();
    }
  }

  sprintf(vtk_file, "Polyr%.4i.vtk", num);

  sortie_vtk = fopen(vtk_file, "w");

  fprintf(sortie_vtk, "# vtk DataFile Version 3.0\n");
  fprintf(sortie_vtk, "RIGID      1\n");
  fprintf(sortie_vtk, "ASCII\n");
  fprintf(sortie_vtk, "DATASET POLYDATA\n");
  fprintf(sortie_vtk, "POINTS %i float\n", nb_points);

  // # Ecriture des Vertices
  for (size_t i = 0; i < nbgrains; i++) {
    Shape* shp = Polyrs[i]->shape;

    size_t nbv = shp->vertex.size();
    for (size_t j = 0; j < nbv; j++) {
      vec3r X = Polyrs[i]->Glob(shp->vertex[j]);
      vec3r n = (X - Polyrs[i]->pos);
      n.normalize();
      X += n * shp->radius * Polyrs[i]->homothety;
      fprintf(sortie_vtk, "%lf %lf %lf\n", X.x, X.y, X.z);
    }
  }

  // Ecriture des Faces (connectivity des vertices)
  fprintf(sortie_vtk, "POLYGONS %i %i\n", nb_faces, nb_face_points + nb_faces);

  for (size_t i = 0; i < nbgrains; i++) {
    Shape* shp = Polyrs[i]->shape;
    size_t nbf = shp->face.size();
    for (size_t j = 0; j < nbf; j++) {
      fprintf(sortie_vtk, "%li ", shp->face[j].size());
      for (size_t k = 0; k < shp->face[j].size(); k++) fprintf(sortie_vtk, "%li ", vertex_compt + shp->face[j][k]);
      fprintf(sortie_vtk, "\n");
    }
    vertex_compt += shp->vertex.size();
  }

  fprintf(sortie_vtk, "CELL_DATA %i \n", nb_faces);
  fprintf(sortie_vtk, "SCALARS Material int 1\n");
  fprintf(sortie_vtk, "LOOKUP_TABLE default\n");
  for (size_t i = 0; i < nbgrains; i++) {
    Shape* shp = Polyrs[i]->shape;
    size_t nbf = shp->face.size();
    for (size_t j = 0; j < nbf; j++) {
      fprintf(sortie_vtk, "%d\n", Polyrs[i]->group);
    }
  }

  fprintf(sortie_vtk, "SCALARS Id int 1\n");
  fprintf(sortie_vtk, "LOOKUP_TABLE default\n");
  for (size_t i = 0; i < nbgrains; i++) {
    Shape* shp = Polyrs[i]->shape;
    size_t nbf = shp->face.size();
    for (size_t j = 0; j < nbf; j++) {
      fprintf(sortie_vtk, "%li\n", PolyrsId[i]);
    }
  }

  fprintf(sortie_vtk, "VECTORS Vel float\n");
  for (size_t i = 0; i < nbgrains; i++) {
    vec3r& V = Polyrs[i]->vel;
    Shape* shp = Polyrs[i]->shape;
    size_t nbf = shp->face.size();
    for (size_t j = 0; j < nbf; j++) {
      fprintf(sortie_vtk, "%lf %lf %lf\n", V.x, V.y, V.z);
    }
  }

  fclose(sortie_vtk);
}

void writeVTKPOLYRVertex(int num) {
  int vertex_compt = 0;
  int nb_points = 0;
  int nb_faces = 0;
  int nb_face_points = 0;
  unsigned int nbgrains = (unsigned int)Polyrs.size();
  for (size_t i = 0; i < nbgrains; i++) {
    nb_points += Polyrs[i]->shape->vertex.size();
    nb_faces += Polyrs[i]->shape->face.size();

    for (size_t j = 0; j < Polyrs[i]->shape->face.size(); j++) {
      nb_face_points += Polyrs[i]->shape->face[j].size();
    }
  }

  char vtk_file[200];
  char vtk_file2[200];

  FILE* sortie_vtk;
  FILE* sortie_vtk2;

  sprintf(vtk_file, "Polyr%.4i.vtk", num);
  sortie_vtk = fopen(vtk_file, "w");
  fprintf(sortie_vtk, "# vtk DataFile Version 3.0\n");
  fprintf(sortie_vtk, "RIGID      1\n");
  fprintf(sortie_vtk, "ASCII\n");
  fprintf(sortie_vtk, "DATASET POLYDATA\n");
  fprintf(sortie_vtk, "POINTS %i float\n", nb_points);

  sprintf(vtk_file2, "PolyrSpheres%.4i.vtk", num);
  sortie_vtk2 = fopen(vtk_file2, "w");
  fprintf(sortie_vtk2, "# vtk DataFile Version 3.0\n");
  fprintf(sortie_vtk2, "Sortie Particles\n");
  fprintf(sortie_vtk2, "ASCII\n");
  fprintf(sortie_vtk2, "DATASET UNSTRUCTURED_GRID\n");
  fprintf(sortie_vtk2, "POINTS %i float\n", nb_points);

  // # Ecriture des Vertices
  for (size_t i = 0; i < nbgrains; i++) {
    Shape* shp = Polyrs[i]->shape;

    size_t nbv = shp->vertex.size();
    for (size_t j = 0; j < nbv; j++) {
      vec3r X = Polyrs[i]->Glob(shp->vertex[j]);
      fprintf(sortie_vtk2, "%lf %lf %lf\n", X.x, X.y, X.z);
      // vec3r n = (X - Polyrs[i]->pos);
      // n.normalize();
      // X += 2.*n*shp->radius*Polyrs[i]->homothety;
      fprintf(sortie_vtk, "%lf %lf %lf\n", X.x, X.y, X.z);
    }
  }

  fprintf(sortie_vtk2, "POINT_DATA %i\n", nb_points);
  fprintf(sortie_vtk2, "VECTORS Radius float\n");
  for (size_t i = 0; i < nbgrains; i++) {
    Shape* shp = Polyrs[i]->shape;

    size_t nbv = shp->vertex.size();
    for (size_t j = 0; j < nbv; j++) {
      fprintf(sortie_vtk2, "0. 0. %lf\n", shp->radius);
    }
  }

  fprintf(sortie_vtk2, "SCALARS Id int 1\n");
  fprintf(sortie_vtk2, "LOOKUP_TABLE default\n");
  for (size_t i = 0; i < nbgrains; i++) {
    Shape* shp = Polyrs[i]->shape;
    size_t nbv = shp->vertex.size();
    for (size_t j = 0; j < nbv; j++) {
      fprintf(sortie_vtk2, "%li\n", PolyrsId[i]);
    }
  }

  fclose(sortie_vtk2);

  // Ecriture des Faces (connectivity des vertices)
  fprintf(sortie_vtk, "POLYGONS %i %i\n", nb_faces, nb_face_points + nb_faces);

  for (size_t i = 0; i < nbgrains; i++) {
    Shape* shp = Polyrs[i]->shape;
    size_t nbf = shp->face.size();
    for (size_t j = 0; j < nbf; j++) {
      fprintf(sortie_vtk, "%li ", shp->face[j].size());
      for (size_t k = 0; k < shp->face[j].size(); k++) fprintf(sortie_vtk, "%li ", vertex_compt + shp->face[j][k]);
      fprintf(sortie_vtk, "\n");
    }
    vertex_compt += shp->vertex.size();
  }

  fprintf(sortie_vtk, "CELL_DATA %i \n", nb_faces);
  fprintf(sortie_vtk, "SCALARS Material int 1\n");
  fprintf(sortie_vtk, "LOOKUP_TABLE default\n");
  for (size_t i = 0; i < nbgrains; i++) {
    Shape* shp = Polyrs[i]->shape;
    size_t nbf = shp->face.size();
    for (size_t j = 0; j < nbf; j++) {
      fprintf(sortie_vtk, "%d\n", Polyrs[i]->group);
    }
  }

  fprintf(sortie_vtk, "SCALARS Id int 1\n");
  fprintf(sortie_vtk, "LOOKUP_TABLE default\n");
  for (size_t i = 0; i < nbgrains; i++) {
    Shape* shp = Polyrs[i]->shape;
    size_t nbf = shp->face.size();
    for (size_t j = 0; j < nbf; j++) {
      fprintf(sortie_vtk, "%li\n", PolyrsId[i]);
    }
  }

  fprintf(sortie_vtk, "VECTORS Vel float\n");
  for (size_t i = 0; i < nbgrains; i++) {
    vec3r& V = Polyrs[i]->vel;
    Shape* shp = Polyrs[i]->shape;
    size_t nbf = shp->face.size();
    for (size_t j = 0; j < nbf; j++) {
      fprintf(sortie_vtk, "%lf %lf %lf\n", V.x, V.y, V.z);
    }
  }

  fclose(sortie_vtk);
}

void writeVTKContactsSpheres(int num) {
  char vtk_file[200];
  FILE* sortie_vtk;

  sprintf(vtk_file, "Forces%.4i.vtk", num);
  sortie_vtk = fopen(vtk_file, "w");

  long nbContacts = 0;
  double fMax = 0;

  for (size_t k = 0; k < box.Interactions.size(); ++k) {
    std::set<Interaction>::iterator it = box.Interactions[k].begin();
    for (; it != box.Interactions[k].end(); ++it) {
      {
        if (it->dn > 0.0 && it->stick == nullptr)
          continue;
        else {
          nbContacts++;
          fMax = std::max(fMax, norm(it->fn * it->n + it->ft));
        }
      }
    }
  }

  fprintf(sortie_vtk, "# vtk DataFile Version 3.0\n");
  fprintf(sortie_vtk, "Sortie Particles\n");
  fprintf(sortie_vtk, "ASCII\n");
  fprintf(sortie_vtk, "DATASET UNSTRUCTURED_GRID\n");
  fprintf(sortie_vtk, "POINTS %li float\n", nbContacts);

  for (size_t k = 0; k < box.Interactions.size(); ++k) {
    std::set<Interaction>::iterator it = box.Interactions[k].begin();
    for (; it != box.Interactions[k].end(); ++it) {
      {
        if (it->dn > 0.0 && it->stick == nullptr)
          continue;
        else {
          vec3r X = it->pos;
          fprintf(sortie_vtk, "%lf %lf %lf\n", X.x, X.y, X.z);
        }
      }
    }
  }

  fprintf(sortie_vtk, "POINT_DATA %li\n", nbContacts);

  fprintf(sortie_vtk, "VECTORS F_Fmax float\n");

  for (size_t k = 0; k < box.Interactions.size(); ++k) {
    std::set<Interaction>::iterator it = box.Interactions[k].begin();
    for (; it != box.Interactions[k].end(); ++it) {
      {
        int i = it->i;
        int j = it->j;

        if (it->dn > 0.0 && it->stick == nullptr)
          continue;
        else {
          fMax = box.Particles[i].MinskowskiRadius() + box.Particles[j].MinskowskiRadius();
          if (fMax != 0.)
            fprintf(sortie_vtk, "0. 0. %lf\n", fMax / 2.);
          else
            fprintf(sortie_vtk, "0. 0. 0.\n");
        }
      }
    }
  }

  fprintf(sortie_vtk, "VECTORS Fn float\n");

  for (size_t k = 0; k < box.Interactions.size(); ++k) {
    std::set<Interaction>::iterator it = box.Interactions[k].begin();
    for (; it != box.Interactions[k].end(); ++it) {
      {
        if (it->dn > 0.0 && it->stick == nullptr)
          continue;
        else {
          vec3r X = it->fn * it->n;
          fprintf(sortie_vtk, "%lf %lf %lf\n", X.x, X.y, X.z);
        }
      }
    }
  }

  fprintf(sortie_vtk, "VECTORS Ft float\n");

  for (size_t k = 0; k < box.Interactions.size(); ++k) {
    std::set<Interaction>::iterator it = box.Interactions[k].begin();
    for (; it != box.Interactions[k].end(); ++it) {
      {
        if (it->dn > 0.0 && it->stick == nullptr)
          continue;
        else {
          vec3r X = it->ft;
          fprintf(sortie_vtk, "%lf %lf %lf\n", X.x, X.y, X.z);
        }
      }
    }
  }

  fprintf(sortie_vtk, "SCALARS Type float\n");
  fprintf(sortie_vtk, "LOOKUP_TABLE default\n");

  for (size_t k = 0; k < box.Interactions.size(); ++k) {
    std::set<Interaction>::iterator it = box.Interactions[k].begin();
    for (; it != box.Interactions[k].end(); ++it) {
      {
        if (it->dn > 0.0 && it->stick == nullptr)
          continue;
        else {
          fprintf(sortie_vtk, "%d\n", it->type);
        }
      }
    }
  }

  fclose(sortie_vtk);
}

void writeVTKContactsLines(int num) {
  char vtk_file[200];
  FILE* sortie_vtk;

  sprintf(vtk_file, "ForceLines%.4i.vtk", num);
  sortie_vtk = fopen(vtk_file, "w");

  size_t nbContacts = 0;

  for (size_t k = 0; k < box.Interactions.size(); ++k) {
    std::set<Interaction>::iterator it = box.Interactions[k].begin();
    for (; it != box.Interactions[k].end(); ++it) {
      {
        if (it->dn > 0.0 && it->stick == nullptr) continue;

        if (it->type == vvType) nbContacts += 2;
        if (it->type == veType) nbContacts++;
        if (it->type == eeType) nbContacts += 2;
        if (it->type == vfType) nbContacts++;
      }
    }
  }

  fprintf(sortie_vtk, "# vtk DataFile Version 3.0\n");
  fprintf(sortie_vtk, "Sortie Foces\n");
  fprintf(sortie_vtk, "ASCII\n");
  fprintf(sortie_vtk, "DATASET POLYDATA\n");

  fprintf(sortie_vtk, "POINTS %li float\n", 2 * nbContacts);
  for (size_t k = 0; k < box.Interactions.size(); ++k) {
    std::set<Interaction>::iterator it = box.Interactions[k].begin();
    for (; it != box.Interactions[k].end(); ++it) {
      {
        int i = it->i;
        int j = it->j;

        if (it->dn > 0.0 && it->stick == nullptr) continue;

        if (it->type == vvType) {
          fprintf(sortie_vtk, "%e %e %e\n", it->pos.x, it->pos.y, it->pos.z);
          fprintf(sortie_vtk, "%e %e %e\n", box.Particles[i].pos.x, box.Particles[i].pos.y, box.Particles[i].pos.z);

          fprintf(sortie_vtk, "%e %e %e\n", it->pos.x, it->pos.y, it->pos.z);
          fprintf(sortie_vtk, "%e %e %e\n", box.Particles[j].pos.x, box.Particles[j].pos.y, box.Particles[j].pos.z);
        }

        if (it->type == eeType) {
          vec3r dnorm = it->n * (box.Particles[i].MinskowskiRadius());
          fprintf(sortie_vtk, "%e %e %e\n", it->pos.x, it->pos.y, it->pos.z);
          fprintf(sortie_vtk, "%e %e %e\n", it->pos.x + dnorm.x, it->pos.y + dnorm.y, it->pos.z + dnorm.z);

          dnorm = it->n * (box.Particles[j].MinskowskiRadius());
          fprintf(sortie_vtk, "%e %e %e\n", it->pos.x, it->pos.y, it->pos.z);
          fprintf(sortie_vtk, "%e %e %e\n", it->pos.x - dnorm.x, it->pos.y - dnorm.y, it->pos.z - dnorm.z);
        }

        if (it->type == veType || it->type == vfType) {
          fprintf(sortie_vtk, "%e %e %e\n", it->pos.x, it->pos.y, it->pos.z);
          fprintf(sortie_vtk, "%e %e %e\n", box.Particles[i].pos.x, box.Particles[i].pos.y, box.Particles[i].pos.z);
        }
      }
    }
  }

  fprintf(sortie_vtk, "LINES %li %li\n", nbContacts, 3 * nbContacts);

  for (size_t i = 0; i < nbContacts; i++) {
    fprintf(sortie_vtk, "2 %li %li\n", 2 * i, 2 * i + 1);
  }

  fprintf(sortie_vtk, "POINT_DATA %li \n", 2 * nbContacts);
  fprintf(sortie_vtk, "SCALARS Fn float\n");
  fprintf(sortie_vtk, "LOOKUP_TABLE default\n");

  for (size_t k = 0; k < box.Interactions.size(); ++k) {
    std::set<Interaction>::iterator it = box.Interactions[k].begin();
    for (; it != box.Interactions[k].end(); ++it) {
      {
        if (it->dn > 0.0 && it->stick == nullptr) continue;

        if (it->type == vvType || it->type == eeType) {
          fprintf(sortie_vtk, "%e\n", fabs(it->fn));
          fprintf(sortie_vtk, "%e\n", fabs(it->fn));
          fprintf(sortie_vtk, "%e\n", fabs(it->fn));
          fprintf(sortie_vtk, "%e\n", fabs(it->fn));
        }
        if (it->type == veType || it->type == vfType) {
          fprintf(sortie_vtk, "%e\n", fabs(it->fn));
          fprintf(sortie_vtk, "%e\n", fabs(it->fn));
        }
      }
    }
  }

  fprintf(sortie_vtk, "SCALARS Ft float\n");
  fprintf(sortie_vtk, "LOOKUP_TABLE default\n");

  for (size_t k = 0; k < box.Interactions.size(); ++k) {
    std::set<Interaction>::iterator it = box.Interactions[k].begin();
    for (; it != box.Interactions[k].end(); ++it) {
      {
        if (it->dn > 0.0 && it->stick == nullptr) continue;
        if (it->type == vvType || it->type == eeType) {
          fprintf(sortie_vtk, "%e\n", norm(it->ft));
          fprintf(sortie_vtk, "%e\n", norm(it->ft));
          fprintf(sortie_vtk, "%e\n", norm(it->ft));
          fprintf(sortie_vtk, "%e\n", norm(it->ft));
        }
        if (it->type == veType || it->type == vfType) {
          fprintf(sortie_vtk, "%e\n", norm(it->ft));
          fprintf(sortie_vtk, "%e\n", norm(it->ft));
        }
      }
    }
  }

  fclose(sortie_vtk);
}

void writeVTKOBB(int num) {
  char vtk_file[200];
  FILE* sortie_vtk;

  sprintf(vtk_file, "OBB%.4i.vtk", num);
  sortie_vtk = fopen(vtk_file, "w");

  size_t nbSommets = 8 * box.Particles.size();

  fprintf(sortie_vtk, "# vtk DataFile Version 3.0\n");
  fprintf(sortie_vtk, "Sortie Foces\n");
  fprintf(sortie_vtk, "ASCII\n");
  fprintf(sortie_vtk, "DATASET POLYDATA\n");

  fprintf(sortie_vtk, "POINTS %li float\n", nbSommets);

  OBB obbi;

  for (size_t i = 0; i < box.Particles.size(); i++) {

    obbi = box.Particles[i].shape->obb;
    obbi.rotate(box.Particles[i].Q);
    obbi.extent *= box.Particles[i].homothety;
    obbi.center *= box.Particles[i].homothety;
    obbi.center += box.Particles[i].pos;
    // if (enlarged_obb) obbi.enlarge(0.5 * box.DVerlet);

    // je vais écrire les 8 sommets
    vec3r corner;

    // les 4 sommets du bas
    corner = obbi.center - obbi.extent[0] * obbi.e[0] - obbi.extent[1] * obbi.e[1] - obbi.extent[2] * obbi.e[2];
    fprintf(sortie_vtk, "%e %e %e\n", corner.x, corner.y, corner.z);
    corner += 2.0 * obbi.extent[0] * obbi.e[0];
    fprintf(sortie_vtk, "%e %e %e\n", corner.x, corner.y, corner.z);
    corner += 2.0 * obbi.extent[1] * obbi.e[1];
    fprintf(sortie_vtk, "%e %e %e\n", corner.x, corner.y, corner.z);
    corner -= 2.0 * obbi.extent[0] * obbi.e[0];
    fprintf(sortie_vtk, "%e %e %e\n", corner.x, corner.y, corner.z);

    // les 4 sommets du haut
    corner = obbi.center - obbi.extent[0] * obbi.e[0] - obbi.extent[1] * obbi.e[1] - obbi.extent[2] * obbi.e[2];
    corner += 2.0 * obbi.extent[2] * obbi.e[2];
    fprintf(sortie_vtk, "%e %e %e\n", corner.x, corner.y, corner.z);
    corner += 2.0 * obbi.extent[0] * obbi.e[0];
    fprintf(sortie_vtk, "%e %e %e\n", corner.x, corner.y, corner.z);
    corner += 2.0 * obbi.extent[1] * obbi.e[1];
    fprintf(sortie_vtk, "%e %e %e\n", corner.x, corner.y, corner.z);
    corner -= 2.0 * obbi.extent[0] * obbi.e[0];
    fprintf(sortie_vtk, "%e %e %e\n", corner.x, corner.y, corner.z);
  }

  fprintf(sortie_vtk, "LINES %li %li\n", 6 * box.Particles.size(), 24 * box.Particles.size());

  long indexParticle = 0;

  for (size_t i = 0; i < box.Particles.size(); i++) {
    fprintf(sortie_vtk, "5 %li %li  %li %li %li\n", indexParticle, indexParticle + 1, indexParticle + 2,
            indexParticle + 3, indexParticle);
    fprintf(sortie_vtk, "5 %li %li  %li %li %li\n", indexParticle + 4, indexParticle + 5, indexParticle + 6,
            indexParticle + 7, indexParticle + 4);

    fprintf(sortie_vtk, "2 %li %li\n", indexParticle, indexParticle + 4);
    fprintf(sortie_vtk, "2 %li %li\n", indexParticle + 1, indexParticle + 5);
    fprintf(sortie_vtk, "2 %li %li\n", indexParticle + 2, indexParticle + 6);
    fprintf(sortie_vtk, "2 %li %li\n", indexParticle + 3, indexParticle + 7);
    indexParticle += 8;
  }

  fprintf(sortie_vtk, "CELL_DATA %li \n", 6 * box.Particles.size());
  fprintf(sortie_vtk, "SCALARS Id int 1\n");
  fprintf(sortie_vtk, "LOOKUP_TABLE default\n");
  for (size_t i = 0; i < box.Particles.size(); i++) {
    for (int j = 0; j < 6; j++) {
      fprintf(sortie_vtk, "%li\n", i);
    }
  }

  fclose(sortie_vtk);
}

void nbConfToVtk() {
  std::string nomfich;
  char nom[25];
  struct dirent* lecture;
  DIR* rep;
  rep = opendir(".");

  int number;

  tabConf.clear();

  while ((lecture = readdir(rep))) {
    strcpy(nom, lecture->d_name);

    nomfich = nom;

    if (nomfich[0] == 'c' && nomfich[1] == 'o' && nomfich[2] == 'n' && nomfich[3] == 'f') {
      std::string temp;

      for (unsigned int i = 0; i < nomfich.size(); i++) {
        if (isdigit(nomfich[i])) {
          temp += nomfich[i];
        }
      }
      std::istringstream stream(temp);
      stream >> number;
      tabConf.push_back(number);
    }
  }
  closedir(rep);

  // tirer :
  for (size_t i = 0; i < tabConf.size(); i++) {
    for (size_t j = i + 1; j < tabConf.size(); j++) {
      if (tabConf[j] < tabConf[i]) {
        int aux = tabConf[i];
        tabConf[i] = tabConf[j];
        tabConf[j] = aux;
      }
    }
  }
}

// =====================================================================
// Main function
// =====================================================================

int main(int argc, char* argv[]) {
  box.initParser();
  box.setInteractive(true);

  std::string confFileName;

  try {
    TCLAP::CmdLine cmd("VTK maker of Rockable simulations", ' ', "0.3");
    TCLAP::UnlabeledValueArg<std::string> nameArg("input", "Name of the conf-file", false, "conf0", "conf-file");

    cmd.add(nameArg);

    cmd.parse(argc, argv);

    confFileName = nameArg.getValue();
  } catch (TCLAP::ArgException& e) {
    std::cerr << "error: " << e.error() << " for argument " << e.argId() << std::endl;
  }

  nbConfToVtk();

  for (size_t i = 0; i < tabConf.size(); i++) {
    box.clearMemory();
    tryToReadConf(tabConf[i]);
    box.computeAABB();
    box.System.read();
    confNum = box.iconf;

    if (box.Particles.empty()) {
      std::cerr << "No particles!" << std::endl;
    } else {
      SeparateParticlesByType();
      if (Sphers.size() > 0) writeVTKSPHER(confNum);
      if (Polyrs.size() > 0) {
        writeVTKPOLYR(confNum);
      }
      writeVTKContactsSpheres(confNum);
      // writeVTKContactsLines(confNum);
      writeVTKOBB(confNum);
    }
  }
  return 0;
}