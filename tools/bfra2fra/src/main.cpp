#include <Eigen/Core>
#include <iostream>
#include <regex>

#include "ReadFrameField.h"
#include "FrameField.h"
#include "TetMeshConnectivity.h"
#include "readMeshFixed.h"
#include "WriteFrameField.h"

bool endsWith(const std::string& str, const std::string& suffix) {
  if (str.length() >= suffix.length()) {
    return str.compare(str.length() - suffix.length(), suffix.length(), suffix) == 0;
  } else {
    return false;
  }
}

int main(int argc, char *argv[]) {
  if (argc != 3) {
	std::cerr << "Usage: bfra2fra (.mesh) (.bfra)" << std::endl;
	return -1;
  }

  std::string meshfile = argv[1];
  std::string framefile = argv[2];
  bool is_bfra = endsWith(framefile, ".bfra");

  std::string permfile;
  Eigen::MatrixXi T, assignments, tet_F;
  Eigen::MatrixXd frames, V;

  if(!is_bfra) {
        std::cout << "not a bfra file" << std::endl;
  }

  if (!CubeCover::readMESH(meshfile, V, T, tet_F)) {
        std::cerr << "could not read .mesh file " << meshfile << std::endl;
        return -1;
  }


  if (CubeCover::readFrameField(framefile, permfile, T, frames, assignments, true)) {
        std::regex txt_regex("\\.bfra$");

        // Replace .txt with .obj
        std::string new_filename = std::regex_replace(framefile, txt_regex, ".fra");
        CubeCover::TetMeshConnectivity mesh(T);

        CubeCover::FrameField* field = CubeCover::fromFramesAndAssignments(mesh, frames, assignments, true);
        field->computeLocalAssignments();
        std::cout << "found " << field->nSingularEdges() << " singular edges" << std::endl;
        field->combAssignments();

        CubeCover::writeFrameField(new_filename, permfile, *field);

        // save bad verts
        std::string bad_verts = std::regex_replace(framefile, std::regex("\\.bfra"), "_badverts.txt");

        std::ofstream ofs(bad_verts);
        if (!ofs)
        {
          return false;
        }
        int nsing = field->nSingularEdges();
        ofs << "ids " << V.rows() << " " << 2 * nsing << std::endl;
        for (int i = 0; i < nsing; i++)
        {
          int edge = field->singularEdge(i);
          int v0 = mesh.edgeVertex(edge, 0);
          int v1 = mesh.edgeVertex(edge, 1);
          ofs << v0 << std::endl << v1 << std::endl;
        }
  }

  return 0;
}
