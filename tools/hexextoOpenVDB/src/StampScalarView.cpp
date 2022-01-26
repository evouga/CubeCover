#ifndef STAMPSCALARVIEW_H
#define STAMPSCALARVIEW_H

#include <Eigen/Core>
#include <vector>
#include <glm/glm.hpp>
#include <openvdb/openvdb.h>
#include "ReadHexEx.h"

#include "SceneInfo.h"
#include "FrameField.h"


#include "TetMeshConnectivity.h"

class SceneInfo;


void stampScalarView(SceneInfo info, Eigen::MatrixXd C_pertet );


#endif
