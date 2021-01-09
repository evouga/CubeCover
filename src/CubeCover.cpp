#include "CubeCover.h"
#include "TetMeshConnectivity.h"
#include "FrameField.h"
#include "Integration.h"
#include <iostream>

namespace CubeCover
{
    bool cubeCover(
        const Eigen::MatrixXd& V,
        const Eigen::MatrixXi& T,
        const Eigen::MatrixXd& frames,
        const Eigen::MatrixXi& assignments,
        Eigen::MatrixXd& parameterization,
        CubeCoverOptions opt
    )
    {
        TetMeshConnectivity mesh(T);
        if (!mesh.isManifold(opt.verbose))
            return false;
        if (!mesh.isFaceConnected())
        {
            if (opt.verbose)
            {
                std::cerr << "input mesh is not face-connected" << std::endl;
                return false;
            }
        }

        
        FrameField* field;
        if (opt.assignmentHandling == CubeCoverOptions::AssignmentHandling::AH_USEPROVIDED)
        {
            field = fromFramesAndAssignments(mesh, frames, assignments, opt.verbose);
        }
        else
        {
            Eigen::MatrixXi identity(0, 2 + frames.rows() / mesh.nTets());
            field = fromFramesAndAssignments(mesh, frames, identity, opt.verbose);
            field->computeLocalAssignments();
        }

        if (!field)
            return false;

        if (!integrate(V, *field, parameterization, opt.scale, opt.MIPtol, opt.parameterizationType,
            opt.boundaryConditions == CubeCoverOptions::BoundaryConditions::BC_FORCEINTEGER, opt.verbose))
        {
            delete field;
            return false;
        }

        delete field;
        return true;
    }

    bool cubeCover(
        const Eigen::MatrixXd& V,
        const Eigen::MatrixXi& T,
        const Eigen::MatrixXd& frames,
        const Eigen::MatrixXi& assignments,
        Eigen::MatrixXd& parameterization
    )
    {
        CubeCoverOptions opt;
        return cubeCover(V, T, frames, assignments, parameterization, opt);
    }

};