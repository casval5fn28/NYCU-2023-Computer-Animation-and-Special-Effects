#include "simulation/kinematics.h"

#include "Eigen/Dense"
#include <iostream>

#include "acclaim/bone.h"
#include "util/helper.h"

namespace kinematics {

    void BFS(const acclaim::Posture& posture, acclaim::Bone* bone, bool visited[]) {
        visited[bone->idx] = true;
        bone->start_position = bone->parent->end_position;

        bone->rotation = bone->parent->rotation * bone->rot_parent_current;
        bone->rotation *= util::rotateDegreeZYX(posture.bone_rotations[bone->idx]);

        bone->end_position = bone->start_position;
        bone->end_position += bone->rotation * (bone->dir.normalized() * bone->length);

        acclaim::Bone* tmp = bone->sibling;
        while (tmp != nullptr) {
            if (!visited[tmp->idx]) {
                BFS(posture, tmp, visited);
            }
            tmp = tmp->sibling;
        }

        if (bone->child != nullptr) {
            if (!visited[bone->child->idx]) {
                BFS(posture, bone->child, visited);
            }
        }
    }

    void forwardSolver(const acclaim::Posture& posture, acclaim::Bone* bone) {
        // TODO (FK)
        // Same as HW2
        // Hint:
        //   1. If you don't use `axis` in this function, you can copy-paste your code

        bool visited[31] = {};
        // root
        bone->start_position = posture.bone_translations[0];

        bone->rotation = bone->rot_parent_current;
        bone->rotation = util::rotateDegreeZYX(posture.bone_rotations[bone->idx]);

        bone->end_position = bone->start_position;
        bone->end_position += bone->rotation * (bone->dir.normalized() * bone->length);
        //cout<<"x:"<<posture.bone_rotations[bone->idx].x()<<" y:"<<posture.bone_rotations[bone->idx].y()<<" z:"<<posture.bone_rotations[bone->idx].z()<<" w:"<<posture.bone_rotations[bone->idx].w()<<endl;
        visited[0] = true;

        BFS(posture, bone->child, visited);
    }


    Eigen::VectorXd pseudoInverseLinearSolver(const Eigen::Matrix4Xd& Jacobian, const Eigen::Vector4d& target) {
        // TODO (find x which min(| jacobian * x - target |))
        // Hint:
        //   1. Linear algebra - least squares solution
        //   2. https://en.wikipedia.org/wiki/Moore%E2%80%93Penrose_inverse#Construction
        // Note:
        //   1. SVD or other pseudo-inverse method is useful
        //   2. Some of them have some limitation, if you use that method you should check it.
        Eigen::VectorXd deltatheta(Jacobian.cols());
        // get pseudo inverse of the Jacobian with SVD
        Eigen::JacobiSVD<Eigen::Matrix4Xd> svd(Jacobian, Eigen::ComputeThinU | Eigen::ComputeThinV);
        // minimum-norm solution
        deltatheta = svd.solve(target);

        //Eigen::Matrix3Xd J = Jacobian.topRows(3);
        //Eigen::Vector3d tgt = target.head<3>();
        //Eigen::VectorXd deltatheta = J.transpose() * (J * J.transpose()).inverse() * tgt;
        return deltatheta;  
    }

    /**
     * @brief Perform inverse kinematics (IK)
     *
     * @param target_pos The position where `end_bone` will move to.
     * @param start_bone This bone is the last bone you can move while doing IK
     * @param end_bone This bone will try to reach `target_pos`
     * @param posture The original AMC motion's reference, you need to modify this
     *
     * @return True if IK is stable (HW3 bonus)
     */
    bool inverseJacobianIKSolver(const Eigen::Vector4d& target_pos, acclaim::Bone* start_bone, acclaim::Bone* end_bone,
        acclaim::Posture& posture) {
        constexpr int max_iteration = 1000;
        constexpr double epsilon = 1E-3;
        constexpr double step = 0.1;
        // Since bone stores in bones[i] that i == bone->idx, we can use bone - bone->idx to find bones[0] which is the root.
        acclaim::Bone* root_bone = start_bone - start_bone->idx;
        // TODO
        // Perform inverse kinematics (IK)
        // HINTs will tell you what should do in that area.
        // Of course you can ignore it (Any code below this line) and write your own code.

        /////////////
        acclaim::Posture original_posture(posture);

        size_t bone_num = 1;// 0
        std::vector<acclaim::Bone*> boneList;

        acclaim::Bone* current = end_bone;
        while (current != start_bone && current != root_bone) {
            boneList.push_back(current);
            bone_num++;
            current = current->parent;
        }

        // TODO
        // Calculate number of bones need to move to perform IK, store in `bone_num` 
        // a.k.a. how may bones from end_bone to its parent then to start_bone (include both start_bone and end_bone)
        // Store the bones need to move to perform IK into boneList
        // Hint:
        //   1. Traverse from end_bone to start_bone is easier than start to end (since there is only 1 parent)
        //   2. If start bone is not reachable from end. Go to root first.
        // Note:
        //   1. Both start_bone and end_bone should be in the list
        //acclaim::Bone* current = end_bone;


        Eigen::Matrix4Xd Jacobian(4, 3 * bone_num);
        Jacobian.setZero();
        for (int iter = 0; iter < max_iteration; ++iter) {
            forwardSolver(posture, root_bone);
            Eigen::Vector4d desiredVector = target_pos - end_bone->end_position;
            if (desiredVector.norm() < epsilon) {
                break;
            }
            // TODO (compute jacobian)
            //   1. Compute arm vectors
            //   2. Compute jacobian columns, store in `Jacobian`
            // Hint:
            //   1. You should not put rotation in jacobian if it doesn't have that DoF.
            //   2. jacobian.col(/* some column index */) = /* jacobian column */

            current = end_bone;
            for (long long i = 0; i < bone_num; i++) {

                Eigen::Matrix3d rot_mat = current->rotation.linear();
                Eigen::Vector4d delta_pos = end_bone->end_position - current->start_position;

                for (int j = 0; j < 3; j++) {
                    Eigen::Vector3d ai = Eigen::Vector3d::Zero();
                    if (j == 0) {// dofrx == true
                        ai = (rot_mat * Eigen::Vector3d(1, 0, 0)).normalized();
                    }
                    else if (j == 1) {// dofry == true
                        ai = (rot_mat * Eigen::Vector3d(0, 1, 0)).normalized();
                    }
                    else {// dofrz == true
                        ai = (rot_mat * Eigen::Vector3d(0, 0, 1)).normalized();
                    }
                    Eigen::Vector3d radius = Eigen::Vector3d::Zero();
                    radius << delta_pos.head(3);
                    Eigen::Vector3d p_dif = ai.cross(radius);
                    //cout<<j<<": "<<par_dif.x()<<", "<<par_dif.y()<<", "<<par_dif.z()<<"\n";
                    Jacobian.col(i * 3 + j) << p_dif, 1;
                }
                current = current->parent;
            }

            Eigen::VectorXd deltatheta = step * pseudoInverseLinearSolver(Jacobian, desiredVector);

            // TODO (update rotation)
            //   Update `posture.bone_rotation` (in euler angle / degrees) using deltaTheta
            // Hint:
            //   1. You can ignore rotation limit of the bone.
            // Bonus:
            //   1. You cannot ignore rotation limit of the bone.

            current = end_bone;
            for (long long i = 0; i < bone_num; i++) {
                for (int j = 0; j < 3; j++) { // x, y, z
                    posture.bone_rotations[current->idx][j] += deltatheta[i * 3 + j];
                }
                current = current->parent;
            }

        }
        // TODO (Bonus)
        // Return whether IK is stable (i.e. whether the ball is reachable) and let the skeleton not swing its hand in the air
        return true;
    }
}  // namespace kinematics
