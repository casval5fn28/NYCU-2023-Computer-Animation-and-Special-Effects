#include "simulation/kinematics.h"

#include "Eigen/Dense"
#include <iostream>
#include "acclaim/bone.h"
#include "util/helper.h"
using namespace std;

void BFS(const acclaim::Posture& posture, acclaim::Bone* bone, bool visited[]){
    visited[bone->idx] = true;
    bone->start_position = bone->parent->end_position;

    bone->rotation = bone->parent->rotation * bone->rot_parent_current;
    bone->rotation *= util::rotateDegreeZYX(posture.bone_rotations[bone->idx]);

    bone->end_position = bone->start_position;
    bone->end_position += bone->rotation * (bone->dir.normalized() * bone->length);

    acclaim::Bone* tmp = bone->sibling;
    while (tmp != nullptr){
        if (!visited[tmp->idx]){
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

namespace kinematics {
void forwardSolver(const acclaim::Posture& posture, acclaim::Bone* bone) {
    // TODO#1 (FK)
    // You should set these variables:
    //     bone->start_position = Eigen::Vector4d::Zero();
    //     bone->end_position = Eigen::Vector4d::Zero();
    //     bone->rotation = Eigen::Matrix4d::Zero();
    // The sample above just set everything to zero
    // Hint:
    //   1. posture.bone_translations, posture.bone_rotations
    // Note:
    //   1. This function will be called with bone == root bone of the skeleton
    //   2. we use 4D vector to represent 3D vector, so keep the last dimension as "0"
    //   3. util::rotate{Degree | Radian} {XYZ | ZYX}
    //      e.g. rotateDegreeXYZ(x, y, z) means:
    //      x, y and z are presented in degree rotate z degrees along z - axis first, then y degrees along y - axis, and then x degrees along x - axis 

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
    //bone->start_position = Eigen::Vector4d::Zero();
    //bone->end_position = Eigen::Vector4d::Zero();
    //bone->rotation = Eigen::Matrix4d::Zero();

}

std::vector<acclaim::Posture> timeWarper(const std::vector<acclaim::Posture>& postures, int allframe_old, int allframe_new) {

    int total_frames = static_cast<int>(postures.size());
    int total_bones = static_cast<int>(postures[0].bone_rotations.size());
    std::vector<acclaim::Posture> new_postures;
    for (int i = 0; i <= allframe_new; ++i) {
        acclaim::Posture new_poseture(total_bones);
        for (int j = 0; j < total_bones; ++j) {

            // TODO#2 (Time warping)
            // original: |--------------|
            // new     : |----------------------|
            // OR
            // original: |--------------|
            // new     : |-------|
            // You should set these variables:
            //     new_postures[i].bone_translations[j] = Eigen::Vector4d::Zero();
            //     new_postures[i].bone_rotations[j] = Eigen::Vector4d::Zero();
            // The sample above just set everything to zero
            // Hint:
            //   1. Scale the frames.
            //   2. You can use linear interpolation with translations.
            //   3. You should use spherical linear interpolation for rotations.
            
            double scaling = double(allframe_old) / double(allframe_new);
            int floor = (i * scaling >= total_frames - 1) ? total_frames - 1 : i * scaling;
            int top = (floor == total_frames - 1) ? floor : floor + 1;
            double ratio = i * scaling - floor;

            // translation
            Eigen::Vector4d translation = postures[floor].bone_translations[j];
            Eigen::Vector4d translation_new = postures[top].bone_translations[j];

            // rotation
            Eigen::Quaterniond rotation;
            rotation.w() = postures[floor].bone_rotations[j].w();
            rotation.x() = postures[floor].bone_rotations[j].x();
            rotation.y() = postures[floor].bone_rotations[j].y();
            rotation.z() = postures[floor].bone_rotations[j].z();
            
            Eigen::Quaterniond rotation_new;
            rotation_new.w() = postures[top].bone_rotations[j].w();
            rotation_new.x() = postures[top].bone_rotations[j].x();
            rotation_new.y() = postures[top].bone_rotations[j].y();
            rotation_new.z() = postures[top].bone_rotations[j].z();
            
            // final input
            Eigen::Vector4d fin_translation = translation + (translation_new - translation) * ratio;
            new_poseture.bone_translations[j] = fin_translation;
            Eigen::Quaterniond fin_rotation = rotation.slerp(ratio, rotation_new);//取1-ratio會打回臉上
            new_poseture.bone_rotations[j] = Eigen::Vector4d(fin_rotation.x(), fin_rotation.y(), fin_rotation.z(), fin_rotation.w());

            //new_poseture.bone_translations[j] = Eigen::Vector4d::Zero();
            //new_poseture.bone_rotations[j] = Eigen::Vector4d::Zero();
        }

        new_postures.push_back(new_poseture);
    }
    return new_postures;
}
}  // namespace kinematics
