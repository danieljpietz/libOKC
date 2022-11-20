//
// Created by danie on 11/17/2022.
//

#ifndef CPP_NE_2_OKC_LINK_H
#define CPP_NE_2_OKC_LINK_H

#include <Eigen/Eigen>
#include <okc_math.h>
#include <okc_differential_types.h>

enum link_type_t {
    Revolute, Prismatic, Fixed
};

template <class T, size_t N>
class OKC;

template<class T, size_t N>
struct Link {

    size_t idx;
    const std::string name;
    const Link* parent = nullptr;
    const link_type_t link_type;
    const Eigen::Vector3<T> axis;


    const Frame<T, N> frame_o;
    Frame<T, N> frame_l;
    Frame<T, N> frame_g;

    DifferentialFrame<T, N> dx_frame_l;
    DifferentialFrame<T, N> ddx_frame_l;
    DifferentialFrame<T, N> dx_frame_g;
    DifferentialFrame<T, N> ddx_frame_g;

    const MassSpec<T> mass_spec;
    Dynamics_State<T, N> dynamics_state;
    Eigen::Matrix<T, 6, 6> local_mass_matrix;
    Eigen::Vector<T, 6> d_link_star;
    Differential_Dynamics_State<T, N> dx_dynamics_state;
    Differential_Dynamics_State<T, N> ddx_dynamics_state;

    void init() {

        frame_l.rotation = frame_o.rotation;
        frame_l.position = frame_o.position;

        switch (link_type) {
            case Revolute:
                frame_l.jacobian.template block<3, 1>(0, idx) = axis;
                break;
            case Prismatic:
                frame_l.jacobian.template block<3, 1>(3, idx) = axis;
                break;
            case Fixed:
                break;
        }

        if (parent == nullptr) {
            frame_g.rotation = frame_o.rotation;
            frame_g.position = frame_o.position;

            // TODO: Bug causes frame_g = frame_l to crash in debug if vectorization is enabled
            for (size_t j = 0; j < N; ++j) {
                for (size_t i = 0; i < 6; ++i) {
                frame_g.jacobian(i, j) = frame_l.jacobian(i, j);
                }
            }

        } else {
            assert(parent->idx < idx);
        }
    }

    void set(const Eigen::Vector<T, N> &x, const Eigen::Vector<T, N> &d_x) {
        local_kinematics(x, d_x);
        local_differential_kinematics(x, d_x);
        kinematics();
        differential_kinematics();
        dynamics(d_x);
        differential_dynamics(d_x);

    }

protected:

    void local_kinematics(const Eigen::Vector<T, N> &x, const Eigen::Vector<T, N> &d_x) {

        auto velocities = this->frame_l.jacobian * d_x;
        this->frame_l.ang_vel = {velocities[0], velocities[1], velocities[2]};
        this->frame_l.lin_vel = {velocities[3], velocities[4], velocities[5]};

        switch (this->link_type) {

            case Revolute:
                this->frame_l.rotation = this->frame_o.rotation * Eigen::AngleAxis<T>(x[this->idx], this->axis);
                break;

            case Prismatic:
                this->frame_l.position = this->frame_o.position + this->axis * d_x[this->idx];
                break;

            case Fixed:
                break;
        }

    }

    void local_differential_kinematics(const Eigen::Vector<T, N> &x, const Eigen::Vector<T, N> &d_x) {

        dx_frame_l.rotation[idx] = skew(Eigen::Vector3<T>(frame_l.jacobian.template block<3, 1>(0, idx))) * frame_l.rotation;
        dx_frame_l.position.col(idx) = frame_l.jacobian.template block<3, 1>(3, idx);

        ddx_frame_l.ang_vel.col(idx) = frame_l.jacobian.template block<3, 1>(0, idx);
        ddx_frame_l.lin_vel.col(idx) = frame_l.jacobian.template block<3, 1>(3, idx);

    }

    void differential_kinematics() {

        if (parent != nullptr) {
            child_differential_kinematics();
        } else {
            root_differential_kinematics();
        }

    }

    void root_differential_kinematics() {

        this->dx_frame_g.ang_vel = this->dx_frame_l.ang_vel;
        this->dx_frame_g.lin_vel = this->dx_frame_l.lin_vel;

        switch (this->link_type) {
            case Revolute:
                dx_frame_g.rotation = dx_frame_l.rotation;
                break;
            case Prismatic:
                assert(0);
                dx_frame_g.position = dx_frame_l.position;
                break;
            case Fixed:
                break;
        }

    }

    void child_differential_kinematics() {

        _differential_kinematics(parent->frame_g, parent->dx_frame_g, frame_l, dx_frame_l, frame_g, dx_frame_g, idx);


#if 0
        const Frame<T, N>& parent = this->parent->frame_g;
        const DifferentialFrame<T, N>& d_parent = this->parent->dx_frame_g;

        dx_frame_g.rotation = d_parent.rotation * frame_l.rotation;
        dx_frame_g.rotation[idx] = parent.rotation * dx_frame_l.rotation[0];

        for (size_t i = 0; i < N; ++i) {
            dx_frame_g.jacobian_pmap[i].template block<3, 3>(3, 0) = d_parent.rotation[i] * skew(frame_l.position);
        }

        if (link_type == Revolute)
            dx_frame_g.jacobian_pmap[idx].template block<3, 3>(0, 0) = dx_frame_l.rotation[0].transpose();
        else if (link_type == Prismatic)
            dx_frame_g.jacobian_pmap[idx].template block<3, 1>(3, idx) = d_parent.rotation[idx] * dx_frame_g.position.col(0);

        dx_frame_g.jacobian = dx_frame_g.jacobian_pmap * parent.jacobian + frame_g.jacobian_pmap * d_parent.jacobian;
        dx_frame_g.jacobian[idx].template block<3, 3>(3, 0) = d_parent.rotation[idx] * skew(Eigen::Vector3<T>(frame_l.jacobian.template block<3, 1>(3, 0)));

        auto j_p2map =  skew(parent.ang_vel) * skew(frame_l.position) + skew(frame_l.lin_vel);
        auto dx_j_p2map = skew(parent.ang_vel) * skew(Eigen::Vector3<T>(frame_l.jacobian.template block<3, 1>(3, 0)));
        auto ddx_j_p2map = skew(Eigen::Vector3<T>(frame_l.jacobian.template block<3, 1>(0, 0)))*skew(frame_l.position) + skew(Eigen::Vector3<T>(frame_l.jacobian.template block<3, 1>(0, 0)));

        dx_frame_g.dot_jacobian_pmap[idx].template block<3, 3>(0, 0) = -skew(frame_l.ang_vel) * (dx_frame_l.rotation[0].transpose());
        ddx_frame_g.dot_jacobian_pmap[idx].template block<3, 3>(0, 0) = -skew(frame_l.jacobian.template block<3, 1>(0, 0)) * (frame_l.rotation.transpose());

        for (size_t i = 0; i < idx; ++i) {
            dx_frame_g.dot_jacobian_pmap[i].template block<3, 3>(3, 0) = -d_parent.rotation[i] * j_p2map;
        }

        dx_frame_g.dot_jacobian_pmap[idx].template block<3, 3>(3, 0) = parent.rotation * dx_j_p2map;
        ddx_frame_g.dot_jacobian_pmap[idx].template block<3, 3>(3, 0) = parent.rotation * ddx_j_p2map;

        dx_frame_g.dot_jacobian = dx_frame_g.dot_jacobian_pmap * parent.jacobian + frame_g.dot_jacobian_pmap * d_parent.jacobian +
                dx_frame_g.jacobian * parent.dot_jacobian + frame_g.jacobian + d_parent.dot_jacobian;
        dx_frame_g.dot_jacobian[idx] = d_parent.rotation * parent.ang_vel;
#endif
    }

    void kinematics() {
        if (parent != nullptr) {
            child_kinematics();
        } else {
            root_kinematics();
        }
    }

    void child_kinematics() {

        auto parent_g = parent->frame_g;

        frame_g.rotation = parent_g.rotation * frame_l.rotation;
        frame_g.position = parent_g.position + parent_g.rotation * frame_l.position;

        frame_g.lin_vel = parent_g.lin_vel + skew(parent_g.ang_vel) * frame_g.position;
        frame_g.ang_vel = parent_g.ang_vel + frame_g.rotation * frame_l.ang_vel;

        frame_g.jacobian_pmap.template block<3, 3>(0, 0) = frame_l.rotation.transpose();
        frame_g.jacobian_pmap.template block<3, 3>(3, 0) = -parent_g.rotation * skew(frame_l.position);

        Eigen::Matrix<T, 6, N> jac_vec;
        jac_vec << frame_l.jacobian.template topRows<3>(), parent_g.rotation *
                                                           frame_l.jacobian.template bottomRows<3>();

        frame_g.jacobian = frame_g.jacobian_pmap * parent_g.jacobian + jac_vec;

        frame_g.dot_jacobian_pmap.template block<3, 3>(0, 0) = -skew(frame_l.ang_vel) * frame_l.rotation.transpose();
        frame_g.dot_jacobian_pmap.template block<3, 3>(3, 0) =
                -parent_g.rotation * (skew(parent_g.ang_vel) * skew(frame_l.position) + skew(frame_l.lin_vel));

        Eigen::Matrix<T, 6, N> d_jac_vec;

        d_jac_vec << frame_l.dot_jacobian.template topRows<3>(),
                parent_g.rotation * frame_l.dot_jacobian.template bottomRows<3>();

        frame_g.dot_jacobian = frame_g.dot_jacobian_pmap * parent_g.jacobian
                               + frame_g.jacobian_pmap * parent_g.dot_jacobian
                               + d_jac_vec;

    }

    void root_kinematics() {
        this->frame_g.ang_vel = this->frame_l.ang_vel;
        this->frame_g.lin_vel = this->frame_l.lin_vel;
        switch (this->link_type) {
            case Revolute:
                frame_g.rotation = frame_l.rotation;
                break;
            case Prismatic:
                frame_g.position = frame_l.position;
                break;
            case Fixed:
                break;
        }
    }

    void dynamics(const Eigen::Vector<T, N> &d_x) {

        const T mass = this->mass_spec.mass;
        const Eigen::Vector3<T> com = this->mass_spec.com;
        const Eigen::Vector3<T> mass_com = mass * com;
        const Eigen::Matrix3<T> inertia = this->mass_spec.inertia;

        Eigen::Matrix<T, 3, 3> m_corner = skew(this->mass_spec.com) * this->frame_g.rotation.transpose();

        local_mass_matrix << inertia, m_corner, m_corner.transpose(), Eigen::Matrix3<T>(
                Eigen::DiagonalMatrix<T, 3>(mass,
                                            mass,
                                            mass));

        this->dynamics_state.mass_matrix =
                this->frame_g.jacobian.transpose() * local_mass_matrix * this->frame_g.jacobian;


        d_link_star << this->frame_g.ang_vel.cross(inertia * this->frame_g.ang_vel),
                this->frame_g.rotation *
                this->frame_g.ang_vel.cross(this->frame_g.ang_vel. cross(mass_com));

        this->dynamics_state.centrifugal =
                this->frame_g.jacobian.transpose() *
                (local_mass_matrix * this->frame_g.dot_jacobian * d_x + d_link_star);
    }

    void differential_dynamics(const Eigen::Vector<T, N> &d_x) {
        _differential_dynamics<T, N>(mass_spec, local_mass_matrix, d_link_star, frame_g, dx_frame_g, dynamics_state, dx_dynamics_state, d_x, Eigen::Matrix<T, N, N>::Zero(), idx);
        _differential_dynamics<T,N>(mass_spec, local_mass_matrix, d_link_star, frame_g, ddx_frame_g, dynamics_state, ddx_dynamics_state, d_x, Eigen::Matrix<T, N, N>::Identity(), idx);
    }



};

#endif //CPP_NE_2_OKC_LINK_H
