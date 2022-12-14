//
// Created by danie on 11/17/2022.
//

#ifndef CPP_NE_2_OKC_SYSTEM_H
#define CPP_NE_2_OKC_SYSTEM_H

#include <map>
#include <Eigen/Eigen>
#include "okc_link.h"
#include "okc_types.h"

template <class T, size_t N>
class OKC {

    T _time = 0;
    std::array<Link<T, N>*, N> _links;

    Eigen::Vector<T, N> _x;
    Eigen::Vector<T, N> _dot_x;
    Eigen::Vector<T, N> _ddot_x;
    Eigen::Matrix<T, N, N> _dx_ddot_x;
    Eigen::Matrix<T, N, N> _ddx_ddot_x;

    Dynamics_State<T, N> dynamics_state;

    Differential_Dynamics_State<T, N> dx_dynamics_state;
    Differential_Dynamics_State<T, N> ddx_dynamics_state;


protected:

    explicit OKC(std::array<Link<T, N>*, N> links) : _links(links) {
        for (auto link : _links) {
            assert(link->parent->idx < link->idx);
            link->init();
        }
    }

    OKC() {}

    void set_links(std::array<Link<T, N>*, N> links) {
        for (size_t i = 0; i < N; ++i) {
            Link<T, N>* link = links[i];
            link->idx = i;
            link->init();
            _links[i] = links[i];
        }
    }

    void _update() {

        dynamics_state = Dynamics_State<T, N>::Zero();
        dx_dynamics_state = Differential_Dynamics_State<T, N>::Zero();
        ddx_dynamics_state = Differential_Dynamics_State<T, N>::Zero();

        for (auto link: _links) {
            link->set(_x, _dot_x);
            dynamics_state += link->dynamics_state;
            dx_dynamics_state += link->dx_dynamics_state;
            ddx_dynamics_state += link->ddx_dynamics_state;
        }

        Eigen::Matrix<T, N, N> system_mass_inverse = dynamics_state.mass_matrix.inverse();

        this->_ddot_x = system_mass_inverse *(-dynamics_state.centrifugal);

        for (size_t i = 0; i < N; ++i) {
            _dx_ddot_x.col(i) = system_mass_inverse * (dx_dynamics_state.mass_matrix[i] * (- _ddot_x));
            _ddx_ddot_x.col(i) = system_mass_inverse * (ddx_dynamics_state.mass_matrix[i] * (- _ddot_x));
        }

        update();

    }



public:

    inline const Eigen::Vector<T, N>& pos() const {
        return _x;
    }

    inline const Eigen::Vector<T, N>& vel() const {
        return _dot_x;
    }

    inline const Eigen::Vector<T, N>& accel() const {
        return _ddot_x;
    }

    inline const T time() const {
        return _time;
    }

    void set(const Eigen::Vector<T, N>& joint_positions, const Eigen::Vector<T, N>& joint_velocities) {
        _x = joint_positions;
        _dot_x = joint_velocities;
        _update();
    }

    void step(T step_size) {
        _dot_x += step_size * this->_ddot_x;
        _x += step_size * _dot_x;
        _update();
    }

    std::map<std::string, Eigen::Vector<T, Eigen::Dynamic>> simulate(T runtime, T step_size) {

        size_t len = runtime / step_size;
        Eigen::Matrix<T, N, Eigen::Dynamic> positions = Eigen::Matrix<T, N, Eigen::Dynamic>(N, len);
        Eigen::Matrix<T, N, Eigen::Dynamic> velocities = Eigen::Matrix<T, N, Eigen::Dynamic>(N, len);

        for (size_t i = 0; i < len; ++i) {
            _time += step_size;
            positions.col(i) = this->pos();
            velocities.col(i) = this->vel();
            this->step(step_size);
        }

        std::map<std::string, Eigen::Vector<T, Eigen::Dynamic>> result;

        for (auto link : _links) {
            std::string col_name = link->name.length() ? link->name : std::to_string(link->idx);
            Eigen::Vector<T, Eigen::Dynamic> link_positions = Eigen::Vector<T, Eigen::Dynamic>(len);
            Eigen::Vector<T, Eigen::Dynamic> link_velocities = Eigen::Vector<T, Eigen::Dynamic>(len);
            link_positions = positions.row(link->idx);
            link_velocities = velocities.row(link->idx);
            result["p " + col_name] =  link_positions;
            result["v " + col_name] = link_positions;
        }

        return result;

    }

    virtual void update() {}

};

#endif //CPP_NE_2_OKC_SYSTEM_H
