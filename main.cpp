#include <iostream>
#include <okc.h>

class MySystem : public OKC<float, 2> {

public:

    using Link = Link<float, 2>;

    Link x = Link {
        .name = "Root",
            .link_type = Revolute,
            .axis = {1, 0, 0},
            .frame_o = {
                    .position = {0, 0, 0},
            },
            .mass_spec = {
                    .mass = 1,
                    .com = {0, 0.5, 0},
                    .inertia = Eigen::Matrix3f::Identity(),
            }
    };

    Link y = Link {
        .parent = &x,
        .link_type = Revolute,
        .axis = {1, 0, 0},
        .frame_o = {
                .position = {0, 1, 0},
        },
        .mass_spec = {
                .mass = 1,
                .com = {0, 0.5, 0},
                .inertia = Eigen::Matrix3f::Identity(),
        }
    };

    Link z = Link {
            .parent = &y,
            .link_type = Revolute,
            .axis = {1, 0, 0},
            .frame_o = {
                    .position = {0, 1, 0},
            },
            .mass_spec = {
                    .mass = 1,
                    .com = {0, 0.5, 0},
                    .inertia = Eigen::Matrix3f::Identity(),
            }
    };

    Link w = Link {
            .parent = &z,
            .link_type = Revolute,
            .axis = {1, 0, 0},
            .frame_o = {
                    .position = {0, 1, 0},
            },
            .mass_spec = {
                    .mass = 1,
                    .com = {0, 0.5, 0},
                    .inertia = Eigen::Matrix3f::Identity(),
            }
    };



    MySystem() {
        set_links({&x, &y});
    }

    void update() override {}

};

#include <chrono>

MySystem my_system;

int main() {


    Eigen::Vector2f x = {1, 2};
    Eigen::Vector2f d_x = {2, 2};
    float h = 0.01;
    my_system.set(x, d_x);

    std::cout << my_system.accel() << std::endl;

    for (size_t i = 0; i < 4; ++i)
        std::cout << my_system.w.dx_frame_g.jacobian[i] << "\n" << std::endl;

   // auto result = my_system.simulate(10, 0.1);

    //std::cout << result["p Root"] << std::endl;

    auto t1 = std::chrono::high_resolution_clock::now();
    for(size_t i = 0; i < 10000; ++i) {
        my_system.step(h);
    }
    auto t2 = std::chrono::high_resolution_clock::now();

    std::cout << my_system.pos() << std::endl;
    std::cout << std::chrono::duration_cast<std::chrono::microseconds>(t2-t1).count() << std::endl;

    return 0;
}
